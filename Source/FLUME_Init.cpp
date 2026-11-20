#include "FLUME_Init.H"
#include "FLUME_State.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Math.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include <cstdint>
#include <AMReX_Tuple.H>
#include <AMReX_TypeList.H>

namespace flume {

void
fill_plane_bed (amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real elevation,
                const amrex::Array<amrex::Real, 2>& slope)
{
    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);
    const amrex::Real sx = slope[0];
    const amrex::Real sy = slope[1];
    const int ilo = geom.Domain().smallEnd(0);
    const int jlo = geom.Domain().smallEnd(1);
    const auto& zb = bed.arrays();
    amrex::ParallelFor(bed,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            // Distance of the cell centre from the lower corner of the domain.
            const amrex::Real x = (static_cast<amrex::Real>(i - ilo) + amrex::Real(0.5)) * dx;
            const amrex::Real y = (static_cast<amrex::Real>(j - jlo) + amrex::Real(0.5)) * dy;
            zb[box_no](i, j, k) = elevation + sx * x + sy * y;
        });
    amrex::Gpu::streamSynchronize();
}

namespace {

//! Distance of a cell centre from the centre of the domain along one
//! direction, in cell sizes: index i of n cells that start at ilo.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cells_from_centre (int i, int ilo, int n) noexcept
{
    // (i - ilo + 1/2) - n/2, as an exact half-integer.
    return amrex::Real(0.5) * static_cast<amrex::Real>(2 * (i - ilo) + 1 - n);
}

} // namespace

void
fill_state_from_surface (amrex::MultiFab& state, const amrex::MultiFab& bed, amrex::Real surface,
                         const amrex::Array<amrex::Real, 2>& velocity)
{
    const amrex::Real u = velocity[0];
    const amrex::Real v = velocity[1];
    const auto& zb = bed.const_arrays();
    const auto& s = state.arrays();
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real h = amrex::max(surface - zb[box_no](i, j, k), amrex::Real(0.0));
            s[box_no](i, j, k, State::h) = h;
            s[box_no](i, j, k, State::hu) = h * u;
            s[box_no](i, j, k, State::hv) = h * v;
        });
    amrex::Gpu::streamSynchronize();
}

void
fill_state_from_curved_surface (amrex::MultiFab& state, const amrex::MultiFab& bed,
                                const amrex::Geometry& geom, amrex::Real surface,
                                const amrex::Array<amrex::Real, 2>& slope,
                                amrex::Real curvature,
                                const amrex::Array<amrex::Real, 2>& velocity)
{
    const amrex::Real u = velocity[0];
    const amrex::Real v = velocity[1];
    const amrex::Real sx = slope[0];
    const amrex::Real sy = slope[1];
    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);
    const int ilo = geom.Domain().smallEnd(0);
    const int jlo = geom.Domain().smallEnd(1);
    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);
    const auto& zb = bed.const_arrays();
    const auto& s = state.arrays();
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real x = cells_from_centre(i, ilo, nx) * dx;
            const amrex::Real y = cells_from_centre(j, jlo, ny) * dy;
            const amrex::Real eta = surface + (sx * x + sy * y) + curvature * (x * x + y * y);
            const amrex::Real h = amrex::max(eta - zb[box_no](i, j, k), amrex::Real(0.0));
            s[box_no](i, j, k, State::h) = h;
            s[box_no](i, j, k, State::hu) = h * u;
            s[box_no](i, j, k, State::hv) = h * v;
        });
    amrex::Gpu::streamSynchronize();
}

void
add_bed_bowl (amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real depth,
              amrex::Real radius)
{
    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);
    const int ilo = geom.Domain().smallEnd(0);
    const int jlo = geom.Domain().smallEnd(1);
    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);
    const amrex::Real inv_r2 = amrex::Real(1.0) / (radius * radius);
    const auto& zb = bed.arrays();
    amrex::ParallelFor(bed,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real x = cells_from_centre(i, ilo, nx) * dx;
            const amrex::Real y = cells_from_centre(j, jlo, ny) * dy;
            zb[box_no](i, j, k) -= depth * (amrex::Real(1.0) - (x * x + y * y) * inv_r2);
        });
    amrex::Gpu::streamSynchronize();
}

void
add_bed_noise (amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real amplitude)
{
    const int ilo = geom.Domain().smallEnd(0);
    const int jlo = geom.Domain().smallEnd(1);
    const auto& zb = bed.arrays();
    amrex::ParallelFor(bed,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            // 32-bit unsigned arithmetic, which wraps.
            std::uint32_t n = (static_cast<std::uint32_t>(i - ilo) * std::uint32_t(73856093))
                              ^ (static_cast<std::uint32_t>(j - jlo) * std::uint32_t(19349663));
            n ^= n >> 13;
            n *= std::uint32_t(0x5bd1e995);
            n ^= n >> 15;
            // The low 11 bits, 0 to 2047, as a number in [-1, 1).
            const amrex::Real r = (static_cast<amrex::Real>(n & std::uint32_t(2047)) - amrex::Real(1024.0))
                                  / amrex::Real(1024.0);
            zb[box_no](i, j, k) += amplitude * r;
        });
    amrex::Gpu::streamSynchronize();
}

void
fill_state_dam (amrex::MultiFab& state, const amrex::MultiFab& bed, const amrex::Geometry& geom,
                int dir, amrex::Real position, amrex::Real surface_lo, amrex::Real surface_hi,
                const amrex::Array<amrex::Real, 2>& velocity_lo,
                const amrex::Array<amrex::Real, 2>& velocity_hi)
{
    const amrex::Real u_lo = velocity_lo[0];
    const amrex::Real v_lo = velocity_lo[1];
    const amrex::Real u_hi = velocity_hi[0];
    const amrex::Real v_hi = velocity_hi[1];
    const amrex::Real x_lo = geom.ProbLo(dir);
    const amrex::Real dx = geom.CellSize(dir);
    const int i_lo = geom.Domain().smallEnd(dir);
    const auto& zb = bed.const_arrays();
    const auto& s = state.arrays();
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::IntVect iv(i, j, k);
            const amrex::Real x = x_lo + (static_cast<amrex::Real>(iv[dir] - i_lo)
                                          + amrex::Real(0.5)) * dx;
            const bool low = x < position;
            const amrex::Real surface = low ? surface_lo : surface_hi;
            const amrex::Real h = amrex::max(surface - zb[box_no](i, j, k), amrex::Real(0.0));
            s[box_no](i, j, k, State::h) = h;
            s[box_no](i, j, k, State::hu) = h * (low ? u_lo : u_hi);
            s[box_no](i, j, k, State::hv) = h * (low ? v_lo : v_hi);
        });
    amrex::Gpu::streamSynchronize();
}

amrex::Long
count_not_finite (const amrex::MultiFab& mf)
{
    const int ncomp = mf.nComp();
    const auto& a = mf.const_arrays();
    amrex::Long n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                     amrex::TypeList<amrex::Long>{}, mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            amrex::Long bad = 0;
            for (int c = 0; c < ncomp; ++c) {
                if (!amrex::Math::isfinite(a[box_no](i, j, k, c))) { bad = 1; }
            }
            return {bad};
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace flume
