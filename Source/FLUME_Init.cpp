#include "FLUME_Init.H"
#include "FLUME_State.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Math.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
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
