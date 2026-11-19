#include "FLUME_Boundary.H"
#include "FLUME_State.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>

namespace flume {

namespace {

//! One sign per component for each direction of the wall normal.
using WallSigns = amrex::Array<amrex::GpuArray<amrex::Real, State::ncomp>, 2>;

/**
 * Mirror the first layer of cells inside each wall into the ghost layer
 * outside it. `sign[d][n]` multiplies component n at a wall whose normal is
 * direction d. The mirror cell of a ghost cell next to the domain is the
 * first cell inside, which is in the same box.
 */
void
fill_walls (amrex::MultiFab& mf, const amrex::Geometry& geom, const WallSigns& sign)
{
    const amrex::Box& domain = geom.Domain();
    const int ncomp = mf.nComp();
    AMREX_ALWAYS_ASSERT(ncomp <= State::ncomp);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& valid = mfi.validbox();
        const auto& a = mf.array(mfi);
        for (int d = 0; d < 2; ++d) {
            if (geom.isPeriodic(d)) { continue; }
            const auto sgn = sign[d];
            for (int side = 0; side < 2; ++side) {
                const int face = (side == 0) ? domain.smallEnd(d) : domain.bigEnd(d);
                const int edge = (side == 0) ? valid.smallEnd(d) : valid.bigEnd(d);
                if (edge != face) { continue; }
                // The ghost strip outside the wall, and the offset to its mirror.
                amrex::Box strip = valid;
                const int ghost = (side == 0) ? face - 1 : face + 1;
                strip.setSmall(d, ghost);
                strip.setBig(d, ghost);
                const amrex::IntVect to_mirror = amrex::IntVect::TheDimensionVector(d)
                                                 * ((side == 0) ? 1 : -1);
                amrex::ParallelFor(strip, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        const amrex::IntVect m = amrex::IntVect(i, j, k) + to_mirror;
                        a(i, j, k, n) = sgn[n] * a(m, n);
                    });
            }
        }
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

void
fill_state_ghost_cells (amrex::MultiFab& state, const amrex::Geometry& geom)
{
    state.FillBoundary(geom.periodicity());
    const amrex::Real p = amrex::Real(1.0);
    const amrex::Real m = amrex::Real(-1.0);
    // Component order: h, hu, hv. A wall across x negates hu; across y, hv.
    const WallSigns sign{{{p, m, p}, {p, p, m}}};
    fill_walls(state, geom, sign);
}

void
fill_bed_ghost_cells (amrex::MultiFab& bed, const amrex::Geometry& geom)
{
    bed.FillBoundary(geom.periodicity());
    const amrex::Real p = amrex::Real(1.0);
    const WallSigns sign{{{p, p, p}, {p, p, p}}};
    fill_walls(bed, geom, sign);
}

} // namespace flume
