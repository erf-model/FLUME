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
 * Mirror the cells inside each wall into the ghost layers outside it: layer m
 * (1 is next to the wall) takes the m-th cell inside. `sign[d][n]` multiplies
 * component n at a wall whose normal is direction d.
 *
 * Every fab whose ghost cells reach outside the wall is filled, not only
 * those whose box touches it: a box one cell short of the wall holds the
 * first layer outside it in its second ghost layer.
 *
 * The mirror cell is in the same fab: a valid cell, or, in a box thinner
 * than the ghost layers, a ghost cell that FillBoundary has filled (call it
 * first). In a domain thinner than the ghost layers the mirror of an outer
 * layer is the inner layer beyond the other wall, which is filled before it:
 * the layers are done in order.
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
            const int n_grow = mf.nGrowVect()[d];
            for (int layer = 1; layer <= n_grow; ++layer) {
                for (int side = 0; side < 2; ++side) {
                    const int face = (side == 0) ? domain.smallEnd(d) : domain.bigEnd(d);
                    // The ghost strip outside the wall, and the offset to its mirror.
                    const int ghost = (side == 0) ? face - layer : face + layer;
                    const bool held = (side == 0) ? ghost >= valid.smallEnd(d) - n_grow
                                                  : ghost <= valid.bigEnd(d) + n_grow;
                    if (!held) { continue; }
                    amrex::Box strip = valid;
                    strip.setSmall(d, ghost);
                    strip.setBig(d, ghost);
                    const amrex::IntVect to_mirror = amrex::IntVect::TheDimensionVector(d)
                                                     * ((side == 0) ? 2 * layer - 1 : 1 - 2 * layer);
                    amrex::ParallelFor(strip, ncomp,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                        {
                            const amrex::IntVect m = amrex::IntVect(i, j, k) + to_mirror;
                            a(i, j, k, n) = sgn[n] * a(m, n);
                        });
                }
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
