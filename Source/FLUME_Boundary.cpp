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

/**
 * Add slope * (distance to the image cell) to every ghost cell of `bed`
 * beyond a periodic side in direction d. The image is the cell inside the
 * domain that FillBoundary copied: p whole periods of the domain away, with
 * p = floor((i - lo) / n) for the index i along d, lo the first index of the
 * domain and n its number of cells. Cells inside the domain have p = 0.
 */
void
continue_slope (amrex::MultiFab& bed, const amrex::Geometry& geom, int d, amrex::Real slope)
{
    const int lo = geom.Domain().smallEnd(d);
    const int n = geom.Domain().length(d);
    const amrex::Real period = slope * geom.CellSize(d) * static_cast<amrex::Real>(n);
    for (amrex::MFIter mfi(bed); mfi.isValid(); ++mfi) {
        // Only cells beyond the periodic sides: a corner ghost cell beyond
        // a periodic side and a wall has no image (fill_walls fills it
        // later, or nothing reads it), so it is not touched.
        const amrex::Box gbx = mfi.fabbox() & geom.growPeriodicDomain(bed.nGrowVect());
        const auto& a = bed.array(mfi);
        amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int c = amrex::IntVect(i, j, k)[d];
                const int p = (c >= lo) ? (c - lo) / n : -((lo - c + n - 1) / n);
                if (p != 0) { a(i, j, k) += static_cast<amrex::Real>(p) * period; }
            });
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
fill_bed_ghost_cells (amrex::MultiFab& bed, const amrex::Geometry& geom,
                      const amrex::Array<amrex::Real, 2>& continued_slope)
{
    bed.FillBoundary(geom.periodicity());
    for (int d = 0; d < 2; ++d) {
        if (geom.isPeriodic(d) && continued_slope[d] != amrex::Real(0.0)) {
            continue_slope(bed, geom, d, continued_slope[d]);
        }
    }
    const amrex::Real p = amrex::Real(1.0);
    const WallSigns sign{{{p, p, p}, {p, p, p}}};
    fill_walls(bed, geom, sign);
}

} // namespace flume
