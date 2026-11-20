#include "FLUME_Boundary.H"
#include "FLUME_State.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <algorithm>

namespace flume {

namespace {

//! One sign per component for each direction of the wall normal.
using WallSigns = amrex::Array<amrex::GpuArray<amrex::Real, State::ncomp>, 2>;

using Pieces = amrex::Array<amrex::Array<std::vector<Boundary::Piece>, 2>, 2>;

/**
 * Fill the ghost layers outside every side that is not periodic, piece by
 * piece. Layer m (1 is next to the side) of a wall piece takes the m-th cell
 * inside, times `sign[d][n]` for component n at a side whose normal is
 * direction d. Every layer of an outflow piece takes the boundary cell.
 *
 * The plain copy keeps still water still to round-off, but where the bed of
 * the boundary cell is lower than that of the cell inside, still water there
 * is weakly unstable: the face to the outside carries the whole discharge of
 * the boundary cell and the face inside, after the hydrostatic
 * reconstruction, only a part of it, and a disturbance grows by a fraction of
 * a percent per step (docs/source/Model.rst, "Sides of the domain").
 *
 * Every fab whose ghost cells reach outside the side is filled, not only
 * those whose box touches it: a box one cell short of the side holds the
 * first layer outside it in its second ghost layer.
 *
 * The cell copied is in the same fab: a valid cell, or, in a box thinner
 * than the ghost layers, a ghost cell that FillBoundary has filled (call it
 * first). In a domain thinner than the ghost layers the mirror of an outer
 * wall layer is the inner layer beyond the other side, which is filled
 * before it: the layers are done in order.
 */
void
fill_sides (amrex::MultiFab& mf, const amrex::Geometry& geom, const Pieces& pieces,
            const WallSigns& sign)
{
    const amrex::Box& domain = geom.Domain();
    const int ncomp = mf.nComp();
    AMREX_ALWAYS_ASSERT(ncomp <= State::ncomp);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& valid = mfi.validbox();
        const auto& a = mf.array(mfi);
        for (int d = 0; d < 2; ++d) {
            if (geom.isPeriodic(d)) { continue; }
            const int t = 1 - d; // along the side
            const auto sgn = sign[d];
            const int n_grow = mf.nGrowVect()[d];
            for (int layer = 1; layer <= n_grow; ++layer) {
                for (int side = 0; side < 2; ++side) {
                    const int face = (side == 0) ? domain.smallEnd(d) : domain.bigEnd(d);
                    // The ghost strip outside the side.
                    const int ghost = (side == 0) ? face - layer : face + layer;
                    const bool held = (side == 0) ? ghost >= valid.smallEnd(d) - n_grow
                                                  : ghost <= valid.bigEnd(d) + n_grow;
                    if (!held) { continue; }
                    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
                    // From a ghost cell to its mirror, and to the boundary cell.
                    const amrex::IntVect to_mirror = e * ((side == 0) ? 2 * layer - 1 : 1 - 2 * layer);
                    const amrex::IntVect to_face = e * (face - ghost);
                    for (const auto& piece : pieces[d][side]) {
                        amrex::Box strip = valid;
                        strip.setSmall(d, ghost);
                        strip.setBig(d, ghost);
                        strip.setSmall(t, std::max(valid.smallEnd(t), piece.lo));
                        strip.setBig(t, std::min(valid.bigEnd(t), piece.hi));
                        if (!strip.ok()) { continue; }
                        if (piece.type == BoundaryType::wall) {
                            amrex::ParallelFor(strip, ncomp,
                                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                                {
                                    const amrex::IntVect m = amrex::IntVect(i, j, k) + to_mirror;
                                    a(i, j, k, n) = sgn[n] * a(m, n);
                                });
                        } else {
                            amrex::ParallelFor(strip, ncomp,
                                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                                {
                                    const amrex::IntVect m = amrex::IntVect(i, j, k) + to_face;
                                    a(i, j, k, n) = a(m, n);
                                });
                        }
                    }
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
        // a periodic side and another side has no image (fill_sides fills
        // it later, or nothing reads it), so it is not touched.
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

const char*
type_name (BoundaryType t)
{
    return (t == BoundaryType::outflow) ? "outflow" : "wall";
}

} // namespace

Boundary::Boundary (const Inputs& in, const amrex::Geometry& geom) : m_geom(geom)
{
    const char* const side_name[2][2] = {{"xlo", "xhi"}, {"ylo", "yhi"}};
    const char* const dir_name[2] = {"x", "y"};
    for (const auto& s : in.bc_segments) {
        amrex::Print() << "FLUME: boundary segment " << s.name << ": " << type_name(s.type)
                       << " on " << side_name[s.direction][s.high] << ", cells " << s.cell_lo
                       << " to " << s.cell_hi << " in " << dir_name[1 - s.direction] << "\n";
    }
    for (int d = 0; d < 2; ++d) {
        if (geom.isPeriodic(d)) { continue; }
        const int t = 1 - d;
        const int n_along = geom.Domain().length(t);
        for (int side = 0; side < 2; ++side) {
            // The segments of this side, in order along it. read_inputs has
            // checked that they lie on the side, hold a cell each and do not
            // overlap.
            std::vector<Piece> seg;
            for (const auto& s : in.bc_segments) {
                if (s.direction == d && s.high == side) {
                    seg.push_back(Piece{s.cell_lo, s.cell_hi, s.type});
                }
            }
            std::sort(seg.begin(), seg.end(),
                      [] (const Piece& a, const Piece& b) { return a.lo < b.lo; });
            // The cells between the segments are walls.
            auto& out = m_pieces[d][side];
            int next = 0;
            for (const auto& p : seg) {
                if (p.lo > next) { out.push_back(Piece{next, p.lo - 1, BoundaryType::wall}); }
                out.push_back(p);
                next = p.hi + 1;
            }
            if (next < n_along) { out.push_back(Piece{next, n_along - 1, BoundaryType::wall}); }
        }
    }
}

void
Boundary::fill_state (amrex::MultiFab& state) const
{
    state.FillBoundary(m_geom.periodicity());
    const amrex::Real p = amrex::Real(1.0);
    const amrex::Real m = amrex::Real(-1.0);
    // Component order: h, hu, hv. A wall across x negates hu; across y, hv.
    const WallSigns sign{{{p, m, p}, {p, p, m}}};
    fill_sides(state, m_geom, m_pieces, sign);
}

void
Boundary::fill_bed (amrex::MultiFab& bed, const amrex::Array<amrex::Real, 2>& continued_slope) const
{
    bed.FillBoundary(m_geom.periodicity());
    for (int d = 0; d < 2; ++d) {
        if (m_geom.isPeriodic(d) && continued_slope[d] != amrex::Real(0.0)) {
            continue_slope(bed, m_geom, d, continued_slope[d]);
        }
    }
    const amrex::Real p = amrex::Real(1.0);
    const WallSigns sign{{{p, p, p}, {p, p, p}}};
    fill_sides(bed, m_geom, m_pieces, sign);
}

} // namespace flume
