#ifndef FLUME_BOUNDARY_H_
#define FLUME_BOUNDARY_H_

#include "FLUME_Inputs.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace flume {

/**
 * The sides of the domain and how their ghost cells are filled.
 *
 * Between boxes and across a periodic side of the domain a ghost cell takes
 * the value of the cell it covers. A side that is not periodic is made of
 * segments (flume.bc.*), each a run of boundary cells with one type; the
 * cells no segment covers are walls. So without flume.bc.segments every side
 * that is not periodic is a wall.
 *
 *   wall      the ghost cells mirror the cells inside (the second layer the
 *             second cell), with the discharge normal to the side negated:
 *             no water crosses it.
 *   outflow   every ghost layer takes the boundary cell next to it (zero
 *             gradient, also of the bed): water leaves, or enters, as the
 *             flow at the side carries it, and a wave leaves with little
 *             reflection. Still water over a boundary cell lower than the
 *             cell inside is weakly unstable (FLUME_Boundary.cpp). At second
 *             order the limited slope of the first ghost layer is zero
 *             (minmod of a zero difference), so what the second layer holds
 *             changes no result; it is a copy so that it is finite.
 *   discharge the ghost cells are filled as for outflow; the fluxes on the
 *             faces of the side then come from the inflow hydrograph of the
 *             segment, not from the ghost cells (FLUME_Inflow.H).
 *
 * No kernel reads a corner ghost cell (one that is outside the box in both x
 * and y). Those outside a side that is not periodic are not filled.
 */
class Boundary
{
public:
    //! A run of boundary cells of one side, cells lo to hi along the side.
    struct Piece
    {
        int lo = 0;
        int hi = -1;
        BoundaryType type = BoundaryType::wall;
    };

    //! Resolve the segments of `in` on the grid of `geom` (flume.bc.*; the
    //! cells of each segment were checked by read_inputs) and print them.
    Boundary (const Inputs& in, const amrex::Geometry& geom);

    //! Fill the ghost cells of the state (h, hu, hv): one layer, or two for
    //! the second-order scheme.
    void fill_state (amrex::MultiFab& state) const;

    /**
     * Fill the ghost cells of the bed elevation: as the state, with no sign
     * change at a wall.
     *
     * Across a periodic side in direction d with continued_slope[d] not 0, a
     * ghost cell then adds continued_slope[d] times its distance from the cell
     * it covers (flume.bed.continue_slope): the plane of the bed goes on beyond
     * the side. Pass 0 0 for a bed that repeats.
     */
    void fill_bed (amrex::MultiFab& bed, const amrex::Array<amrex::Real, 2>& continued_slope) const;

    //! The pieces of side `side` (0 low, 1 high) across direction d, in order
    //! along the side, covering every boundary cell. Empty for a periodic
    //! direction.
    [[nodiscard]] const std::vector<Piece>& pieces (int d, int side) const
    {
        return m_pieces[d][side];
    }

private:
    amrex::Geometry m_geom;
    amrex::Array<amrex::Array<std::vector<Piece>, 2>, 2> m_pieces;
};

} // namespace flume

#endif
