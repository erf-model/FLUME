#ifndef FLUME_BOUNDARY_H_
#define FLUME_BOUNDARY_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace flume {

/**
 * Fill the ghost cells of the state (h, hu, hv): one layer, or two for the
 * second-order scheme.
 *
 * Between boxes and across a periodic side of the domain the ghost cell takes
 * the value of the cell it covers. At a side of the domain that is not
 * periodic there is a wall: the ghost cells mirror the cells inside it (the
 * second layer the second cell), with the discharge normal to the wall
 * negated, so no water crosses the wall.
 * No kernel reads a corner ghost cell (one that is outside the box in both x
 * and y). Those outside a wall are not filled.
 */
void fill_state_ghost_cells (amrex::MultiFab& state, const amrex::Geometry& geom);

/**
 * Fill the ghost cells of the bed elevation: as the state, and at a wall
 * each ghost cell takes the bed of the cell it mirrors.
 *
 * Across a periodic side in direction d with continued_slope[d] not 0, a
 * ghost cell then adds continued_slope[d] times its distance from the cell
 * it covers (flume.bed.continue_slope): the plane of the bed goes on beyond
 * the side. Pass 0 0 for a bed that repeats.
 */
void fill_bed_ghost_cells (amrex::MultiFab& bed, const amrex::Geometry& geom,
                           const amrex::Array<amrex::Real, 2>& continued_slope);

} // namespace flume

#endif
