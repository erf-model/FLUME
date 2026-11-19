#ifndef FLUME_BOUNDARY_H_
#define FLUME_BOUNDARY_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace flume {

/**
 * Fill the one layer of ghost cells of the state (h, hu, hv).
 *
 * Between boxes and across a periodic side of the domain the ghost cell takes
 * the value of the cell it covers. At a side of the domain that is not
 * periodic there is a wall: the ghost cell mirrors the cell inside it, with
 * the discharge normal to the wall negated, so no water crosses the wall.
 * No kernel reads a corner ghost cell (one that is outside the box in both x
 * and y). Those outside a wall are not filled.
 */
void fill_state_ghost_cells (amrex::MultiFab& state, const amrex::Geometry& geom);

/**
 * Fill the one layer of ghost cells of the bed elevation: as the state, and
 * at a wall the ghost cell takes the bed of the cell inside it.
 */
void fill_bed_ghost_cells (amrex::MultiFab& bed, const amrex::Geometry& geom);

} // namespace flume

#endif
