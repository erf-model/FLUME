#ifndef FLUME_INIT_H_
#define FLUME_INIT_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_INT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * Bed elevation as a plane:
 * z_b = elevation + slope_x (x - prob_lo_x) + slope_y (y - prob_lo_y),
 * with x and y at the cell centre. Fills the valid cells.
 */
void fill_plane_bed (amrex::MultiFab& bed, const amrex::Geometry& geom,
                     amrex::Real elevation, const amrex::Array<amrex::Real, 2>& slope);

/**
 * State from a level water surface at elevation `surface` and a uniform
 * velocity: h = max(surface - z_b, 0), hu = h u, hv = h v. A cell whose bed is
 * at or above the surface is dry. Fills the valid cells.
 */
void fill_state_from_surface (amrex::MultiFab& state, const amrex::MultiFab& bed,
                              amrex::Real surface, const amrex::Array<amrex::Real, 2>& velocity);

//! Number of valid cells, over all ranks, in which any component is not finite.
amrex::Long count_not_finite (const amrex::MultiFab& mf);

} // namespace flume

#endif
