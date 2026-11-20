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

/**
 * As fill_state_from_surface, under a tilted and curved surface:
 * eta = surface + slope_x X + slope_y Y + curvature (X^2 + Y^2),
 * with X, Y the distance of the cell centre from the centre of the domain.
 */
void fill_state_from_curved_surface (amrex::MultiFab& state, const amrex::MultiFab& bed,
                                     const amrex::Geometry& geom, amrex::Real surface,
                                     const amrex::Array<amrex::Real, 2>& slope,
                                     amrex::Real curvature,
                                     const amrex::Array<amrex::Real, 2>& velocity);

/**
 * Add a paraboloid bowl to the valid cells: the bed moves by
 * -depth (1 - r^2 / radius^2), with r the distance of the cell centre from
 * the centre of the domain. A negative `depth` makes a mound.
 */
void add_bed_bowl (amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real depth,
                   amrex::Real radius);

/**
 * Add a rough-bed perturbation to the valid cells: cell (i, j), counted from
 * the lower corner of the domain, moves by amplitude * r(i, j), where r is a
 * fixed pseudo-random number in [-1, 1) (a multiple of 1/1024) from an
 * integer hash of the two indices. The bed does not depend on the boxes or
 * the number of ranks. Tests/check_plotfile.py holds the same hash.
 */
void add_bed_noise (amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real amplitude);

/**
 * State on the two sides of a dam: as fill_state_from_surface, with the
 * surface at `surface_lo` and the velocity `velocity_lo` in the cells whose
 * centre along direction `dir` (0 for x, 1 for y) is below `position`, and
 * `surface_hi` and `velocity_hi` in the others. Fills the valid cells.
 */
void fill_state_dam (amrex::MultiFab& state, const amrex::MultiFab& bed,
                     const amrex::Geometry& geom, int dir, amrex::Real position,
                     amrex::Real surface_lo, amrex::Real surface_hi,
                     const amrex::Array<amrex::Real, 2>& velocity_lo,
                     const amrex::Array<amrex::Real, 2>& velocity_hi);

//! Number of valid cells, over all ranks, in which any component is not finite.
amrex::Long count_not_finite (const amrex::MultiFab& mf);

} // namespace flume

#endif
