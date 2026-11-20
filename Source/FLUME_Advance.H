#ifndef FLUME_ADVANCE_H_
#define FLUME_ADVANCE_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * HLLC flux on every face of every box, in x (flux[0]) and y (flux[1]).
 *
 * Each face takes the two cells next to it, with the depths of the
 * hydrostatic reconstruction (face_depths) and the cells' own velocities.
 * A cell whose depth is at or below `dry_depth` has no velocity, and presents
 * no depth to its faces: nothing flows out of it, and water can flow in.
 * `state` and `bed` need one layer of filled ghost cells. `flux[d]` lives on
 * the faces across direction d and holds State::ncomp components in the order
 * of the state: mass, x momentum, y momentum.
 */
void compute_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                     const amrex::MultiFab& bed, amrex::Real gravity, amrex::Real dry_depth);

/**
 * As compute_fluxes, second order in space: each face takes the two cells
 * next to it as their piecewise-linear reconstruction presents them
 * (FLUME_Reconstruct.H): depth, bed and velocities on the side of the face.
 * The hydrostatic reconstruction then acts on these in place of the cell
 * values. `state` and `bed` need two layers of filled ghost cells.
 */
void compute_fluxes_second_order (amrex::Array<amrex::MultiFab, 2>& flux,
                                  const amrex::MultiFab& state, const amrex::MultiFab& bed,
                                  amrex::Real gravity, amrex::Real dry_depth);

/**
 * Scale the fluxes so that no cell gives away more water in the step `dt`
 * than it holds.
 *
 * For each valid cell the mass that leaves through its four faces in the step
 * is added up. Where that is more than the depth h, the cell's factor is
 * h / (mass leaving), below 1; elsewhere it is 1. The three fluxes of a face
 * are then multiplied by the factor of the cell the mass flux leaves. A face
 * has one flux for its two cells, so the step still conserves mass, and the
 * new depth of every cell is >= 0 up to rounding. A face with no mass flux is
 * not scaled, so water at rest is not touched.
 *
 * `factor` is work space: one component and one layer of ghost cells on the
 * boxes of `state`. The ghost cells of `state` are not read.
 */
void limit_outflow (amrex::Array<amrex::MultiFab, 2>& flux, amrex::MultiFab& factor,
                    const amrex::MultiFab& state, const amrex::Geometry& geom, amrex::Real dt);

/**
 * Bed-slope term of the hydrostatic reconstruction in every valid cell:
 * g/2 (h_hi^2 - h_lo^2) in x (component 0) and in y (component 1), with h_hi
 * and h_lo the depths the cell presents to its two faces (face_depths, as in
 * compute_fluxes). For water at rest under a level surface the term balances
 * the difference of the momentum fluxes.
 *
 * `source` holds two components on the boxes of `state` and needs no ghost
 * cells. `state` and `bed` need one layer of filled ghost cells.
 */
void compute_bed_source (amrex::MultiFab& source, const amrex::MultiFab& state,
                         const amrex::MultiFab& bed, amrex::Real gravity, amrex::Real dry_depth);

/**
 * As compute_bed_source, for compute_fluxes_second_order:
 *
 *   g/2 (h_hi^2 - h_lo^2) - g h d_eta ,
 *
 * with h_hi and h_lo the depths the reconstructed sides of the cell present
 * to its two faces, and d_eta the limited difference of the surface elevation
 * over the cell. The second part is the slope of the bed inside the cell.
 * Under a level surface d_eta is zero (to rounding) and the first part
 * balances the fluxes as at first order.
 *
 * `state` and `bed` need two layers of filled ghost cells.
 */
void compute_bed_source_second_order (amrex::MultiFab& source, const amrex::MultiFab& state,
                                      const amrex::MultiFab& bed, amrex::Real gravity,
                                      amrex::Real dry_depth);

/**
 * One forward-Euler step of the valid cells:
 *
 *   U <- U - dt/dx (Fx_hi - Fx_lo) - dt/dy (Fy_hi - Fy_lo) + dt S ,
 *
 * where S is the bed-slope term: `source` (compute_bed_source, from the state
 * before the step) over the cell size, on the momentum of each direction.
 *
 * A new depth below zero by rounding only (limit_outflow leaves up to a few
 * machine epsilons of the old depth) is set to zero. A cell whose new depth
 * is at or below `dry_depth` keeps its water and loses its momentum.
 * The ghost cells of `state` are not read and are out of date afterwards.
 */
void update_state (amrex::MultiFab& state, const amrex::Array<amrex::MultiFab, 2>& flux,
                   const amrex::MultiFab& source, const amrex::Geometry& geom,
                   amrex::Real dry_depth, amrex::Real dt);

/**
 * Last part of the two-stage step: the valid cells of `state` (the result of
 * two forward-Euler steps) become the mean of `old` (the state before them)
 * and `state`. A cell whose depth is then at or below `dry_depth` loses its
 * momentum. The mean of two depths >= 0 is >= 0.
 */
void average_with (amrex::MultiFab& state, const amrex::MultiFab& old, amrex::Real dry_depth);

} // namespace flume

#endif
