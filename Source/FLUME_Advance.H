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
 * hydrostatic reconstruction (face_depth) and the cells' own velocities.
 * `state` and `bed` need one layer of filled ghost cells. `flux[d]` lives on
 * the faces across direction d and holds State::ncomp components in the order
 * of the state: mass, x momentum, y momentum.
 */
void compute_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                     const amrex::MultiFab& bed, amrex::Real gravity);

/**
 * One forward-Euler step of the valid cells:
 *
 *   U <- U - dt/dx (Fx_hi - Fx_lo) - dt/dy (Fy_hi - Fy_lo) + dt S ,
 *
 * where S is the bed-slope term of the hydrostatic reconstruction: in each
 * direction, g/2 (h_hi^2 - h_lo^2) / dx on the momentum of that direction,
 * with h_hi and h_lo the depths the cell presents to its two faces. For water
 * at rest under a level surface the term balances the flux difference.
 * `bed` needs one layer of filled ghost cells; the ghost cells of `state` are
 * not read and are out of date afterwards.
 */
void update_state (amrex::MultiFab& state, const amrex::Array<amrex::MultiFab, 2>& flux,
                   const amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real gravity,
                   amrex::Real dt);

} // namespace flume

#endif
