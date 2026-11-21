#ifndef FLUME_TIMESTEP_H_
#define FLUME_TIMESTEP_H_

#include <AMReX_Geometry.H>
#include <AMReX_INT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

//! Result of max_wave_rate().
struct WaveRate
{
    //! Largest (|u| + c) / dx + (|v| + c) / dy over the wet cells, 1/s
    amrex::Real rate = amrex::Real(0.0);
    //! Cells with a negative or not-finite depth, or a not-finite rate
    amrex::Long bad_cells = 0;
};

/**
 * The largest rate at which gravity waves cross a cell, over all ranks:
 * max over cells of (|u| + c) / dx + (|v| + c) / dy, with u = hu / h,
 * v = hv / h and c = sqrt(g h). The two directions are added, which is the
 * limit an unsplit two-dimensional scheme needs. A cell with h <= dry_depth
 * is dry: it has no wave and counts as zero, so a thin film does not set the
 * time step. The CFL time step is cfl / rate.
 */
WaveRate max_wave_rate (const amrex::MultiFab& state, const amrex::Geometry& geom,
                        amrex::Real gravity, amrex::Real dry_depth);

/**
 * The time step with rain falling, s: the largest dt with
 *
 *   dt (R + S sqrt(r_max dt)) <= cfl,   S = sqrt(g) (1/dx + 1/dy),
 *   R = max(rate, sqrt(g dry_depth) (1/dx + 1/dy)),
 *
 * with `rate` the largest wave rate of the state (max_wave_rate) and `r_max`
 * (m/s, > 0) the largest rain rate over the step. The rain of the step adds
 * at most r_max dt to a depth, which lowers |u| and |v| and raises c by at
 * most sqrt(g r_max dt); R covers the cells at or below dry_depth, which
 * count as zero in `rate`. So the state at the end of the step also meets
 * the CFL condition with this step. The stages use the state at the start of
 * the step, so the stability of this step needs only the usual CFL step; the
 * bound gives a dry domain under rain a step that is finite, and keeps the
 * rain of one step small next to the waves it makes on a thin film. With
 * y = sqrt(dt) the bound is the root of the cubic
 * S sqrt(r_max) y^3 + R y^2 - cfl = 0, found by Newton's method from the
 * smaller of the roots of its two terms alone, which is above it (the cubic
 * is increasing and convex for y > 0, so the iterates fall monotonically to
 * the root). A run aborts if the iteration ends away from the root.
 */
double rain_limited_step (double rate, double r_max, double gravity, double dry_depth,
                          double dxi, double dyi, double cfl);

} // namespace flume

#endif
