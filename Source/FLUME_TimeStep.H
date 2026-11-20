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

} // namespace flume

#endif
