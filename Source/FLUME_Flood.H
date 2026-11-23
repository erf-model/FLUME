#ifndef FLUME_FLOOD_H_
#define FLUME_FLOOD_H_

#include "FLUME_Boundary.H"
#include "FLUME_FloodMaps.H"
#include "FLUME_Forcing.H"
#include "FLUME_Gauges.H"
#include "FLUME_Inputs.H"
#include "FLUME_Source.H"
#include "FLUME_WaterBalance.H"

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>
#include <vector>

namespace flume {

/**
 * The flood model on one level: grid, state, bed, time loop and plotfiles.
 *
 * The grid is two-dimensional: nx by ny cells and one cell in z. The z extent
 * of the domain is 0 to 1 m and carries no meaning.
 *
 * The state (FLUME_State.H) and the bed elevation carry ghost cells in x and
 * y (FLUME_Boundary.H): one layer at first order, two at second order. The
 * bed does not change, so its ghost cells are filled once; those of the
 * state are filled at the start of each stage. A stage is HLLC fluxes on the
 * faces from the hydrostatic reconstruction, scaled where a cell would give
 * away more water than it holds, then a forward-Euler update
 * (FLUME_Advance.H). With flume.order = 1 a step is one stage from the cell
 * values. With flume.order = 2 the faces take a piecewise-linear
 * reconstruction of the cells, and a step is the mean of the state and of
 * two stages from it. The step ends with the source terms
 * (FLUME_Source.H), each over the whole step, in order: rain from the
 * forcing provider (FLUME_Forcing.H, FLUME_Rain.H), then infiltration
 * (FLUME_Infiltration.H), then Manning friction, implicit in the discharge
 * (FLUME_Friction.H).
 *
 * After each step the flood maps take the new state (FLUME_FloodMaps.H),
 * and every flume.log_int steps, and at the last step, the water balance
 * (FLUME_WaterBalance.H) and the gauges (FLUME_Gauges.H) write a line.
 */
class Flood
{
public:
    //! Build the grid and fill the bed and the initial state. Aborts if any
    //! initial value is not finite.
    explicit Flood (const Inputs& in);

    //! Write the initial plotfile, run the time loop, write the last plotfile.
    void evolve ();

    /**
     * The next time step, s: cfl / max_wave_rate(), or with rain falling
     * during it rain_limited_step() (FLUME_TimeStep.H). A step that reaches
     * flume.stop_time, or falls short of it by round-off only (4 machine
     * epsilons of the stop time), is set to end on it. With no wave anywhere
     * (a dry domain) `dry` is set: the step is the rain-limited one if rain
     * falls before flume.stop_time (or ever, without a stop time), and
     * otherwise goes straight to flume.stop_time; with neither rain nor a
     * stop time it aborts. Aborts as check_state() does, and on a step that
     * is zero, not finite or too small to advance the time. `lands` is set
     * when the step ends at flume.stop_time.
     */
    double compute_dt (bool& lands, bool& dry) const;

    //! Abort if a depth is negative or not finite, or a wave speed is not
    //! finite, in any cell. Returns the largest wave rate (FLUME_TimeStep.H).
    amrex::Real check_state () const;

    //! Write <plot_file><step>: depth, discharge_x, discharge_y,
    //! bed_elevation, surface_elevation (bed + depth), then the fields of the
    //! sources (Source::plot_fields; infiltrated_depth with infiltration),
    //! then the flood maps (with flume.flood_maps = 1).
    void write_plotfile () const;

    [[nodiscard]] double time () const { return m_time; }
    [[nodiscard]] int step () const { return m_step; }
    [[nodiscard]] const amrex::MultiFab& state () const { return m_state; }
    [[nodiscard]] const amrex::MultiFab& bed () const { return m_bed; }
    [[nodiscard]] const amrex::Geometry& geom () const { return m_geom; }
    [[nodiscard]] const WaterBalance& balance () const { return *m_balance; }

private:
    //! One forward-Euler stage of `dt` on m_state, at the order of the run.
    //! The water balance counts the outflow of the stage over `dt_weight`, s:
    //! the step, or half of it for each stage of the two-stage step.
    void advance_stage (amrex::Real dt, double dt_weight);

    //! The water-balance line and the gauge lines of the current step.
    void write_logs ();

    Inputs m_in;
    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::MultiFab m_state; //!< State::ncomp components, ghost layers in x and y
    amrex::MultiFab m_bed;   //!< bed elevation, m, ghost layers in x and y
    //! The state before the step (flume.order = 2 only), no ghost cells
    amrex::MultiFab m_state_old;
    //! Fluxes on the faces across x and across y, State::ncomp components
    amrex::Array<amrex::MultiFab, 2> m_flux;
    //! Work space of limit_outflow: one component, one ghost layer in x and y
    amrex::MultiFab m_outflow_factor;
    //! Bed-slope term in x and in y (compute_bed_source), no ghost cells
    amrex::MultiFab m_bed_source;
    //! Building height of each cell, m (flume.buildings.file), one
    //! component, no ghost cells; not defined without buildings
    amrex::MultiFab m_buildings;
    //! The sides that are not periodic (flume.periodic, flume.bc.*)
    Boundary m_boundary;
    //! What drives the flood from outside (flume.rain.*), or null; declared
    //! before m_sources, which may hold a reference to it
    std::unique_ptr<Forcing> m_forcing;
    //! The source terms, applied in order at the end of each step
    std::vector<std::unique_ptr<Source>> m_sources;
    //! The water-balance log (flume.log_int)
    std::unique_ptr<WaterBalance> m_balance;
    //! The flood maps (flume.flood_maps = 1), or null
    std::unique_ptr<FloodMaps> m_maps;
    //! The gauge points (flume.gauges.*)
    std::unique_ptr<Gauges> m_gauges;
    //! Time, s. In double precision in every build: a single-precision sum of
    //! steps stops growing after about 2^24 steps.
    double m_time = 0.0;
    int m_step = 0;
};

} // namespace flume

#endif
