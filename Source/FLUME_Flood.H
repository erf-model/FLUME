#ifndef FLUME_FLOOD_H_
#define FLUME_FLOOD_H_

#include "FLUME_Boundary.H"
#include "FLUME_FloodMaps.H"
#include "FLUME_Forcing.H"
#include "FLUME_Gauges.H"
#include "FLUME_Inflow.H"
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
#include <optional>
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
 * (FLUME_Infiltration.H), then the storm-drain inlets (FLUME_Drains.H),
 * then Manning friction, implicit in the discharge
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
    //! initial value is not finite. With flume.restart, then continue from
    //! that checkpoint (restart()).
    explicit Flood (const Inputs& in);

    /**
     * Build on the boxes and owners of a host model (ERF,
     * FLUME_ERFFloodEngine.cpp) instead of chopping the domain by
     * flume.max_grid_size. The boxes must cover the domain of `in` once, one
     * cell thick (aborts otherwise). No restart: the host calls restart().
     */
    Flood (const Inputs& in, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! Write the initial plotfile (not on a restart), run the time loop,
    //! write the last plotfile, and the checkpoints with flume.check_int.
    //! The stand-alone run: start_maps(), write_plotfile(), start_logs(),
    //! step_to(flume.stop_time) until the end, then finish(true).
    void evolve ();

    //! The flood maps take the initial state (with flume.flood_maps = 1).
    void start_maps ();

    //! The water balance takes the initial state, and the gauge files are
    //! started with the line of the current step (with flume.log_int > 0).
    void start_logs ();

    /**
     * Take one step, set by compute_dt() with `target` as the time to land
     * on (flume.stop_time in a stand-alone run, the end of the host's step
     * inside a host model). The maps take the new state; the logs, plotfile
     * and checkpoint are written every flume.log_int, plot_int and check_int
     * steps.
     */
    void step_to (const std::optional<double>& target);

    //! The end of the run: the last plotfile (with `write_plot`), log lines
    //! and checkpoint where the last step did not write them, and the
    //! water-balance summary.
    void finish (bool write_plot);

    //! Set the time, s, to the host model's (finite, >= 0; aborts otherwise).
    void set_time (double t);

    //! Whether the run continues from a checkpoint (restart()).
    [[nodiscard]] bool restarted () const { return m_restart_step >= 0; }

    /**
     * Continue from the checkpoint `dir`: the time, step, state, sources,
     * balance and maps; the logs are cut back to its step and appended to.
     * Aborts, naming flume.restart, when the checkpoint does not match the
     * inputs: another precision, grid, periodicity, bed, list of sources
     * or of outflow and discharge segments, or no maps with flume.flood_maps = 1.
     */
    void restart (const std::string& dir);

    /**
     * Write <check_file><step>: the time and step, the state, the bed, what
     * each source carries from step to step, the water-balance totals and
     * the flood maps. Collective.
     */
    void write_checkpoint () const;

    //! write_checkpoint() into the directory `name`.
    void write_checkpoint (const std::string& name) const;

    /**
     * The next time step, s: cfl / max_wave_rate(), or with rain falling
     * during it rain_limited_step() (FLUME_TimeStep.H), and no longer than
     * flume.max_dt. With discharge
     * segments the rate is at least that of their face states
     * (Inflow::max_rate), with the largest discharge from now to the end of
     * the longest step allowed: the target time or flume.max_dt, or sooner
     * the CFL step of the
     * cells; over a dry domain with neither, the rest of the hydrograph,
     * so that a later peak sets the step from the start (safe, and smaller
     * than needed). A domain with a discharge above 0 to come is not dry. A step that reaches
     * `target` (flume.stop_time in a stand-alone run, the end of the host's
     * step inside a host model), or falls short of it by round-off only (4
     * machine epsilons of the target), is set to end on it. With no wave anywhere
     * (a dry domain) `dry` is set: the step is the rain-limited one if rain
     * falls before the target (or ever, without one), and
     * otherwise goes straight to the target, or takes flume.max_dt; with
     * neither rain, a target nor flume.max_dt it aborts. Aborts as
     * check_state() does, and on a step that
     * is zero, not finite or too small to advance the time. `lands` is set
     * when the step ends at the target.
     */
    double compute_dt (bool& lands, bool& dry, const std::optional<double>& target) const;

    //! Abort if a depth is negative or not finite, or a wave speed is not
    //! finite, in any cell. Returns the largest wave rate (FLUME_TimeStep.H).
    amrex::Real check_state () const;

    //! Write <plot_file><step>: depth, discharge_x, discharge_y,
    //! bed_elevation, surface_elevation (bed + depth), then the fields of the
    //! sources (Source::plot_fields; infiltrated_depth with infiltration,
    //! drained_depth with drain inlets),
    //! then the flood maps (with flume.flood_maps = 1).
    void write_plotfile () const;

    //! write_plotfile() into the plotfile `name`.
    void write_plotfile (const std::string& name) const;

    [[nodiscard]] double time () const { return m_time; }
    [[nodiscard]] int step () const { return m_step; }
    [[nodiscard]] const amrex::MultiFab& state () const { return m_state; }
    [[nodiscard]] const amrex::MultiFab& bed () const { return m_bed; }
    [[nodiscard]] const amrex::Geometry& geom () const { return m_geom; }
    [[nodiscard]] const WaterBalance& balance () const { return *m_balance; }
    //! The rain a host model hands over (flume.rain.type = microphysics), or
    //! nullptr with any other rain.
    [[nodiscard]] HostRain* host_rain () { return dynamic_cast<HostRain*>(m_forcing.get()); }

private:
    //! One forward-Euler stage of `dt` on m_state, at the order of the run,
    //! in the step that ends at `time_end`, s (the inflow hydrographs give the
    //! mean discharge of the step). The water balance counts the outflow of
    //! the stage over `dt_weight`, s: the step, or half of it for each stage
    //! of the two-stage step.
    void advance_stage (amrex::Real dt, double dt_weight, double time_end);

    //! The water-balance line and the gauge lines of the current step.
    void write_logs ();

    //! Allocate the fields, fill the bed and the initial state, and set up
    //! the forcing, sources, inflow, balance, maps and gauges on m_ba, m_dm.
    void build ();

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
    //! The inflow hydrographs of the discharge segments (flume.bc.*)
    std::unique_ptr<Inflow> m_inflow;
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
    //! The step of the checkpoint the run continues from, or -1
    int m_restart_step = -1;
    //! The steps of the last plotfile, log line and checkpoint (-1: none)
    int m_last_plot = 0;
    int m_last_log = 0;
    int m_last_check = -1;
    //! Whether the last step found the domain dry (the message for a dry
    //! domain is printed once per dry spell)
    bool m_was_dry = false;
};

} // namespace flume

#endif
