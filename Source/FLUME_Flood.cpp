#include "FLUME_Flood.H"
#include "FLUME_Advance.H"
#include "FLUME_Boundary.H"
#include "FLUME_Buildings.H"
#include "FLUME_Init.H"
#include "FLUME_Raster.H"
#include "FLUME_State.H"
#include "FLUME_TimeStep.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_CoordSys.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace flume {

namespace {

//! The grid of the run. One cell in z: the flood grid is two-dimensional.
amrex::Geometry
make_geometry (const Inputs& in)
{
    const amrex::Box domain(amrex::IntVect(0, 0, 0),
                            amrex::IntVect(in.n_cell[0] - 1, in.n_cell[1] - 1, 0));
    const amrex::RealBox real_box({in.prob_lo[0], in.prob_lo[1], amrex::Real(0.0)},
                                  {in.prob_hi[0], in.prob_hi[1], amrex::Real(1.0)});
    const amrex::Array<int, 3> is_periodic{{in.periodic[0], in.periodic[1], 0}};
    return amrex::Geometry(domain, real_box, amrex::CoordSys::cartesian, is_periodic);
}

} // namespace

Flood::Flood (const Inputs& in) : m_in(in), m_boundary(in, make_geometry(in))
{
    m_geom = make_geometry(in);
    const amrex::Box& domain = m_geom.Domain();

    m_ba.define(domain);
    m_ba.maxSize(amrex::IntVect(in.max_grid_size, in.max_grid_size, 1));
    m_dm.define(m_ba);

    // Ghost cells in x and y, none in z: a face reads the cell on each side
    // of it, and at second order also the neighbours of those two cells.
    const amrex::IntVect n_ghost(in.order, in.order, 0);
    m_bed.define(m_ba, m_dm, 1, n_ghost);
    m_state.define(m_ba, m_dm, State::ncomp, n_ghost);
    if (in.order == 2) { m_state_old.define(m_ba, m_dm, State::ncomp, 0); }
    m_outflow_factor.define(m_ba, m_dm, 1, amrex::IntVect(1, 1, 0));
    m_bed_source.define(m_ba, m_dm, 2, 0);
    for (int d = 0; d < 2; ++d) {
        m_flux[d].define(amrex::convert(m_ba, amrex::IntVect::TheDimensionVector(d)), m_dm,
                         State::ncomp, 0);
    }

    // The bed is checked before the state is built from it: a depth from a
    // bed that is not finite gives inf * 0 in the discharge.
    if (!in.bed_file.empty()) {
        const Raster dem = read_raster("bed.file", in.bed_file, m_geom);
        sample_bilinear(m_bed, 0, dem, m_geom);
        amrex::Print() << "FLUME: bed from " << in.bed_file << " (" << dem.nx << " x " << dem.ny
                       << " raster cells read)\n";
    } else {
        fill_plane_bed(m_bed, m_geom, in.bed_elevation, in.bed_slope);
        if (in.bed_bowl_depth != amrex::Real(0.0)) {
            add_bed_bowl(m_bed, m_geom, in.bed_bowl_depth, in.bed_bowl_radius);
        }
        if (in.bed_noise > amrex::Real(0.0)) { add_bed_noise(m_bed, m_geom, in.bed_noise); }
    }
    if (!in.buildings_file.empty()) {
        const Raster b = read_raster("buildings.file", in.buildings_file, m_geom);
        m_buildings.define(m_ba, m_dm, 1, 0);
        amrex::MultiFab& height = m_buildings;
        sample_nearest(height, 0, b, m_geom);
        const amrex::Long below = count_below_zero(height);
        if (below > 0) {
            amrex::Abort("flume.buildings.file (" + in.buildings_file + "): "
                         + std::to_string(below) + " cell(s) of the domain take a building height"
                         " below 0");
        }
        const amrex::Long n_building = in.buildings_height.has_value()
            ? add_building_footprints(m_bed, height, *in.buildings_height)
            : add_building_heights(m_bed, height);
        amrex::Print() << "FLUME: buildings from " << in.buildings_file << ": " << n_building
                       << " cell(s) raised\n";
    }
    const amrex::Long bad_bed = count_not_finite(m_bed);
    if (bad_bed > 0) {
        std::string msg = "The bed elevation is not finite in ";
        msg += std::to_string(bad_bed);
        msg += " cell(s): check flume.bed.elevation, flume.bed.slope, flume.bed.bowl_depth,"
               " flume.bed.bowl_radius and flume.bed.noise against the size of the domain"
               " (or the values of flume.bed.file and flume.buildings.file)";
        amrex::Abort(msg);
    }
    // The bed does not change: its ghost cells are filled once.
    const amrex::Real zero = amrex::Real(0.0);
    const amrex::Array<amrex::Real, 2> continued_slope =
        (in.bed_continue_slope == 1) ? in.bed_slope : amrex::Array<amrex::Real, 2>{{zero, zero}};
    m_boundary.fill_bed(m_bed, continued_slope);

    if (in.init_dam_position.has_value()) {
        fill_state_dam(m_state, m_bed, m_geom, in.init_dam_direction, *in.init_dam_position,
                       *in.init_surface_elevation, in.init_dam_surface_elevation,
                       in.init_velocity, in.init_dam_velocity);
    } else if (in.init_wave_amplitude != amrex::Real(0.0)) {
        fill_state_wave(m_state, m_bed, m_geom, in.init_wave_direction,
                        *in.init_surface_elevation, in.init_wave_amplitude, in.init_velocity,
                        in.gravity);
    } else if (in.init_surface_elevation.has_value()) {
        const bool level = in.init_surface_slope[0] == amrex::Real(0.0)
                           && in.init_surface_slope[1] == amrex::Real(0.0)
                           && in.init_surface_curvature == amrex::Real(0.0);
        if (level) {
            fill_state_from_surface(m_state, m_bed, *in.init_surface_elevation, in.init_velocity);
        } else {
            fill_state_from_curved_surface(m_state, m_bed, m_geom, *in.init_surface_elevation,
                                           in.init_surface_slope, in.init_surface_curvature,
                                           in.init_velocity);
        }
    } else {
        m_state.setVal(amrex::Real(0.0)); // dry
    }
    const amrex::Long bad_state = count_not_finite(m_state);
    if (bad_state > 0) {
        std::string msg = "The initial depth or discharge is not finite in ";
        msg += std::to_string(bad_state);
        msg += " cell(s): check flume.init.surface_elevation, flume.init.surface_slope"
               ", flume.init.surface_curvature, flume.init.wave_amplitude"
               ", flume.init.dam_surface_elevation"
               ", flume.init.velocity and flume.init.dam_velocity against the bed"
               " (flume.bed.elevation, flume.bed.slope, flume.bed.bowl_depth,"
               " flume.bed.bowl_radius, flume.bed.noise)";
        amrex::Abort(msg);
    }

    amrex::Print() << "FLUME: " << in.n_cell[0] << " x " << in.n_cell[1] << " cells in "
                   << m_ba.size() << " box(es)\n";

    m_forcing = make_forcing(in, m_geom, m_ba, m_dm);
    m_sources = make_sources(in, m_geom, m_ba, m_dm, m_forcing.get(),
                             m_buildings.ok() ? &m_buildings : nullptr);
    m_inflow = std::make_unique<Inflow>(in, m_geom);
    m_balance = std::make_unique<WaterBalance>(in, m_geom);
    if (in.flood_maps == 1) { m_maps = std::make_unique<FloodMaps>(in, m_ba, m_dm); }
    m_gauges = std::make_unique<Gauges>(in, m_geom);
    if (!in.restart.empty()) { restart(in.restart); }
}

amrex::Real
Flood::check_state () const
{
    const WaveRate w = max_wave_rate(m_state, m_geom, m_in.gravity, m_in.dry_depth);
    if (w.bad_cells > 0) {
        std::string msg = "The depth is negative or not finite, or the wave speed is not finite, in ";
        msg += std::to_string(w.bad_cells);
        msg += " cell(s) at step ";
        msg += std::to_string(m_step);
        msg += ": at step 0 the initial depth or velocity (flume.init.surface_elevation,"
               " flume.init.dam_surface_elevation, flume.init.velocity, flume.init.dam_velocity)"
               " or flume.gravity is too large (with a tilted or curved surface also"
               " flume.init.surface_slope, flume.init.surface_curvature, and with a wave"
               " flume.init.wave_amplitude); later, the outflow"
               " limit should have kept the depth at zero or above: lower flume.cfl and"
               " report the deck";
        amrex::Abort(msg);
    }
    return w.rate;
}

void
Flood::advance_stage (amrex::Real dt, double dt_weight, double time_end)
{
    // The update leaves the ghost cells of the state out of date.
    m_boundary.fill_state(m_state);
    if (m_in.order == 2) {
        compute_fluxes_second_order(m_flux, m_state, m_bed, m_in.gravity, m_in.dry_depth);
    } else {
        compute_fluxes(m_flux, m_state, m_bed, m_in.gravity, m_in.dry_depth);
    }
    // Before the outflow limit, which leaves water that enters as it is.
    if (!m_inflow->empty()) {
        m_inflow->set_fluxes(m_flux, m_state, m_time, time_end, static_cast<double>(dt));
    }
    limit_outflow(m_flux, m_outflow_factor, m_state, m_geom, dt);
    m_balance->add_outflow(m_flux, dt_weight);
    if (m_in.order == 2) {
        compute_bed_source_second_order(m_bed_source, m_state, m_bed, m_in.gravity,
                                        m_in.dry_depth);
    } else {
        compute_bed_source(m_bed_source, m_state, m_bed, m_in.gravity, m_in.dry_depth);
    }
    update_state(m_state, m_flux, m_bed_source, m_geom, m_in.dry_depth, dt);
}

double
Flood::compute_dt (bool& lands, bool& dry) const
{
    lands = false;
    dry = false;
    amrex::Real rate = check_state();
    // The longest the step can be without rain: to the stop time, and no
    // longer than the CFL step.
    double dt_max = std::numeric_limits<double>::infinity();
    if (m_in.stop_time.has_value()) { dt_max = *m_in.stop_time - m_time; }
    // The face states of the discharge segments, with the largest discharge
    // up to then, take part as cells do: a face state is deeper and faster
    // with more discharge, so the bound holds over any shorter step. A
    // discharge to come makes the domain wet.
    if (!m_inflow->empty()) {
        if (rate > amrex::Real(0.0)) {
            dt_max = std::min(dt_max, static_cast<double>(m_in.cfl / rate));
        }
        rate = std::max(rate, m_inflow->max_rate(m_state, m_time, m_time + dt_max));
    }
    dry = !(rate > amrex::Real(0.0));
    if (!dry) { dt_max = std::min(dt_max, static_cast<double>(m_in.cfl / rate)); }
    // The largest rain rate over that span bounds the rain of the step, which
    // is no longer.
    const double r_max = (m_forcing != nullptr)
        ? m_forcing->max_rain_rate(m_time, m_time + dt_max) : 0.0;
    double dt = dt_max;
    if (r_max > 0.0) {
        dt = std::min(dt_max, rain_limited_step(static_cast<double>(rate), r_max,
                                                static_cast<double>(m_in.gravity),
                                                static_cast<double>(m_in.dry_depth),
                                                static_cast<double>(m_geom.InvCellSize(0)),
                                                static_cast<double>(m_geom.InvCellSize(1)),
                                                static_cast<double>(m_in.cfl)));
    } else if (dry) {
        if (!m_in.stop_time.has_value()) {
            amrex::Abort("No cell has a wave speed above zero (no depth is above flume.dry_depth),"
                         " and no rain falls and no water enters through a discharge segment"
                         " from this time on, so nothing sets the time step: set flume.stop_time");
        }
        lands = true;
        return *m_in.stop_time - m_time;
    }
    if (!(std::isfinite(dt) && m_time + dt > m_time)) {
        std::ostringstream msg;
        msg << std::setprecision(17) << "The time step, " << dt
            << " s, is zero, not finite or too small to advance the time, " << m_time
            << " s, at step " << m_step
            << ": check flume.cfl, the cell size, flume.gravity and the initial state";
        amrex::Abort(msg.str());
    }
    if (m_in.stop_time.has_value()) {
        // Land on the stop time when the step reaches it, or falls short by
        // round-off of the summed steps only: otherwise a last step of the
        // size of that round-off would follow.
        const double stop = *m_in.stop_time;
        const double tol = 4.0 * std::numeric_limits<double>::epsilon() * stop;
        if (m_time + dt >= stop - tol) {
            dt = stop - m_time;
            lands = true;
        }
    }
    return dt;
}

void
Flood::write_logs ()
{
    m_balance->log(m_state, m_sources, m_step, m_time);
    m_gauges->log(m_state, m_bed, m_step, m_time);
}

void
Flood::evolve ()
{
    // A restart has written its logs, plotfile and checkpoint at this step.
    const bool restarted = m_restart_step >= 0;
    if (!restarted) {
        if (m_maps) { m_maps->update(m_state, m_time); }
        write_plotfile();
        m_balance->start(m_state, m_sources);
        if (m_in.log_int > 0) {
            m_gauges->start();
            m_gauges->log(m_state, m_bed, m_step, m_time);
        }
    }
    int last_plot = m_step;
    int last_log = m_step;
    int last_check = restarted ? m_step : -1;

    const bool limited = m_in.max_step.has_value() || m_in.stop_time.has_value();
    if (!limited) {
        amrex::Print() << "Neither flume.stop_time nor flume.max_step is set: no step is taken\n";
    }
    bool was_dry = false;  // the message for a dry domain is printed once per dry spell
    while (limited && (!m_in.max_step.has_value() || m_step < *m_in.max_step)
           && (!m_in.stop_time.has_value() || m_time < *m_in.stop_time)) {
        bool lands = false;
        bool dry = false;
        const double dt = compute_dt(lands, dry);

        // The fluxes leave a dry domain as it is: only the sources act on it.
        // Otherwise the state advances with the step in its own precision;
        // in single precision that differs from dt by rounding.
        const auto dt_real = static_cast<amrex::Real>(dt);
        const double time_end = lands ? *m_in.stop_time : m_time + dt;
        if (!dry) {
            if (m_in.order == 2) {
                // Two stages from the state, then the mean with it.
                // The step in the precision of the state, as the update
                // uses it, for the outflow of the balance.
                const double half = 0.5 * static_cast<double>(dt_real);
                amrex::MultiFab::Copy(m_state_old, m_state, 0, 0, State::ncomp, 0);
                advance_stage(dt_real, half, time_end);
                // The second stage does not start from a state that a step
                // would stop on.
                check_state();
                advance_stage(dt_real, half, time_end);
                average_with(m_state, m_state_old, m_in.dry_depth);
            } else {
                advance_stage(dt_real, static_cast<double>(dt_real), time_end);
            }
        } else {
            // Every cell is dry: zero its discharges, as a step with fluxes
            // does in dry cells, so that a cell the rain makes wet starts
            // still (the initial state may hold a discharge in a dry cell).
            m_state.setVal(amrex::Real(0.0), State::hu, 2, 0);
            if (!was_dry) {
                amrex::Print() << "No depth is above flume.dry_depth (" << m_in.dry_depth
                               << " m): the fluxes are not computed until one is\n";
            }
        }
        was_dry = dry;
        // The sources act once, after the explicit terms of the whole step,
        // in order, from the start to the end of the step. Friction is
        // implicit in the discharge, with the depth at the end of the step
        // (FLUME_Friction.H): acting on each stage instead, and then taking
        // the mean, would remove only half of a stiff film's excess
        // discharge per step. They act on a dry domain too: rain falls on
        // it, and friction leaves a dry cell as it is.
        const SourceContext ctx{m_bed, m_geom, m_time, time_end};
        for (auto& s : m_sources) { s->apply(m_state, ctx, dt_real); }

        m_time = time_end;
        ++m_step;
        amrex::Print() << "Step " << m_step << ": time " << m_time << " s, dt " << dt << " s\n";

        // A state that the next step would stop on is not read: the maps and
        // logs would take the square root of, or divide, an overflowed depth.
        const bool log_now = m_in.log_int > 0 && m_step % m_in.log_int == 0;
        const bool plot_now = m_in.plot_int > 0 && m_step % m_in.plot_int == 0;
        const bool check_now = m_in.check_int > 0 && m_step % m_in.check_int == 0;
        if (m_maps || log_now || plot_now || check_now) { check_state(); }
        if (m_maps) { m_maps->update(m_state, m_time); }
        if (log_now) {
            write_logs();
            last_log = m_step;
        }

        if (plot_now) {
            write_plotfile();
            last_plot = m_step;
        }
        if (check_now) {
            write_checkpoint();
            last_check = m_step;
        }
    }

    if (last_plot != m_step) {
        check_state();
        write_plotfile();
    }
    if (m_in.log_int > 0 && last_log != m_step) { write_logs(); }
    if (m_in.check_int > 0 && last_check != m_step) { write_checkpoint(); }
    m_balance->print_summary();
    amrex::Print() << "FLUME: " << m_step << " step(s), time " << m_time << " s\n";
}

} // namespace flume
