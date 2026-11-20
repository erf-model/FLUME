#include "FLUME_Flood.H"
#include "FLUME_Advance.H"
#include "FLUME_Boundary.H"
#include "FLUME_Init.H"
#include "FLUME_State.H"
#include "FLUME_TimeStep.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_CoordSys.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace flume {

Flood::Flood (const Inputs& in) : m_in(in)
{
    // One cell in z: the flood grid is two-dimensional.
    const amrex::Box domain(amrex::IntVect(0, 0, 0),
                            amrex::IntVect(in.n_cell[0] - 1, in.n_cell[1] - 1, 0));
    const amrex::RealBox real_box({in.prob_lo[0], in.prob_lo[1], amrex::Real(0.0)},
                                  {in.prob_hi[0], in.prob_hi[1], amrex::Real(1.0)});
    const amrex::Array<int, 3> is_periodic{{in.periodic[0], in.periodic[1], 0}};
    m_geom.define(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    m_ba.define(domain);
    m_ba.maxSize(amrex::IntVect(in.max_grid_size, in.max_grid_size, 1));
    m_dm.define(m_ba);

    // One layer of ghost cells in x and y, none in z.
    const amrex::IntVect n_ghost(1, 1, 0);
    m_bed.define(m_ba, m_dm, 1, n_ghost);
    m_state.define(m_ba, m_dm, State::ncomp, n_ghost);
    m_outflow_factor.define(m_ba, m_dm, 1, n_ghost);
    m_bed_source.define(m_ba, m_dm, 2, 0);
    for (int d = 0; d < 2; ++d) {
        m_flux[d].define(amrex::convert(m_ba, amrex::IntVect::TheDimensionVector(d)), m_dm,
                         State::ncomp, 0);
    }

    // The bed is checked before the state is built from it: a depth from a
    // bed that is not finite gives inf * 0 in the discharge.
    fill_plane_bed(m_bed, m_geom, in.bed_elevation, in.bed_slope);
    if (in.bed_bowl_depth != amrex::Real(0.0)) {
        add_bed_bowl(m_bed, m_geom, in.bed_bowl_depth, in.bed_bowl_radius);
    }
    if (in.bed_noise > amrex::Real(0.0)) { add_bed_noise(m_bed, m_geom, in.bed_noise); }
    const amrex::Long bad_bed = count_not_finite(m_bed);
    if (bad_bed > 0) {
        std::string msg = "The bed elevation is not finite in ";
        msg += std::to_string(bad_bed);
        msg += " cell(s): check flume.bed.elevation, flume.bed.slope, flume.bed.bowl_depth,"
               " flume.bed.bowl_radius and flume.bed.noise against the size of the domain";
        amrex::Abort(msg);
    }
    // The bed does not change: its ghost cells are filled once.
    fill_bed_ghost_cells(m_bed, m_geom);

    if (in.init_dam_position.has_value()) {
        fill_state_dam(m_state, m_bed, m_geom, in.init_dam_direction, *in.init_dam_position,
                       *in.init_surface_elevation, in.init_dam_surface_elevation,
                       in.init_velocity, in.init_dam_velocity);
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
               ", flume.init.surface_curvature, flume.init.dam_surface_elevation"
               ", flume.init.velocity and flume.init.dam_velocity against the bed"
               " (flume.bed.elevation, flume.bed.slope, flume.bed.bowl_depth,"
               " flume.bed.bowl_radius, flume.bed.noise)";
        amrex::Abort(msg);
    }

    amrex::Print() << "FLUME: " << in.n_cell[0] << " x " << in.n_cell[1] << " cells in "
                   << m_ba.size() << " box(es)\n";
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
               " flume.init.surface_slope, flume.init.surface_curvature); later, the outflow"
               " limit should have kept the depth at zero or above: lower flume.cfl and"
               " report the deck";
        amrex::Abort(msg);
    }
    return w.rate;
}

double
Flood::compute_dt (bool& lands, bool& dry) const
{
    lands = false;
    dry = false;
    const amrex::Real rate = check_state();
    if (!(rate > amrex::Real(0.0))) {
        if (!m_in.stop_time.has_value()) {
            amrex::Abort("No cell has a wave speed above zero (no depth is above flume.dry_depth),"
                         " so nothing sets the time step: set flume.stop_time");
        }
        lands = true;
        dry = true;
        return *m_in.stop_time - m_time;
    }
    double dt = static_cast<double>(m_in.cfl / rate);
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
Flood::evolve ()
{
    write_plotfile();
    int last_plot = m_step;

    const bool limited = m_in.max_step.has_value() || m_in.stop_time.has_value();
    if (!limited) {
        amrex::Print() << "Neither flume.stop_time nor flume.max_step is set: no step is taken\n";
    }
    while (limited && (!m_in.max_step.has_value() || m_step < *m_in.max_step)
           && (!m_in.stop_time.has_value() || m_time < *m_in.stop_time)) {
        bool lands = false;
        bool dry = false;
        const double dt = compute_dt(lands, dry);

        // A dry domain stays as it is, and its step (the whole time left) is
        // bounded by nothing: only the time advances. Otherwise the state
        // advances with the step in its own precision; in single precision
        // that differs from dt by rounding.
        if (!dry) {
            // The update leaves the ghost cells of the state out of date.
            fill_state_ghost_cells(m_state, m_geom);
            const auto dt_real = static_cast<amrex::Real>(dt);
            compute_fluxes(m_flux, m_state, m_bed, m_in.gravity, m_in.dry_depth);
            limit_outflow(m_flux, m_outflow_factor, m_state, m_geom, dt_real);
            compute_bed_source(m_bed_source, m_state, m_bed, m_in.gravity, m_in.dry_depth);
            update_state(m_state, m_flux, m_bed_source, m_geom, m_in.dry_depth, dt_real);
        } else {
            amrex::Print() << "No depth is above flume.dry_depth (" << m_in.dry_depth
                           << " m): the state is not advanced\n";
        }

        m_time = lands ? *m_in.stop_time : m_time + dt;
        ++m_step;
        amrex::Print() << "Step " << m_step << ": time " << m_time << " s, dt " << dt << " s\n";

        if (m_in.plot_int > 0 && m_step % m_in.plot_int == 0) {
            // A state that the next step would stop on is not written.
            check_state();
            write_plotfile();
            last_plot = m_step;
        }
    }

    if (last_plot != m_step) {
        check_state();
        write_plotfile();
    }
    amrex::Print() << "FLUME: " << m_step << " step(s), time " << m_time << " s\n";
}

} // namespace flume
