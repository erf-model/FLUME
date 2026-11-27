/**
 * @file FLUME_ERFFloodEngine.cpp
 * @brief FLUME as a flood engine of ERF (erf.flood.engine = flume)
 *
 * ERF builds the flood grid (ERF_FloodGrid.H): the region its anchor level's
 * boxes cover (the whole domain on level 0, the refined region, such as a
 * catchment, on a finer level), refined in x and y by erf.flood.grid_ratio.
 * Its boxes are those of the anchor level and have the same owners, so a
 * flood box sits on the rank of the atmospheric columns above it. FLUME runs
 * on that grid with its own inputs (flume.*, as in a stand-alone run, less
 * the keys ERF sets: read_inputs, FLUME_Inputs.H).
 *
 * Each ERF step of the anchor level, advance() takes FLUME steps, each set by
 * FLUME's CFL number, rain and flume.max_dt, until FLUME lands on the end of
 * the ERF step. The rain is FLUME's own (flume.rain.*), or, with
 * flume.rain.type = microphysics, the surface precipitation of ERF's
 * microphysics over the ERF step: liquid and frozen (snow, graupel, hail) as
 * liquid water, each column's on the flood cells under it, at a rate held
 * constant over the ERF step (HostRain, FLUME_Forcing.H). With ERF's fixed step equal to
 * flume.max_dt, and FLUME's own step not shorter, every ERF step is one
 * FLUME step, and with prescribed rain the run is that of the stand-alone run
 * with the same flume.max_dt, to the last bit.
 *
 * ERF writes the flood plotfiles (erf.flood.plot_file, plot_int, plot_per):
 * FLUME's plotfile, named with ERF's step. ERF's checkpoints carry FLUME's
 * (under Level_<lev>/FLUME), and erf.restart continues from them. The ERF
 * plotfile can carry flood_depth, flood_discharge_x, flood_discharge_y,
 * flood_bed_elevation and flood_surface_elevation, averaged over each
 * atmospheric column. FLUME's water-balance and gauge files are written as in
 * a stand-alone run, every flume.log_int FLUME steps.
 */
#include <ERF_FloodEngineRegistry.H>
#include <ERF_FloodGrid.H>

#include "FLUME_Flood.H"
#include "FLUME_Inputs.H"
#include "FLUME_State.H"

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

//! The fields the ERF plotfile can carry, and where each comes from:
//! a component of the state, the bed, or the bed plus the depth.
enum class Field : int { depth, discharge_x, discharge_y, bed, surface };

const std::map<std::string, Field>& plot_fields ()
{
    static const std::map<std::string, Field> fields{
        {"flood_depth", Field::depth},
        {"flood_discharge_x", Field::discharge_x},
        {"flood_discharge_y", Field::discharge_y},
        {"flood_bed_elevation", Field::bed},
        {"flood_surface_elevation", Field::surface},
    };
    return fields;
}

class FLUMEFloodEngine : public FloodEngine
{
public:
    const char* name () const override { return "flume"; }

    // FLUME needs no more of ERF than the flood grid, which the coupling has
    // checked (verify_flood_prerequisites).
    void verify_prerequisites (const ERF& /*erf*/, int /*lev*/) override {}

    void initialize (const ERF& /*erf*/, int lev, const FloodGrid& fg,
                     const amrex::MultiFab* /*z_phys_nd_atm*/) override
    {
        m_lev = lev;
        m_grid = fg;
        m_in = flume::read_inputs(&fg.geom);
        m_flood = std::make_unique<flume::Flood>(m_in, fg.ba, fg.dm);
        // Before the plotfile ERF writes at step 0, as a stand-alone run
        // updates the maps before its first plotfile. A restart replaces them.
        m_flood->start_maps();
    }

    void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step) override
    {
        if (stop_time < std::numeric_limits<amrex::Real>::max()) {
            m_stop_time = static_cast<double>(stop_time);
        }
        if (max_step < std::numeric_limits<int>::max()) {
            m_max_steps = static_cast<long long>(max_step) * steps_per_coarse_step;
        }
        m_steps_per_coarse_step = steps_per_coarse_step;
    }

    // ERF's precipitation with flume.rain.type = microphysics (ERF then stops
    // at start-up unless its microphysics accumulates surface precipitation);
    // otherwise FLUME's own rain (flume.rain.*) drives it.
    bool wants_surface_precip () const override
    {
        return m_in.rain_type == flume::RainType::microphysics;
    }

    void advance (const FloodForcing& forcing) override
    {
        if (!m_started) {
            // A clean start, or a restart from a checkpoint without flood
            // state: FLUME starts from its initial state at ERF's time.
            if (!m_flood->restarted()) {
                m_flood->set_time(forcing.time);
                m_flood->start_logs();
            }
            m_started = true;
        }
        follow_clock(forcing);

        // ERF's time at the end of its step (in double, as ERF keeps it).
        const double t_end = forcing.time + static_cast<double>(forcing.dt);
        if (flume::HostRain* rain = m_flood->host_rain()) {
            set_rain(*rain, forcing, t_end);
        }
        while (m_flood->time() < t_end) { m_flood->step_to(t_end); }
        ++m_steps;

        const bool last_step = (m_max_steps >= 0 && m_steps >= m_max_steps)
            || (m_stop_time >= 0.0 && t_end >= m_stop_time - stop_tolerance(forcing));
        if (last_step && !m_finished) {
            m_flood->finish(false);
            m_finished = true;
        }
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        m_flood->write_checkpoint(checkpoint_dir(checkpointname));
    }

    void read_checkpoint (const std::string& restart_chkfile, int step) override
    {
        m_flood->restart(checkpoint_dir(restart_chkfile));
        m_started = true;
        m_steps = step;
    }

    void write_plotfile (const std::string& plotfile_prefix, amrex::Real /*time*/,
                         int step) const override
    {
        // FLUME's time and step go into the header; the name carries ERF's step.
        m_flood->check_state();
        m_flood->write_plotfile(amrex::Concatenate(plotfile_prefix, step, 5));
    }

    std::vector<std::string> plot_var_names () const override
    {
        std::vector<std::string> names;
        for (const auto& f : plot_fields()) { names.push_back(f.first); }
        return names;
    }

    const amrex::MultiFab* plot_var (const std::string& name) const override
    {
        const auto it = plot_fields().find(name);
        if (it == plot_fields().end() || !m_flood) { return nullptr; }
        const amrex::MultiFab& state = m_flood->state();
        const amrex::MultiFab& bed = m_flood->bed();
        auto& out = m_plot[name];
        if (!out) { out = std::make_unique<amrex::MultiFab>(state.boxArray(), state.DistributionMap(), 1, 0); }
        switch (it->second) {
        case Field::depth:
            amrex::MultiFab::Copy(*out, state, flume::State::h, 0, 1, 0);
            break;
        case Field::discharge_x:
            amrex::MultiFab::Copy(*out, state, flume::State::hu, 0, 1, 0);
            break;
        case Field::discharge_y:
            amrex::MultiFab::Copy(*out, state, flume::State::hv, 0, 1, 0);
            break;
        case Field::bed:
            amrex::MultiFab::Copy(*out, bed, 0, 0, 1, 0);
            break;
        case Field::surface:
            amrex::MultiFab::LinComb(*out, amrex::Real(1.0), bed, 0, amrex::Real(1.0), state,
                                     flume::State::h, 0, 1, 0);
            break;
        }
        return out.get();
    }

private:
    /**
     * The rain of this ERF step on the flood grid: the precipitation of each
     * column (kg/m2, the liquid and the frozen part) on every flood cell under
     * it, as a depth of liquid water, m (1 kg/m2 is 1 mm).
     */
    void set_rain (flume::HostRain& rain, const FloodForcing& forcing, double t_end)
    {
        if (forcing.precip_k0 == nullptr) {
            amrex::Abort("[FLOOD] flume: flume.rain.type = microphysics, but ERF handed over no"
                         " surface precipitation");
        }
        if (forcing.precip_k0->nComp() != 2) {
            amrex::Abort("[FLOOD] flume: ERF's surface precipitation has "
                         + std::to_string(forcing.precip_k0->nComp())
                         + " component(s); 2 (liquid, frozen) expected");
        }
        amrex::MultiFab& depth = rain.depth();
        if (!m_precip) {
            m_precip = std::make_unique<amrex::MultiFab>(depth.boxArray(), depth.DistributionMap(), 2, 0);
        }
        fill_flood_from_atm_k0(*m_precip, *forcing.precip_k0, m_grid, 0, 0);
        fill_flood_from_atm_k0(*m_precip, *forcing.precip_k0, m_grid, 1, 1);
        const amrex::Real kg_per_m2_to_m = amrex::Real(1.0e-3);
        const auto& p = m_precip->const_arrays();
        const auto& d = depth.arrays();
        amrex::ParallelFor(depth, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            d[b](i, j, k) = (p[b](i, j, k, 0) + p[b](i, j, k, 1)) * kg_per_m2_to_m;
        });
        amrex::Gpu::streamSynchronize();
        rain.set_step(forcing.time, t_end, "ERF's microphysics");
    }

    //! FLUME's checkpoint inside ERF's, on the flood grid's level.
    std::string checkpoint_dir (const std::string& checkpointname) const
    {
        return checkpointname + "/Level_" + std::to_string(m_lev) + "/FLUME";
    }

    /**
     * How far short of the stop time ERF's last step can end. ERF stops once
     * its coarse level is within 1e-6 of its step of the stop time, so a sum
     * of steps that do not add exactly in binary (0.1 s) can end that much
     * short. ERF hands over its stop time as Real: in single precision it is
     * rounded to float, by at most eps * stop, which is kept below half a
     * step so that the step before the last one does not end the run.
     */
    double stop_tolerance (const FloodForcing& forcing) const
    {
        const double dt = static_cast<double>(forcing.dt);
        const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
        return 1.0e-6 * dt * m_steps_per_coarse_step + std::min(eps * m_stop_time, 0.5 * dt);
    }

    /**
     * FLUME's time must be ERF's at the start of each ERF step: FLUME landed
     * on it at the end of the last one. In single precision ERF hands over
     * its step rounded to float, so the end FLUME landed on can differ from
     * ERF's by the rounding of the last step, which can be much longer than
     * this one (ERF shortens its last step onto the stop time); FLUME then
     * takes ERF's time. A larger difference means a step was missed or
     * repeated, and stops the run.
     */
    void follow_clock (const FloodForcing& forcing)
    {
        const double t = m_flood->time();
        const double dt_last = m_last_dt;
        m_last_dt = static_cast<double>(forcing.dt);
        if (t == forcing.time) { return; }
        const double tol = 1.0e-6 * std::max(dt_last, m_last_dt)
                           + 4.0 * std::numeric_limits<double>::epsilon() * std::abs(forcing.time);
        if (std::abs(t - forcing.time) > tol) {
            std::ostringstream msg;
            msg << std::setprecision(17) << "[FLOOD] flume: ERF's step starts at " << forcing.time
                << " s, but FLUME's time is " << t
                << " s: the flood engine missed or repeated a step of ERF";
            amrex::Abort(msg.str());
        }
        m_flood->set_time(forcing.time);
    }

    int m_lev = 0;
    flume::Inputs m_in;
    std::unique_ptr<flume::Flood> m_flood;
    bool m_started = false;   //!< whether the first advance (or a restart) has happened
    bool m_finished = false;  //!< whether the end of the run has been written
    long long m_steps = 0;    //!< steps of the flood grid's level so far
    long long m_max_steps = -1;  //!< the last of them, from ERF's max_step, or -1
    double m_stop_time = -1.0;   //!< ERF's stop time, elapsed, s, or -1
    int m_steps_per_coarse_step = 1;  //!< steps of the flood grid's level in one of level 0
    double m_last_dt = 0.0;      //!< ERF's step before this one, s
    //! The flood grid ERF built (flume.rain.type = microphysics reads through it)
    FloodGrid m_grid;
    //! This ERF step's precipitation on the flood grid, kg/m2: liquid, frozen
    std::unique_ptr<amrex::MultiFab> m_precip;
    //! The fields plot_var hands out, one per name, valid until the next advance
    mutable std::map<std::string, std::unique_ptr<amrex::MultiFab>> m_plot;
};

} // namespace

ERF_REGISTER_FLOOD_ENGINE(FLUMEFloodEngine, "flume")
