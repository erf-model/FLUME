#include "FLUME_Flood.H"
#include "FLUME_Checkpoint.H"
#include "FLUME_Init.H"
#include "FLUME_State.H"

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace flume {

void
Flood::write_plotfile () const
{
    // depth, discharge_x, discharge_y, bed_elevation, surface_elevation
    amrex::Vector<std::string> names{"depth", "discharge_x", "discharge_y", "bed_elevation",
                                   "surface_elevation"};
    std::vector<PlotField> extra;
    for (const auto& s : m_sources) {
        for (const auto& f : s->plot_fields()) { extra.push_back(f); }
    }
    if (m_maps) {
        for (const auto& f : m_maps->plot_fields()) { extra.push_back(f); }
    }
    const int i_bed = State::ncomp;
    const int i_surface = State::ncomp + 1;
    const int nvar = State::ncomp + 2 + static_cast<int>(extra.size());
    amrex::MultiFab out(m_ba, m_dm, nvar, 0);
    amrex::MultiFab::Copy(out, m_state, 0, 0, State::ncomp, 0);
    amrex::MultiFab::Copy(out, m_bed, 0, i_bed, 1, 0);
    // The plotfile holds the time in the precision of the build.
    // Water surface on the absolute datum: bed + depth. In a dry cell it is
    // the bed elevation.
    amrex::MultiFab::LinComb(out, amrex::Real(1.0), m_bed, 0, amrex::Real(1.0), m_state, State::h,
                             i_surface, 1, 0);
    for (std::size_t n = 0; n < extra.size(); ++n) {
        amrex::MultiFab::Copy(out, *extra[n].second, 0, i_surface + 1 + static_cast<int>(n), 1, 0);
        names.push_back(extra[n].first);
    }

    const std::string name = amrex::Concatenate(m_in.plot_file, m_step, 5);
    amrex::WriteSingleLevelPlotfile(name, out, names, m_geom, static_cast<amrex::Real>(m_time), m_step);
    amrex::Print() << "Wrote " << name << " at time " << m_time << " s\n";
}

namespace {

std::vector<std::string> source_names (const std::vector<std::unique_ptr<Source>>& sources)
{
    std::vector<std::string> v;
    for (const auto& s : sources) { v.push_back(s->name()); }
    return v;
}

std::string list (const std::vector<std::string>& v)
{
    std::string s;
    for (const auto& w : v) { s += (s.empty() ? "" : " ") + w; }
    return s.empty() ? std::string("none") : s;
}

} // namespace

void
Flood::write_checkpoint () const
{
    const std::string name = amrex::Concatenate(m_in.check_file, m_step, 5);
    CheckpointWriter chk(name);
    chk.put_ints("step", {m_step});
    chk.put("time", {m_time});
    const amrex::Box& domain = m_geom.Domain();
    chk.put_ints("n_cell", {domain.length(0), domain.length(1)});
    chk.put("prob_lo", {static_cast<double>(m_geom.ProbLo(0)),
                        static_cast<double>(m_geom.ProbLo(1))});
    chk.put("prob_hi", {static_cast<double>(m_geom.ProbHi(0)),
                        static_cast<double>(m_geom.ProbHi(1))});
    chk.put_ints("periodic", {m_geom.isPeriodic(0) ? std::int64_t{1} : std::int64_t{0}, m_geom.isPeriodic(1) ? std::int64_t{1} : std::int64_t{0}});
    chk.put_words("sources", source_names(m_sources));
    chk.put_multifab("state", m_state);
    chk.put_multifab("bed", m_bed);
    for (const auto& s : m_sources) { s->write_checkpoint(chk); }
    m_balance->write_checkpoint(chk, m_sources);
    chk.put_words("gauges", m_gauges->words());
    if (m_maps) { m_maps->write_checkpoint(chk); }
    chk.finish();
    amrex::Print() << "Wrote " << name << " at time " << m_time << " s\n";
}

void
Flood::restart (const std::string& dir)
{
    const CheckpointReader chk(dir);
    const amrex::Box& domain = m_geom.Domain();
    const std::vector<std::int64_t> n_cell = chk.get_ints("n_cell", 2);
    if (n_cell[0] != domain.length(0) || n_cell[1] != domain.length(1)) {
        chk.fail("its grid is " + std::to_string(n_cell[0]) + " x " + std::to_string(n_cell[1])
                 + " cells, not the " + std::to_string(domain.length(0)) + " x "
                 + std::to_string(domain.length(1)) + " of flume.n_cell");
    }
    const std::vector<double> lo = chk.get("prob_lo", 2);
    const std::vector<double> hi = chk.get("prob_hi", 2);
    for (int d = 0; d < 2; ++d) {
        if (lo[d] != static_cast<double>(m_geom.ProbLo(d))
            || hi[d] != static_cast<double>(m_geom.ProbHi(d))) {
            chk.fail("its domain differs from that of flume.prob_lo and flume.prob_hi");
        }
    }
    const std::vector<std::int64_t> periodic = chk.get_ints("periodic", 2);
    for (int d = 0; d < 2; ++d) {
        if (periodic[d] != (m_geom.isPeriodic(d) ? std::int64_t{1} : std::int64_t{0})) {
            chk.fail("its periodic sides differ from those of flume.periodic");
        }
    }
    const std::vector<std::string> names = source_names(m_sources);
    if (chk.get_words("sources") != names) {
        chk.fail("its sources (" + list(chk.get_words("sources")) + ") differ from those of the"
                 " inputs (" + list(names) + ")");
    }

    const std::vector<std::int64_t> step = chk.get_ints("step", 1);
    if (step[0] < 0 || step[0] > std::numeric_limits<int>::max()) { chk.fail("step is not valid"); }
    m_step = static_cast<int>(step[0]);
    m_time = chk.get("time", 1)[0];
    if (!std::isfinite(m_time) || m_time < 0.0) { chk.fail("time is not valid"); }

    // The bed is built from the inputs again; it must be the one the
    // checkpoint was run on, to the bit.
    {
        amrex::MultiFab bed(m_ba, m_dm, 1, 0);
        chk.get_multifab("bed", bed);
        amrex::MultiFab::Subtract(bed, m_bed, 0, 0, 1, 0);
        if (bed.norm0(0, 0) != amrex::Real(0.0)) {
            chk.fail("its bed differs from the bed of the inputs (flume.bed.*,"
                     " flume.buildings.*)");
        }
    }
    // A gauge file continues only at the same cell (a restart without the
    // log writes no gauge file).
    if (!m_gauges->empty() && chk.get_words("gauges") != m_gauges->words()) {
        chk.fail("its gauges (" + list(chk.get_words("gauges")) + ") differ from those of the"
                 " inputs (" + list(m_gauges->words()) + "): name:i:j of the cell of each,"
                 " from flume.gauges.*");
    }
    chk.get_multifab("state", m_state);
    {
        const amrex::Long bad = count_not_finite(m_state);
        if (bad > 0) {
            chk.fail("its depth or discharge is not finite in " + std::to_string(bad)
                     + " cell(s)");
        }
    }
    // The first stage fills the ghost cells again; they are filled here so
    // that no ghost cell holds what the state had before the restart.
    m_boundary.fill_state(m_state);
    for (auto& s : m_sources) { s->read_checkpoint(chk); }
    // The run without a break has a log line at this step only when it is a
    // multiple of flume.log_int; the first part also wrote one at its last
    // step, which may not be.
    const bool keep_step = m_in.log_int > 0 && m_step % m_in.log_int == 0;
    m_balance->resume(chk, m_state, m_sources, m_step, keep_step);
    if (m_maps) { m_maps->read_checkpoint(chk); }
    if (m_in.log_int > 0) { m_gauges->resume(m_step, keep_step); }
    m_restart_step = m_step;
    std::ostringstream s;
    s << std::setprecision(17) << "FLUME: restarted from " << dir << " at step " << m_step
      << ", time " << m_time << " s\n";
    amrex::Print() << s.str();
}

} // namespace flume
