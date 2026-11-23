#include "FLUME_Flood.H"
#include "FLUME_State.H"

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <cstddef>
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

} // namespace flume
