#ifndef FLUME_GAUGES_H_
#define FLUME_GAUGES_H_

#include "FLUME_Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace flume {

/**
 * Time series at gauge points (flume.gauges.*), one file per gauge,
 * gauge_<name>.txt (Docs/sphinx_doc/Running.rst).
 *
 * A gauge reads the cell that holds its point; there is no interpolation.
 * A point on a face between two cells (within 4 machine epsilons of the
 * largest coordinate, and at most 1e-3 of a cell) takes the cell above it in
 * that direction, and a point on
 * the high edge of the domain the last cell. Each
 * line holds the step, the time (s), the depth (m), the surface elevation
 * (m), and the velocity in x and y (m/s; 0 in a dry cell, depth at or below
 * flume.dry_depth). The lines are written with those of the water-balance
 * log (flume.log_int).
 */
class Gauges
{
public:
    //! Find the cell of each gauge of `in`; open the files on the I/O rank.
    Gauges (const Inputs& in, const amrex::Geometry& geom);

    //! One line per gauge at step `step`, time `time`, s. Every rank calls
    //! it (one reduction).
    void log (const amrex::MultiFab& state, const amrex::MultiFab& bed, int step, double time);

    [[nodiscard]] bool empty () const { return m_names.empty(); }

    //! The cell of gauge n.
    [[nodiscard]] const amrex::IntVect& cell (std::size_t n) const { return m_cells[n]; }

private:
    std::vector<std::string> m_names;
    std::vector<amrex::IntVect> m_cells;
    amrex::Real m_dry_depth = amrex::Real(0.0);
    std::vector<std::unique_ptr<std::ofstream>> m_files;  //!< open on the I/O rank only
};

//! The cell that holds the point `x`, m, on the grid of `geom` (as Gauges).
amrex::IntVect gauge_cell (const amrex::Array<amrex::Real, 2>& x, const amrex::Geometry& geom);

} // namespace flume

#endif
