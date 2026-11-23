#ifndef FLUME_WATERBALANCE_H_
#define FLUME_WATERBALANCE_H_

#include "FLUME_Inputs.H"
#include "FLUME_Source.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace flume {

//! The area of a cell, m2, in double precision.
double cell_area (const amrex::Geometry& geom);

//! The sum of component 0 of the valid cells of `a` over this rank's boxes,
//! in double precision. No MPI reduction.
double local_sum (const amrex::MultiFab& a);

//! The volume, m3, of a sum of depths `depth_sum`, m, over cells of area
//! `area`, m2: 0 for a sum of 0, also when the area overflows (inf * 0 would
//! trap with amrex.fpe_trap_invalid = 1).
inline double volume (double area, double depth_sum)
{
    return (depth_sum == 0.0) ? 0.0 : area * depth_sum;
}

//! Add `x` to the compensated (Kahan) sum `sum`, whose rounding is carried
//! in `carry`: the rounding of a total built from many steps then does not
//! grow with the number of steps.
inline void compensated_add (double& sum, double& carry, double x)
{
    const double y = x - carry;
    const double t = sum + y;
    carry = (t - sum) - y;
    sum = t;
}

/**
 * The water balance of the run, written to water_balance.txt
 * (Docs/sphinx_doc/Running.rst).
 *
 * The volume of water on the surface, V, changes only through the sources
 * (rain in, infiltration out) and through the faces of the outflow segments.
 * With V_0 the volume at the start, S_n the volume source n has added since
 * the start (below 0 for a sink) and Q_m the volume that has left through
 * segment m (below 0 for water that entered), the residual
 *
 *   R = V - V_0 - sum_n S_n + sum_m Q_m
 *
 * is zero in exact arithmetic. Each term is counted where it happens, not
 * from V: S_n by the source itself (Source::balance_volume), Q_m from the
 * mass flux the update uses on each face of the segment. So R measures how
 * well the update, the outflow limit and the sources keep the water they
 * account for: a residual above round-off is a defect. Walls are not
 * counted: the mass flux through a wall face is zero, up to the rounding of
 * one product (limit_outflow), and what it carries shows in R.
 *
 * Each line of the log holds the step, the time, V, each source's column
 * (water added, or water removed for a sink), each segment's column (water
 * out) and R, all in m3 and cumulative from the start.
 */
class WaterBalance
{
public:
    //! The outflow segments of `in` on the grid of `geom`. Opens the log on
    //! the I/O rank when `in.log_int` is above 0.
    WaterBalance (const Inputs& in, const amrex::Geometry& geom);

    //! Add the water that the fluxes `flux` (after limit_outflow) carry out
    //! through each outflow segment in `dt`, s, to this rank's totals: the
    //! mass flux on each boundary face of the segment, times the face area
    //! and `dt`. A stage of the two-stage step passes half the step. No MPI.
    void add_outflow (const amrex::Array<amrex::MultiFab, 2>& flux, double dt);

    //! The volume on the surface at the start (V_0), the header and the line
    //! of step 0.
    void start (const amrex::MultiFab& state,
                const std::vector<std::unique_ptr<Source>>& sources);

    //! One line of the log at the end of step `step`, time `time`, s. Sums
    //! over the ranks (one reduction).
    void log (const amrex::MultiFab& state, const std::vector<std::unique_ptr<Source>>& sources,
              int step, double time);

    //! Print the totals of the last line logged.
    void print_summary () const;

    //! The residual of the last line logged, m3.
    [[nodiscard]] double residual () const { return m_residual; }

private:
    //! One outflow segment: the boundary faces of cells lo to hi along the
    //! side, on the low (high = 0) or high side across direction d.
    struct Segment
    {
        std::string name;
        int d = 0;
        int high = 0;
        int lo = 0;
        int hi = -1;
    };

    //! The totals over all ranks: V, the sources, the segments.
    std::vector<double> totals (const amrex::MultiFab& state,
                                const std::vector<std::unique_ptr<Source>>& sources) const;
    void write_line (int step, double time, const std::vector<double>& t);

    amrex::Geometry m_geom;
    double m_cell_area = 0.0;
    std::vector<Segment> m_segments;
    std::vector<double> m_out;  //!< volume out through each segment, this rank, m3
    std::vector<double> m_out_carry;  //!< rounding carried in the sums m_out
    std::vector<std::string> m_source_names;
    std::vector<bool> m_source_removes;
    double m_v0 = 0.0;
    double m_residual = 0.0;
    std::vector<double> m_last;  //!< totals of the last line logged
    bool m_on = false;
    std::ofstream m_file;        //!< open on the I/O rank only
};

} // namespace flume

#endif
