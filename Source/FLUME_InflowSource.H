#ifndef FLUME_INFLOWSOURCE_H_
#define FLUME_INFLOWSOURCE_H_

#include "FLUME_Inputs.H"
#include "FLUME_Source.H"
#include "FLUME_TimeTable.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace flume {

//! A cell of an inflow and the part of the inflow it takes.
struct InflowCell
{
    int i = 0;
    int j = 0;
    double fraction = 0.0;  //!< in (0, 1]; the fractions of an inflow sum to 1
};

/**
 * The cells of the line from (x0, y0) to (x1, y1), m, on the grid of `geom`,
 * each with the part of the line's length inside it: the line is cut where
 * it crosses the lines of the grid, and each piece goes to the cell that
 * holds its middle (a piece along a cell face, to rounding, goes to the
 * cell above or to the right of it, or to the last cell at the high edge of
 * the domain, as a point on that face would). A piece shorter than 1e-12 of
 * the line, rounding where it passes through a corner of the grid, is
 * dropped. The fractions are divided by their sum, so they sum to 1 to round-off. In
 * double precision; the ends must be inside the domain.
 */
std::vector<InflowCell> line_cells (double x0, double y0, double x1, double y1,
                                    const amrex::Geometry& geom);

/**
 * An inflow inside the domain (flume.inflows.<name>.*,
 * Docs/sphinx_doc/theory/SourceTerms.rst): the discharge Q(t), m3/s, of a
 * hydrograph file, added as depth to the cell that holds a point, or to the
 * cells a line crosses, each taking the part f_c of the line's length inside
 * it (line_cells). Over a step from t0 to t1 cell c gains the depth
 *
 *   dh_c = f_c (I(t1) - I(t0)) / (dx dy),
 *
 * with I the integral of the hydrograph (TimeTable), so the volume added over
 * the run is the integral of the hydrograph over the run's time, to
 * round-off. The water enters with no momentum, as rain does: the discharge
 * of the cell is not changed (a manhole that overflows, a pipe that ends
 * upward). The depth is added with a compensated sum, as the rain's.
 *
 * The source acts after the rain and before the infiltration (make_sources),
 * so the soil and the drains see its water in the same step. Its column in
 * the water balance is inflow_<name>.
 */
class InflowSource : public Source
{
public:
    //! Read the hydrograph (collective) and find the cells. Aborts, naming
    //! the inflow's location or line key, when a cell of the inflow is a
    //! building cell (`buildings` above 0; may be null).
    InflowSource (const InflowSpec& spec, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings);

    [[nodiscard]] std::string name () const override { return "inflow_" + m_name; }

    //! Add the inflow from ctx.time to ctx.time_end to the depth of its cells.
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

    [[nodiscard]] std::string balance_name () const override { return "inflow_" + m_name; }

    //! The water added to this rank's cells since the start, m3, less the
    //! part the compensated sum still holds back, as for the rain.
    [[nodiscard]] double balance_volume () const override;

    //! The largest fraction of a cell times the largest discharge from t0 to
    //! t1, over the cell area: the fastest the inflow raises a depth, m/s.
    [[nodiscard]] double max_depth_rate (double t0, double t1) const override;

    //! The cells of the inflow, the rounding carried in the depth and this
    //! rank's volume.
    void write_checkpoint (CheckpointWriter& chk) const override;
    void read_checkpoint (const CheckpointReader& chk) override;

    //! The cells of the inflow, on every rank.
    [[nodiscard]] const std::vector<InflowCell>& cells () const { return m_cells; }

private:
    std::string m_name;
    TimeTable m_q;                    //!< discharge, m3/s
    std::vector<InflowCell> m_cells;  //!< every cell of the inflow
    double m_fraction_max = 0.0;      //!< the largest fraction of one cell
    //! The cells of each local box (index of the MFIter, no tiling), on the device
    std::vector<amrex::Gpu::DeviceVector<InflowCell>> m_local;
    amrex::MultiFab m_carry;          //!< depth rounded away so far, m
    double m_cell_area = 0.0;         //!< m2
    double m_volume = 0.0;            //!< water added to this rank's cells, m3
    double m_volume_carry = 0.0;      //!< rounding carried in the sum m_volume
    //! name:point:i:j or name:line:x0:y0:x1:y1, in the checkpoint; a restart
    //! with another inflow stops.
    std::vector<std::string> m_words;
};

} // namespace flume

#endif
