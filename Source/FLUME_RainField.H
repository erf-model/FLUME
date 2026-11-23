#ifndef FLUME_RAIN_FIELD_H_
#define FLUME_RAIN_FIELD_H_

#include "FLUME_Forcing.H"
#include "FLUME_Inputs.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace flume {

/**
 * The weights of the entries of a table of times T_0 < ... < T_N in the
 * integral of a quantity over [t0, t1] (s, t1 >= t0): the integral of a
 * quantity whose value at T_k is M_k is the sum of w_k M_k. With
 * RainInterpolation::step the quantity is M_k from T_k to T_{k+1}; with
 * RainInterpolation::linear it is linear in time between the entries. Before
 * T_0 it is M_0, and from T_N on it is M_N, in both modes. Each weight is the
 * difference of a function of the time at the two ends, so the weights of
 * consecutive intervals add up to the weight of their union to round-off.
 * Returns the entries with a weight above 0, in order, as (k, w_k in s).
 */
std::vector<std::pair<int, double>> entry_weights (const std::vector<double>& time,
                                                   RainInterpolation mode, double t0, double t1);

/**
 * Read an index file of gridded rain (flume.rain.file): lines of a time in s
 * and the path of an ESRI ASCII grid of the rain rate in mm/h, with the
 * times strictly increasing; `#` starts a comment that runs to the end of
 * the line, and blank lines are skipped. A relative path is taken from the
 * directory that holds the index file. The I/O rank reads the file and sends
 * the table to every rank. Every rank aborts with the same message, naming
 * the input, the file and the line, when the file cannot be opened or is not
 * a regular file, a line does not hold exactly two items, a time is not a
 * finite number or not above the one before, or the file holds no entry.
 * The grids themselves are not opened here.
 */
void read_rain_index (const std::string& path, std::vector<double>& time,
                      std::vector<std::string>& grid);

/**
 * Rain from a time series of grids of the rain rate (flume.rain.type =
 * gridded), each on its own grid, which need not be the solver's or the same
 * for every entry. Each grid is mapped to the solver's cells by the mean over
 * each cell (sample_conservative), so the rain on the domain is the integral
 * of the grid over the domain; in time the rate is held or interpolated
 * between entries as by entry_weights. The depth over an interval is the
 * sum of the weights times the mapped grids, so the rain a run receives is
 * the integral of the input over the domain and over the run's time, to
 * round-off.
 *
 * At start-up every grid is read and mapped once: each must hold the domain,
 * hold no NODATA value that a cell reads and no value below 0, and its
 * largest mapped rate is kept for max_rain_rate. During the run a grid is
 * read and mapped again when a step needs it. The grids a later step can use
 * (one with step interpolation, two with linear) are kept; a grid that only
 * the current step uses is mapped, added and dropped, so a step that spans
 * many entries holds one extra grid at a time. A grid of no rain is not read
 * again.
 */
class GriddedRain : public Forcing
{
public:
    GriddedRain (const std::string& index, RainInterpolation mode, const amrex::Geometry& geom,
                 const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    [[nodiscard]] std::string name () const override { return "gridded rain"; }
    void rain_depth (amrex::MultiFab& dst, double t0, double t1) const override;
    [[nodiscard]] double max_rain_rate (double t0, double t1) const override;

    //! Number of grids in the index.
    [[nodiscard]] int size () const { return static_cast<int>(m_time.size()); }

private:
    //! Read entry k's grid and map it to the solver's cells, in m/s.
    void load (int k, amrex::MultiFab& dst) const;
    //! Entry k's mapped grid, read and mapped when not held.
    const amrex::MultiFab& frame (int k) const;

    std::string m_index;
    RainInterpolation m_mode = RainInterpolation::step;
    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    std::vector<double> m_time;
    std::vector<std::string> m_grid;
    std::vector<double> m_max;   //!< largest mapped rate of each grid, m/s
    mutable std::map<int, amrex::MultiFab> m_held;
};

/**
 * A rain cell whose rate is a Gaussian of the distance from a centre that
 * moves at a constant velocity (flume.rain.type = moving_cell):
 *
 *   r(x, y, t) = P exp(-((x - x_c(t))^2 + (y - y_c(t))^2) / (2 sigma^2)),
 *   x_c(t) = x_0 + u t,  y_c(t) = y_0 + v t,
 *
 * with P the peak rate. The mean over a cell at one time is exact (a product
 * of differences of error functions); the integral over an interval of time
 * is the 4-point Gauss-Legendre rule on sub-intervals over which the centre
 * moves at most sigma / 4 (one sub-interval for a cell at rest, which is then
 * exact). For tests of the mapping and of domains that cover part of a rain
 * field.
 */
class MovingCellRain : public Forcing
{
public:
    //! `peak` in m/s (>= 0), `sigma` in m (> 0), `center` (m) and `velocity`
    //! (m/s) finite.
    MovingCellRain (double peak, double sigma, amrex::Array<double, 2> center,
                    amrex::Array<double, 2> velocity, const amrex::Geometry& geom);

    [[nodiscard]] std::string name () const override { return "moving rain cell"; }
    void rain_depth (amrex::MultiFab& dst, double t0, double t1) const override;
    [[nodiscard]] double max_rain_rate (double t0, double t1) const override;

private:
    double m_peak;
    double m_sigma;
    amrex::Array<double, 2> m_center;
    amrex::Array<double, 2> m_velocity;
    amrex::Geometry m_geom;
};

} // namespace flume

#endif
