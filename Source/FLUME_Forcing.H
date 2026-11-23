#ifndef FLUME_FORCING_H_
#define FLUME_FORCING_H_

#include "FLUME_Inputs.H"
#include "FLUME_TimeTable.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

namespace flume {

/**
 * What drives the flood from outside the water: today the rain.
 *
 * The solver asks a provider for the rain that falls on each cell of its grid
 * over an interval of time, and for the largest rain rate over an interval
 * (which bounds the time step, FLUME_TimeStep.H). A provider gives the depth
 * of rain over the interval, the integral of the rate, not the rate at one
 * time: the water the run receives is then the integral of the forcing to
 * round-off, whatever the steps are.
 *
 * Rain the same in every cell (UniformRain), rain from a time series of
 * grids and a moving rain cell (FLUME_RainField.H) are providers today. A
 * provider for ERF's microphysics, and the wind and pressure that a storm
 * surge needs, are added behind this interface later (flume.rain.type
 * reserves the name microphysics).
 */
class Forcing
{
public:
    Forcing () = default;
    virtual ~Forcing () = default;
    Forcing (const Forcing&) = delete;
    Forcing& operator= (const Forcing&) = delete;
    Forcing (Forcing&&) = delete;
    Forcing& operator= (Forcing&&) = delete;

    //! The name printed at start-up, for example "constant rain".
    [[nodiscard]] virtual std::string name () const = 0;

    //! Fill component 0 of the valid cells of `dst` with the depth of rain,
    //! m, that falls on each cell from `t0` to `t1` (s, t1 >= t0).
    virtual void rain_depth (amrex::MultiFab& dst, double t0, double t1) const = 0;

    //! The largest rain rate, m/s, over every cell and every time from `t0`
    //! to `t1` (s; t1 may be +infinity). An upper bound is enough.
    [[nodiscard]] virtual double max_rain_rate (double t0, double t1) const = 0;
};

/**
 * Rain that is the same in every cell, with a rate that is a function of
 * time given by a table: times t_0 < t_1 < ... < t_N (s) and rates r_0 ... r_N
 * (m/s, >= 0). With RainInterpolation::step the rate is r_k from t_k to
 * t_{k+1}; with RainInterpolation::linear it is linear in time between the
 * entries. Before t_0 the rate is r_0, and from t_N on it is r_N, in both
 * modes. A table of one entry is constant rain (flume.rain.type = constant).
 *
 * The depth over an interval is I(t1) - I(t0), with I the integral of the
 * rate from t_0, exact for the piecewise-constant or piecewise-linear rate
 * (TimeTable, FLUME_TimeTable.H).
 */
class UniformRain : public Forcing
{
public:
    //! `time` strictly increasing, `rate` (m/s) >= 0 and finite, the same
    //! size, at least one entry. `label` is printed by name().
    UniformRain (std::vector<double> time, std::vector<double> rate, RainInterpolation mode,
                 std::string label);

    [[nodiscard]] std::string name () const override { return m_label; }
    void rain_depth (amrex::MultiFab& dst, double t0, double t1) const override;
    [[nodiscard]] double max_rain_rate (double t0, double t1) const override;

    //! The rate at time t, m/s.
    [[nodiscard]] double rate (double t) const { return m_table.value(t); }
    //! The integral of the rate from t_0 to t, m (negative before t_0).
    [[nodiscard]] double integral (double t) const { return m_table.integral(t); }

private:
    TimeTable m_table;  //!< the rate, m/s
    std::string m_label;
};

//! Rain rates in the inputs and in a hyetograph file are in mm/h; the solver
//! works in m/s.
inline constexpr double mm_per_h = 1.0 / 3.6e6;

/**
 * Read a hyetograph file (flume.rain.file; read_time_table of
 * FLUME_TimeTable.H): lines of two numbers, a time in
 * s and a rain rate in mm/h, with the times strictly increasing; `#` starts a
 * comment that runs to the end of the line, and blank lines are skipped. The
 * I/O rank reads the file and sends the table to every rank. Every rank
 * aborts with the same message, naming the input, the file and the line,
 * when the file cannot be opened or is not a regular file, a line does not
 * hold exactly two finite numbers, a rate is below 0, a time is not above the
 * one before, or the file holds no entry. Returns the times (s) and the rates
 * (mm/h, as written).
 */
void read_hyetograph (const std::string& path, std::vector<double>& time,
                      std::vector<double>& rate_mm_per_h);

//! The provider the inputs ask for (flume.rain.*), or nullptr with
//! flume.rain.type = none, for the cells of `geom` in the boxes `ba` on the
//! ranks `dm`. Prints what it is at start-up.
std::unique_ptr<Forcing> make_forcing (const Inputs& in, const amrex::Geometry& geom,
                                       const amrex::BoxArray& ba,
                                       const amrex::DistributionMapping& dm);

} // namespace flume

#endif
