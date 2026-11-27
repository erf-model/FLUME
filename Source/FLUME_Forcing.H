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
 * grids and a moving rain cell (FLUME_RainField.H) are providers, and so is
 * the rain a host model hands over step by step (HostRain: ERF's
 * microphysics). The wind and pressure that a storm surge needs are added
 * behind this interface later.
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

/**
 * Rain a host model hands over one host step at a time (flume.rain.type =
 * microphysics: the surface precipitation of ERF's microphysics, set by
 * FLUME_ERFFloodEngine.cpp). Before FLUME takes its steps over a host step
 * from T0 to T1, the host fills depth() with the depth of water, m, that fell
 * on each cell over that step, and calls set_step(T0, T1). The rate is held
 * constant over the host step: a FLUME step from t0 to t1 inside it receives
 *
 *     D * (f(t1) - f(t0)),   f(t) = (t - T0) / (T1 - T0), clamped to [0, 1],
 *
 * with f(T0) = 0 and f(T1) = 1 exactly, so a host step taken in one FLUME step
 * receives D to the last bit, and the steps that cover a host step receive D
 * to round-off. Outside the host step no rain falls.
 *
 * D can be slightly negative: ERF's step precipitation is a difference of
 * accumulators, and a correction from a finer level can take some back. It is
 * passed on as it is; the rain source (add_rain_depth, FLUME_Rain.H) holds a
 * negative depth back until the cell has the water, so the depth never goes
 * below 0 and the rain the run receives is the host's to round-off.
 */
class HostRain : public Forcing
{
public:
    HostRain (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    [[nodiscard]] std::string name () const override { return "rain from the host model"; }
    void rain_depth (amrex::MultiFab& dst, double t0, double t1) const override;
    [[nodiscard]] double max_rain_rate (double t0, double t1) const override;

    //! The depth of water, m, on each cell over the next host step (one
    //! component, no ghost cells), for the host to fill before set_step().
    [[nodiscard]] amrex::MultiFab& depth () { return m_depth; }

    //! The host step that depth() covers, from `t0` to `t1` (s, t1 > t0).
    //! Aborts, naming `host`, when a depth is not finite. Collective.
    void set_step (double t0, double t1, const std::string& host);

private:
    amrex::MultiFab m_depth;  //!< the host step's depth, m
    double m_t0 = 0.0;        //!< start of the host step, s
    double m_t1 = 0.0;        //!< end of the host step, s; m_t1 == m_t0 before set_step()
    double m_max_rate = 0.0;  //!< the largest depth over the step's length, m/s, >= 0
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
