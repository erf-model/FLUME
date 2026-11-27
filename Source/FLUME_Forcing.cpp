#include "FLUME_Forcing.H"
#include "FLUME_Init.H"
#include "FLUME_RainField.H"
#include "FLUME_TimeTable.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace flume {

UniformRain::UniformRain (std::vector<double> time, std::vector<double> rate,
                          RainInterpolation mode, std::string label)
    : m_table(std::move(time), std::move(rate), mode), m_label(std::move(label))
{}

void
UniformRain::rain_depth (amrex::MultiFab& dst, double t0, double t1) const
{
    // Never below 0: the integral can only grow, but its two values carry
    // their own rounding.
    const double d = std::max(integral(t1) - integral(t0), 0.0);
    dst.setVal(static_cast<amrex::Real>(d), 0, 1, 0);
}

double
UniformRain::max_rain_rate (double t0, double t1) const
{
    return m_table.max_value(t0, t1);
}

HostRain::HostRain (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : m_depth(ba, dm, 1, 0)
{
    m_depth.setVal(amrex::Real(0.0));
}

void
HostRain::set_step (double t0, double t1, const std::string& host)
{
    // A host step of length 0 has no rate to spread its rain over; ERF never
    // advances the flood by 0 s.
    if (!(std::isfinite(t0) && std::isfinite(t1) && t1 > t0)) {
        std::ostringstream msg;
        msg << std::setprecision(17) << "FLUME: the step of " << host << " from " << t0
            << " s to " << t1 << " s is empty or not finite";
        amrex::Abort(msg.str());
    }
    const amrex::Long bad = count_not_finite(m_depth);
    if (bad > 0) {
        amrex::Abort("FLUME: the rain from " + host + " is not finite in " + std::to_string(bad)
                     + " cell(s)");
    }
    m_t0 = t0;
    m_t1 = t1;
    // max() is collective; a step of no rain, or of a correction only, has a
    // rate of 0.
    const auto d_max = static_cast<double>(m_depth.max(0));
    m_max_rate = std::max(d_max, 0.0) / (t1 - t0);
}

void
HostRain::rain_depth (amrex::MultiFab& dst, double t0, double t1) const
{
    // The fraction of the host step from its start to t, exactly 0 and 1 at
    // its ends.
    const auto f = [this] (double t) {
        if (t <= m_t0) { return 0.0; }
        if (t >= m_t1) { return 1.0; }
        return (t - m_t0) / (m_t1 - m_t0);
    };
    const double w = (m_t1 > m_t0) ? f(t1) - f(t0) : 0.0;
    amrex::MultiFab::Copy(dst, m_depth, 0, 0, 1, 0);
    if (w != 1.0) { dst.mult(static_cast<amrex::Real>(w), 0, 1, 0); }
}

double
HostRain::max_rain_rate (double t0, double t1) const
{
    return (m_t1 > m_t0 && t0 < m_t1 && t1 > m_t0) ? m_max_rate : 0.0;
}

void
read_hyetograph (const std::string& path, std::vector<double>& time,
                 std::vector<double>& rate_mm_per_h)
{
    read_time_table("flume.rain.file", path, "rain rate", "mm/h", time, rate_mm_per_h);
}

std::unique_ptr<Forcing>
make_forcing (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
              const amrex::DistributionMapping& dm)
{
    if (in.rain_type == RainType::constant) {
        amrex::Print() << "FLUME: rain: constant, " << in.rain_rate << " mm/h\n";
        return std::make_unique<UniformRain>(std::vector<double>{0.0},
                                             std::vector<double>{in.rain_rate * mm_per_h},
                                             RainInterpolation::step, "constant rain");
    }
    if (in.rain_type == RainType::hyetograph) {
        std::vector<double> time;
        std::vector<double> rate;
        read_hyetograph(in.rain_file, time, rate);
        const double r_max = *std::max_element(rate.begin(), rate.end());
        if (time.front() > 0.0 && rate.front() > 0.0) {
            amrex::Print() << "FLUME: rain: note: the first time of " << in.rain_file << " is "
                           << time.front() << " s; its rate, " << rate.front()
                           << " mm/h, falls from the start of the run (a line \"0 0\" first"
                              " starts the storm at " << time.front() << " s)\n";
        }
        for (auto& r : rate) { r *= mm_per_h; }
        const bool linear = in.rain_interpolation == RainInterpolation::linear;
        auto rain = std::make_unique<UniformRain>(std::move(time), std::move(rate),
                                                  in.rain_interpolation, "hyetograph rain");
        amrex::Print() << "FLUME: rain: hyetograph from " << in.rain_file << ", "
                       << (linear ? "linear" : "step") << " in time, up to " << r_max << " mm/h\n";
        return rain;
    }
    if (in.rain_type == RainType::gridded) {
        auto rain = std::make_unique<GriddedRain>(in.rain_file, in.rain_interpolation, geom, ba, dm);
        const bool linear = in.rain_interpolation == RainInterpolation::linear;
        amrex::Print() << "FLUME: rain: gridded from " << in.rain_file << ", " << rain->size()
                       << " grid(s), " << (linear ? "linear" : "step") << " in time, up to "
                       << rain->max_rain_rate(-std::numeric_limits<double>::infinity(),
                                              std::numeric_limits<double>::infinity()) / mm_per_h
                       << " mm/h on the solver's cells\n";
        return rain;
    }
    if (in.rain_type == RainType::moving_cell) {
        amrex::Print() << "FLUME: rain: moving cell, peak " << in.rain_cell_peak_rate
                       << " mm/h, sigma " << in.rain_cell_sigma << " m, centre ("
                       << in.rain_cell_center[0] << ", " << in.rain_cell_center[1]
                       << ") m at t = 0, velocity (" << in.rain_cell_velocity[0] << ", "
                       << in.rain_cell_velocity[1] << ") m/s\n";
        return std::make_unique<MovingCellRain>(in.rain_cell_peak_rate * mm_per_h,
                                                in.rain_cell_sigma, in.rain_cell_center,
                                                in.rain_cell_velocity, geom);
    }
    if (in.rain_type == RainType::microphysics) {
        amrex::Print() << "FLUME: rain: the surface precipitation of ERF's microphysics, liquid and"
                          " frozen as liquid water, held constant over each ERF step\n";
        return std::make_unique<HostRain>(ba, dm);
    }
    return nullptr;
}

} // namespace flume
