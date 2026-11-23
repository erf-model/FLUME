#include "FLUME_Forcing.H"
#include "FLUME_RainField.H"
#include "FLUME_TimeTable.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>
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
    return nullptr;
}

} // namespace flume
