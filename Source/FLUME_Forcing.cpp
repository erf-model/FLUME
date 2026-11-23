#include "FLUME_Forcing.H"
#include "FLUME_RainField.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>

namespace flume {

namespace {

//! The whole token as a finite double, or false.
bool
to_double (const std::string& t, double& v)
{
    if (t.empty()) { return false; }
    errno = 0;
    char* end = nullptr;
    v = std::strtod(t.c_str(), &end);
    return end == t.c_str() + t.size() && errno != ERANGE && std::isfinite(v);
}

//! A double as text with every digit.
std::string
all_digits (double v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
    return s.str();
}

//! Parse the file on the I/O rank. Returns the message of the first failed
//! check, or an empty string.
std::string
parse_hyetograph (const std::string& path, std::vector<double>& time, std::vector<double>& rate)
{
    const std::string w = "flume.rain.file (" + path + ")";
    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) {
        std::string msg = w + ": the file cannot be opened";
        if (std::filesystem::path(path).is_relative()) {
            const auto cwd = std::filesystem::current_path(ec);
            if (!ec) { msg += " (a relative path, taken from " + cwd.string() + ")"; }
        }
        return msg;
    }
    if (!std::filesystem::is_regular_file(path, ec)) { return w + ": the path is not a regular file"; }
    std::ifstream f(path);
    if (!f.good()) { return w + ": the file cannot be opened"; }
    std::string line;
    long long line_no = 0;
    while (std::getline(f, line)) {
        ++line_no;
        const std::string where = w + ", line " + std::to_string(line_no);
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) { line.erase(hash); }
        std::istringstream in(line);
        std::vector<std::string> tok;
        std::string t;
        while (in >> t) { tok.push_back(t); }
        if (tok.empty()) { continue; }
        if (tok.size() != 2) {
            return where + ": a line holds two numbers, a time in s and a rain rate in mm/h;"
                   " it holds " + std::to_string(tok.size()) + " item(s)";
        }
        double tv = 0.0;
        double rv = 0.0;
        if (!to_double(tok[0], tv)) {
            return where + ": the time must be a finite number; it is \"" + tok[0] + "\"";
        }
        if (!to_double(tok[1], rv)) {
            return where + ": the rain rate must be a finite number; it is \"" + tok[1] + "\"";
        }
        if (!(rv >= 0.0)) {
            return where + ": the rain rate must be >= 0 mm/h; it is " + all_digits(rv);
        }
        if (!time.empty() && !(tv > time.back())) {
            return where + ": the time must be above the time of the entry before ("
                   + all_digits(time.back()) + " s); it is " + all_digits(tv);
        }
        time.push_back(tv);
        rate.push_back(rv);
    }
    if (f.bad()) { return w + ": the file cannot be read"; }
    if (time.empty()) {
        return w + ": the file holds no entry: it needs at least one line with a time in s and a"
               " rain rate in mm/h";
    }
    return {};
}

} // namespace

UniformRain::UniformRain (std::vector<double> time, std::vector<double> rate,
                          RainInterpolation mode, std::string label)
    : m_time(std::move(time)), m_rate(std::move(rate)), m_mode(mode), m_label(std::move(label))
{
    AMREX_ALWAYS_ASSERT(!m_time.empty() && m_time.size() == m_rate.size());
    m_cumulative.assign(m_time.size(), 0.0);
    for (std::size_t k = 1; k < m_time.size(); ++k) {
        const double dt = m_time[k] - m_time[k - 1];
        const double r = (m_mode == RainInterpolation::step) ? m_rate[k - 1]
                                                              : 0.5 * (m_rate[k - 1] + m_rate[k]);
        m_cumulative[k] = m_cumulative[k - 1] + r * dt;
    }
}

double
UniformRain::rate (double t) const
{
    if (t < m_time.front()) { return m_rate.front(); }
    if (t >= m_time.back()) { return m_rate.back(); }
    // The last entry at or before t; t is below the last time, so k + 1 exists.
    const auto k = static_cast<std::size_t>(
        std::upper_bound(m_time.begin(), m_time.end(), t) - m_time.begin()) - 1;
    if (m_mode == RainInterpolation::step) { return m_rate[k]; }
    const double w = (t - m_time[k]) / (m_time[k + 1] - m_time[k]);
    return m_rate[k] + w * (m_rate[k + 1] - m_rate[k]);
}

double
UniformRain::integral (double t) const
{
    if (t < m_time.front()) { return m_rate.front() * (t - m_time.front()); }
    if (t >= m_time.back()) { return m_cumulative.back() + m_rate.back() * (t - m_time.back()); }
    const auto k = static_cast<std::size_t>(
        std::upper_bound(m_time.begin(), m_time.end(), t) - m_time.begin()) - 1;
    const double s = t - m_time[k];
    if (m_mode == RainInterpolation::step) { return m_cumulative[k] + m_rate[k] * s; }
    // The mean of the rates at t_k and at t, times the time since t_k.
    const double slope = (m_rate[k + 1] - m_rate[k]) / (m_time[k + 1] - m_time[k]);
    return m_cumulative[k] + s * (m_rate[k] + 0.5 * slope * s);
}

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
    // The rate is largest at an end of the interval or at an entry inside it.
    double r = std::max(rate(t0), (t1 >= m_time.back()) ? m_rate.back() : rate(t1));
    const auto first = std::upper_bound(m_time.begin(), m_time.end(), t0);
    for (auto it = first; it != m_time.end() && *it < t1; ++it) {
        r = std::max(r, m_rate[static_cast<std::size_t>(it - m_time.begin())]);
    }
    return r;
}

void
read_hyetograph (const std::string& path, std::vector<double>& time,
                 std::vector<double>& rate_mm_per_h)
{
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    time.clear();
    rate_mm_per_h.clear();
    std::string error;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        error = parse_hyetograph(path, time, rate_mm_per_h);
    }
    // The message (empty when the file is good), then the table.
    long long len = static_cast<long long>(error.size());
    amrex::ParallelDescriptor::Bcast(&len, 1, io);
    if (len > 0) {
        error.resize(static_cast<std::size_t>(len));
        amrex::ParallelDescriptor::Bcast(error.data(), static_cast<std::size_t>(len), io);
        amrex::Abort(error);
    }
    long long n = static_cast<long long>(time.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, io);
    time.resize(static_cast<std::size_t>(n));
    rate_mm_per_h.resize(static_cast<std::size_t>(n));
    amrex::ParallelDescriptor::Bcast(time.data(), static_cast<std::size_t>(n), io);
    amrex::ParallelDescriptor::Bcast(rate_mm_per_h.data(), static_cast<std::size_t>(n), io);
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
