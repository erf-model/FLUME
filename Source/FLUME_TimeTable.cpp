#include "FLUME_TimeTable.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
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
parse_table (const std::string& key, const std::string& path, const std::string& quantity,
             const std::string& unit, std::vector<double>& time, std::vector<double>& value)
{
    const std::string w = key + " (" + path + ")";
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
            return where + ": a line holds two numbers, a time in s and a " + quantity + " in "
                   + unit + "; it holds " + std::to_string(tok.size()) + " item(s)";
        }
        double tv = 0.0;
        double vv = 0.0;
        if (!to_double(tok[0], tv)) {
            return where + ": the time must be a finite number; it is \"" + tok[0] + "\"";
        }
        if (!to_double(tok[1], vv)) {
            return where + ": the " + quantity + " must be a finite number; it is \"" + tok[1]
                   + "\"";
        }
        if (!(vv >= 0.0)) {
            return where + ": the " + quantity + " must be >= 0 " + unit + "; it is "
                   + all_digits(vv);
        }
        if (!time.empty() && !(tv > time.back())) {
            return where + ": the time must be above the time of the entry before ("
                   + all_digits(time.back()) + " s); it is " + all_digits(tv);
        }
        time.push_back(tv);
        value.push_back(vv);
    }
    if (f.bad()) { return w + ": the file cannot be read"; }
    if (time.empty()) {
        return w + ": the file holds no entry: it needs at least one line with a time in s and a "
               + quantity + " in " + unit;
    }
    return {};
}

} // namespace

TimeTable::TimeTable (std::vector<double> time, std::vector<double> value, RainInterpolation mode)
    : m_time(std::move(time)), m_value(std::move(value)), m_mode(mode)
{
    AMREX_ALWAYS_ASSERT(!m_time.empty() && m_time.size() == m_value.size());
    m_cumulative.assign(m_time.size(), 0.0);
    for (std::size_t k = 1; k < m_time.size(); ++k) {
        const double dt = m_time[k] - m_time[k - 1];
        const double r = (m_mode == RainInterpolation::step) ? m_value[k - 1]
                                                              : 0.5 * (m_value[k - 1] + m_value[k]);
        m_cumulative[k] = m_cumulative[k - 1] + r * dt;
    }
}

double
TimeTable::value (double t) const
{
    if (t < m_time.front()) { return m_value.front(); }
    if (t >= m_time.back()) { return m_value.back(); }
    // The last entry at or before t; t is below the last time, so k + 1 exists.
    const auto k = static_cast<std::size_t>(
        std::upper_bound(m_time.begin(), m_time.end(), t) - m_time.begin()) - 1;
    if (m_mode == RainInterpolation::step) { return m_value[k]; }
    const double w = (t - m_time[k]) / (m_time[k + 1] - m_time[k]);
    return m_value[k] + w * (m_value[k + 1] - m_value[k]);
}

double
TimeTable::integral (double t) const
{
    if (t < m_time.front()) { return m_value.front() * (t - m_time.front()); }
    if (t >= m_time.back()) { return m_cumulative.back() + m_value.back() * (t - m_time.back()); }
    const auto k = static_cast<std::size_t>(
        std::upper_bound(m_time.begin(), m_time.end(), t) - m_time.begin()) - 1;
    const double s = t - m_time[k];
    if (m_mode == RainInterpolation::step) { return m_cumulative[k] + m_value[k] * s; }
    // The mean of the values at t_k and at t, times the time since t_k.
    const double slope = (m_value[k + 1] - m_value[k]) / (m_time[k + 1] - m_time[k]);
    return m_cumulative[k] + s * (m_value[k] + 0.5 * slope * s);
}

double
TimeTable::max_value (double t0, double t1) const
{
    // The value is largest at an end of the interval or at an entry inside it.
    double r = std::max(value(t0), (t1 >= m_time.back()) ? m_value.back() : value(t1));
    const auto first = std::upper_bound(m_time.begin(), m_time.end(), t0);
    for (auto it = first; it != m_time.end() && *it < t1; ++it) {
        r = std::max(r, m_value[static_cast<std::size_t>(it - m_time.begin())]);
    }
    return r;
}

void
read_time_table (const std::string& key, const std::string& path, const std::string& quantity,
                 const std::string& unit, std::vector<double>& time, std::vector<double>& value)
{
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    time.clear();
    value.clear();
    std::string error;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        error = parse_table(key, path, quantity, unit, time, value);
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
    value.resize(static_cast<std::size_t>(n));
    amrex::ParallelDescriptor::Bcast(time.data(), static_cast<std::size_t>(n), io);
    amrex::ParallelDescriptor::Bcast(value.data(), static_cast<std::size_t>(n), io);
}

} // namespace flume
