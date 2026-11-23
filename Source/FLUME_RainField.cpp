#include "FLUME_RainField.H"
#include "FLUME_Raster.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Math.H>
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

//! Parse the index on the I/O rank. Returns the message of the first failed
//! check, or an empty string.
std::string
parse_index (const std::string& path, std::vector<double>& time, std::vector<std::string>& grid)
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
    const std::filesystem::path dir = std::filesystem::path(path).parent_path();
    std::string line;
    int line_no = 0;
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
            return where + ": a line holds a time in s and the path of an ESRI ASCII grid of the"
                   " rain rate in mm/h; it holds " + std::to_string(tok.size()) + " item(s)";
        }
        double tv = 0.0;
        if (!to_double(tok[0], tv)) {
            return where + ": the time must be a finite number; it is \"" + tok[0] + "\"";
        }
        if (!time.empty() && !(tv > time.back())) {
            return where + ": the time must be above the time of the entry before ("
                   + all_digits(time.back()) + " s); it is " + all_digits(tv);
        }
        const std::filesystem::path g(tok[1]);
        time.push_back(tv);
        grid.push_back(g.is_relative() ? (dir / g).string() : g.string());
    }
    if (f.bad()) { return w + ": the file cannot be read"; }
    if (time.empty()) {
        return w + ": the file holds no entry: it needs at least one line with a time in s and"
               " the path of an ESRI ASCII grid of the rain rate in mm/h";
    }
    return {};
}

//! Send a string from the I/O rank to every rank.
void
bcast_string (std::string& s, int io)
{
    long long len = static_cast<long long>(s.size());
    amrex::ParallelDescriptor::Bcast(&len, 1, io);
    s.resize(static_cast<std::size_t>(len));
    if (len > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(len), io); }
}

//! The integrals from T_j to t, s, of the weights of the two entries of
//! linear interpolation on [T_j, T_j + d]: the left entry's weight is 1 - u
//! and the right entry's u, with u = (t - T_j) / d held in [0, 1].
struct Pair
{
    double left;
    double right;
};

Pair
linear_part (double tj, double d, double t)
{
    const double u = std::clamp((t - tj) / d, 0.0, 1.0);
    return {d * (u - 0.5 * u * u), d * 0.5 * u * u};
}

//! The integral from T_0 to t (negative before T_0) of the weight of entry k
//! in the interpolation of the table.
double
cumulative_weight (const std::vector<double>& T, RainInterpolation mode, int k, double t)
{
    const int n = static_cast<int>(T.size());
    const auto ku = static_cast<std::size_t>(k);
    if (n == 1) { return t - T[0]; }
    if (mode == RainInterpolation::step) {
        if (k == 0) { return std::min(t, T[1]) - T[0]; }
        if (k == n - 1) { return std::max(t, T[ku]) - T[ku]; }
        return std::clamp(t, T[ku], T[ku + 1]) - T[ku];
    }
    double c = 0.0;
    if (k == 0) { c += std::min(t, T[0]) - T[0]; }
    if (k == n - 1) { c += std::max(t, T[ku]) - T[ku]; }
    if (k > 0) { c += linear_part(T[ku - 1], T[ku] - T[ku - 1], t).right; }
    if (k < n - 1) { c += linear_part(T[ku], T[ku + 1] - T[ku], t).left; }
    return c;
}

//! The first and last entries whose weight can be above 0 over [t0, t1].
std::pair<int, int>
entry_range (const std::vector<double>& T, RainInterpolation mode, double t0, double t1)
{
    const int n = static_cast<int>(T.size());
    // The last entry at or before t, or 0 before T_0.
    const auto at = [&] (double t) {
        const auto k = static_cast<int>(std::upper_bound(T.begin(), T.end(), t) - T.begin()) - 1;
        return std::max(k, 0);
    };
    const int lo = at(t0);
    const int hi = (mode == RainInterpolation::linear) ? std::min(at(t1) + 1, n - 1) : at(t1);
    return {lo, hi};
}

//! erf(b) - erf(a) for a <= b, from erfc on the positive side so that a
//! cell far from the centre keeps its digits.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
erf_difference (amrex::Real a, amrex::Real b) noexcept
{
    if (a >= amrex::Real(0.0)) { return std::erfc(a) - std::erfc(b); }
    if (b <= amrex::Real(0.0)) { return std::erfc(-b) - std::erfc(-a); }
    return std::erf(b) - std::erf(a);
}

} // namespace

std::vector<std::pair<int, double>>
entry_weights (const std::vector<double>& time, RainInterpolation mode, double t0, double t1)
{
    std::vector<std::pair<int, double>> w;
    const auto [lo, hi] = entry_range(time, mode, t0, t1);
    for (int k = lo; k <= hi; ++k) {
        // Each cumulative weight grows with the time; the difference can
        // only fall below 0 by rounding.
        const double wk = cumulative_weight(time, mode, k, t1) - cumulative_weight(time, mode, k, t0);
        if (wk > 0.0) { w.emplace_back(k, wk); }
    }
    return w;
}

void
read_rain_index (const std::string& path, std::vector<double>& time,
                 std::vector<std::string>& grid)
{
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    time.clear();
    grid.clear();
    std::string error;
    if (amrex::ParallelDescriptor::IOProcessor()) { error = parse_index(path, time, grid); }
    bcast_string(error, io);
    if (!error.empty()) { amrex::Abort(error); }
    long long n = static_cast<long long>(time.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, io);
    time.resize(static_cast<std::size_t>(n));
    amrex::ParallelDescriptor::Bcast(time.data(), static_cast<std::size_t>(n), io);
    // The paths, one per line (a path holds no line break).
    std::string joined;
    for (const auto& g : grid) { joined += g + "\n"; }
    bcast_string(joined, io);
    grid.clear();
    std::istringstream in(joined);
    std::string g;
    while (std::getline(in, g)) { grid.push_back(g); }
    AMREX_ALWAYS_ASSERT(static_cast<long long>(grid.size()) == n);
}

GriddedRain::GriddedRain (const std::string& index, RainInterpolation mode,
                          const amrex::Geometry& geom, const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm)
    : m_index(index), m_mode(mode), m_geom(geom), m_ba(ba), m_dm(dm)
{
    read_rain_index(index, m_time, m_grid);
    // Every grid once: the checks, and the largest rate for the time step.
    amrex::MultiFab tmp(m_ba, m_dm, 1, 0);
    m_max.resize(m_time.size());
    for (int k = 0; k < size(); ++k) {
        load(k, tmp);
        m_max[static_cast<std::size_t>(k)] = static_cast<double>(tmp.max(0));
    }
}

void
GriddedRain::load (int k, amrex::MultiFab& dst) const
{
    const auto ku = static_cast<std::size_t>(k);
    const std::string key = "rain.file (" + m_index + "), entry " + std::to_string(k + 1)
                            + " at " + all_digits(m_time[ku]) + " s: grid";
    const Raster r = read_raster(key, m_grid[ku], m_geom);
    // A rain rate below 0 would take water away. Every rank holds the same
    // values, so every rank aborts with the same message.
    for (std::size_t m = 0; m < r.value.size(); ++m) {
        if (r.nodata[m] == 0 && r.value[m] < 0.0) {
            const auto i = static_cast<long long>(m % static_cast<std::size_t>(r.nx));
            const auto j = static_cast<long long>(m / static_cast<std::size_t>(r.nx));
            amrex::Abort("flume." + key + " (" + m_grid[ku] + "): the rain rate must be >= 0 mm/h;"
                         " it is " + all_digits(r.value[m]) + " in the cell whose west edge is at x = "
                         + all_digits(r.x0 + static_cast<double>(i) * r.dx)
                         + " and south edge at y = "
                         + all_digits(r.y0 + static_cast<double>(j) * r.dy));
        }
    }
    sample_conservative(dst, 0, r, m_geom);
    dst.mult(static_cast<amrex::Real>(mm_per_h), 0, 1, 0);
}

const amrex::MultiFab&
GriddedRain::frame (int k) const
{
    auto it = m_held.find(k);
    if (it == m_held.end()) {
        it = m_held.emplace(std::piecewise_construct, std::forward_as_tuple(k),
                            std::forward_as_tuple(m_ba, m_dm, 1, 0)).first;
        load(k, it->second);
    }
    return it->second;
}

void
GriddedRain::rain_depth (amrex::MultiFab& dst, double t0, double t1) const
{
    const auto w = entry_weights(m_time, m_mode, t0, t1);
    // The first grid an interval from t1 on can use: the run's time only
    // moves forward, so grids before it are not needed again. A grid this
    // interval uses before it is mapped into `once` and not kept, so a step
    // that spans many entries holds one such grid at a time.
    const int keep = entry_range(m_time, m_mode, t1, t1).first;
    std::unique_ptr<amrex::MultiFab> once;
    dst.setVal(amrex::Real(0.0), 0, 1, 0);
    for (const auto& [k, wk] : w) {
        // A grid of no rain adds nothing (its rates are >= 0).
        if (m_max[static_cast<std::size_t>(k)] == 0.0) { continue; }
        const amrex::MultiFab* mf = nullptr;
        if (k >= keep || m_held.count(k) != 0) {
            mf = &frame(k);
        } else {
            if (!once) { once = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0); }
            load(k, *once);
            mf = once.get();
        }
        const auto& src = mf->const_arrays();
        const auto& d = dst.arrays();
        const auto a = static_cast<amrex::Real>(wk);
        amrex::ParallelFor(dst,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int kk) noexcept
            {
                d[box_no](i, j, kk) += a * src[box_no](i, j, kk);
            });
        // `once` is overwritten by the next grid: finish with it first.
        amrex::Gpu::streamSynchronize();
    }
    for (auto it = m_held.begin(); it != m_held.end() && it->first < keep;) {
        it = m_held.erase(it);
    }
}

double
GriddedRain::max_rain_rate (double t0, double t1) const
{
    const auto [lo, hi] = entry_range(m_time, m_mode, t0, t1);
    double r = 0.0;
    for (int k = lo; k <= hi; ++k) { r = std::max(r, m_max[static_cast<std::size_t>(k)]); }
    return r;
}

MovingCellRain::MovingCellRain (double peak, double sigma, amrex::Array<double, 2> center,
                                amrex::Array<double, 2> velocity, const amrex::Geometry& geom)
    : m_peak(peak), m_sigma(sigma), m_center(center), m_velocity(velocity), m_geom(geom)
{
    AMREX_ALWAYS_ASSERT(m_peak >= 0.0 && m_sigma > 0.0);
}

void
MovingCellRain::rain_depth (amrex::MultiFab& dst, double t0, double t1) const
{
    // No rain, or no time: nothing to integrate (and no sub-intervals to
    // count over a long step of a dry domain under a cell of no rain).
    if (m_peak == 0.0 || !(t1 > t0)) {
        dst.setVal(amrex::Real(0.0), 0, 1, 0);
        return;
    }
    // 4-point Gauss-Legendre nodes and weights on [-1, 1].
    constexpr double gx[4] = {-0.86113631159405257522, -0.33998104358485626480,
                              0.33998104358485626480, 0.86113631159405257522};
    constexpr double gw[4] = {0.34785484513745385737, 0.65214515486254614263,
                              0.65214515486254614263, 0.34785484513745385737};
    const double speed = std::hypot(m_velocity[0], m_velocity[1]);
    const double span = t1 - t0;
    // At most sigma / 4 of travel per sub-interval. A step over which the
    // cell crosses more than 2500 sigma would cost every cell more than 4e4
    // error-function pairs: stop with a message rather than run that long.
    const double want = std::ceil(speed * span / (0.25 * m_sigma));
    constexpr double max_sub = 1.0e4;
    if (want > max_sub) {
        amrex::Abort("flume.rain.cell: the cell moves " + all_digits(speed * span)
                     + " m in one step of " + all_digits(span) + " s, more than 2500 times"
                     " flume.rain.cell.sigma (" + all_digits(m_sigma) + " m); shorten the step"
                     " (flume.cfl, or flume.dry_depth on a dry domain) or widen the cell");
    }
    const auto n_sub = static_cast<long long>(std::max(1.0, want));
    // Per node: the centre relative to the domain's lower corner (m) and
    // the weight times the peak rate (m). Positions in double, then in the
    // build's precision: relative to the corner they lose no digits that
    // matter at the scale of a cell.
    const auto n_node = static_cast<std::size_t>(4 * n_sub);
    std::vector<amrex::Real> node(3 * n_node);
    const double h = span / static_cast<double>(n_sub);
    const double lo_x = static_cast<double>(m_geom.ProbLo(0));
    const double lo_y = static_cast<double>(m_geom.ProbLo(1));
    for (long long s = 0; s < n_sub; ++s) {
        const double mid = t0 + (static_cast<double>(s) + 0.5) * h;
        for (int q = 0; q < 4; ++q) {
            const double t = mid + 0.5 * h * gx[q];
            const auto m = 3 * (static_cast<std::size_t>(4 * s) + static_cast<std::size_t>(q));
            node[m] = static_cast<amrex::Real>(m_center[0] + m_velocity[0] * t - lo_x);
            node[m + 1] = static_cast<amrex::Real>(m_center[1] + m_velocity[1] * t - lo_y);
            node[m + 2] = static_cast<amrex::Real>(0.5 * h * gw[q] * m_peak);
        }
    }
    amrex::Gpu::DeviceVector<amrex::Real> dnode(node.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, node.begin(), node.end(), dnode.begin());
    const amrex::Real* nodes = dnode.data();
    const auto nn = static_cast<int>(n_node);

    const amrex::Real gdx = m_geom.CellSize(0);
    const amrex::Real gdy = m_geom.CellSize(1);
    const int ilo = m_geom.Domain().smallEnd(0);
    const int jlo = m_geom.Domain().smallEnd(1);
    const auto inv = static_cast<amrex::Real>(1.0 / (std::sqrt(2.0) * m_sigma));
    // The mean over [a, b] of exp(-(x - c)^2 / (2 sigma^2)) is
    // sigma sqrt(pi / 2) (erf((b - c) / (sigma sqrt 2)) - erf((a - c) / (sigma sqrt 2))) / (b - a).
    const auto fx = static_cast<amrex::Real>(m_sigma * std::sqrt(0.5 * amrex::Math::pi<double>())
                                             / static_cast<double>(gdx));
    const auto fy = static_cast<amrex::Real>(m_sigma * std::sqrt(0.5 * amrex::Math::pi<double>())
                                             / static_cast<double>(gdy));
    const auto& d = dst.arrays();
    amrex::ParallelFor(dst,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real xa = static_cast<amrex::Real>(i - ilo) * gdx;
            const amrex::Real xb = static_cast<amrex::Real>(i + 1 - ilo) * gdx;
            const amrex::Real ya = static_cast<amrex::Real>(j - jlo) * gdy;
            const amrex::Real yb = static_cast<amrex::Real>(j + 1 - jlo) * gdy;
            amrex::Real sum = amrex::Real(0.0);
            for (int q = 0; q < nn; ++q) {
                const amrex::Real cx = nodes[3 * q];
                const amrex::Real cy = nodes[3 * q + 1];
                const amrex::Real ex = fx * erf_difference((xa - cx) * inv, (xb - cx) * inv);
                const amrex::Real ey = fy * erf_difference((ya - cy) * inv, (yb - cy) * inv);
                sum += nodes[3 * q + 2] * ex * ey;
            }
            d[box_no](i, j, k) = amrex::max(sum, amrex::Real(0.0));
        });
    amrex::Gpu::streamSynchronize();
}

double
MovingCellRain::max_rain_rate (double /*t0*/, double /*t1*/) const
{
    return m_peak;
}

} // namespace flume
