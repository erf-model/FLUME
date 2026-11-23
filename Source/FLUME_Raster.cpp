#include "FLUME_Raster.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace flume {

namespace {

//! "flume.<key> (<path>)", the start of every message about a file.
std::string
where (const std::string& key, const std::string& path)
{
    return "flume." + key + " (" + path + ")";
}

//! A double as text with every digit.
std::string
all_digits (double v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
    return s.str();
}

//! A failed check on the I/O rank: the message is sent to every rank, which
//! aborts with it.
struct RasterError
{
    std::string message;
};

[[noreturn]] void
fail (const std::string& message)
{
    throw RasterError{message};
}

//! Splits a buffer into whitespace-separated tokens, one at a time.
class Tokens
{
public:
    explicit Tokens (const std::string& buf) : m_buf(buf) {}

    //! The next token, or an empty string at the end.
    std::string next ()
    {
        const std::size_t n = m_buf.size();
        while (m_pos < n && is_space(m_buf[m_pos])) { ++m_pos; }
        const std::size_t start = m_pos;
        while (m_pos < n && !is_space(m_buf[m_pos])) { ++m_pos; }
        return m_buf.substr(start, m_pos - start);
    }

private:
    static bool is_space (char c)
    {
        return c == '\0' || std::isspace(static_cast<unsigned char>(c)) != 0;
    }

    const std::string& m_buf;
    std::size_t m_pos = 0;
};

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

std::string
lower (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [] (unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

//! Read and parse the file, keeping the part the domain reads. Runs on the
//! I/O rank only; every failed check throws RasterError.
Raster
parse (const std::string& key, const std::string& path, const amrex::Geometry& geom)
{
    const std::string w = where(key, path);

    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) { fail(w + ": the file cannot be opened"); }
    if (!std::filesystem::is_regular_file(path, ec)) { fail(w + ": the path is not a regular file"); }
    std::ifstream f(path, std::ios::binary);
    if (!f.good()) { fail(w + ": the file cannot be opened"); }
    const std::string buf((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    if (f.bad()) { fail(w + ": the file cannot be read"); }

    // Header: keyword-value pairs up to the first token that is not a keyword.
    const std::vector<std::string> known{"ncols", "nrows", "xllcorner", "xllcenter", "yllcorner",
                                         "yllcenter", "cellsize", "dx", "dy", "nodata_value"};
    std::vector<int> seen(known.size(), 0);
    std::vector<double> hv(known.size(), 0.0);
    Tokens tok(buf);
    std::string t = tok.next();
    while (!t.empty() && std::isalpha(static_cast<unsigned char>(t[0])) != 0) {
        const std::string k = lower(t);
        const auto it = std::find(known.begin(), known.end(), k);
        if (it == known.end()) {
            fail(w + ": unknown header keyword " + t + "; an ESRI ASCII grid has ncols, nrows,"
                 " xllcorner or xllcenter, yllcorner or yllcenter, cellsize (or dx and dy) and an"
                 " optional NODATA_value, and its values are finite numbers");
        }
        const auto m = static_cast<std::size_t>(it - known.begin());
        if (seen[m] != 0) { fail(w + ": header keyword " + t + " is given twice"); }
        const std::string v = tok.next();
        if (!to_double(v, hv[m])) {
            fail(w + ": header keyword " + t + " needs a finite number; it has \"" + v + "\"");
        }
        seen[m] = 1;
        t = tok.next();
    }
    auto has = [&] (const char* k) {
        return seen[static_cast<std::size_t>(std::find(known.begin(), known.end(), k) - known.begin())] != 0;
    };
    auto get = [&] (const char* k) {
        return hv[static_cast<std::size_t>(std::find(known.begin(), known.end(), k) - known.begin())];
    };
    for (const char* k : {"ncols", "nrows"}) {
        if (!has(k)) { fail(w + ": the header has no " + std::string(k)); }
        const double v = get(k);
        if (!(v >= 1.0 && v <= 1073741824.0 && v == std::floor(v))) {
            fail(w + ": " + std::string(k) + " must be a whole number from 1 to 2^30; it is "
                 + all_digits(v));
        }
    }
    if (has("xllcorner") == has("xllcenter")) {
        fail(w + ": the header needs exactly one of xllcorner and xllcenter");
    }
    if (has("yllcorner") == has("yllcenter")) {
        fail(w + ": the header needs exactly one of yllcorner and yllcenter");
    }
    const bool cs = has("cellsize");
    if (cs == (has("dx") || has("dy")) || has("dx") != has("dy")) {
        fail(w + ": the header needs cellsize, or dx and dy, but not both");
    }
    const auto nx = static_cast<long long>(get("ncols"));
    const auto ny = static_cast<long long>(get("nrows"));
    const double dx = cs ? get("cellsize") : get("dx");
    const double dy = cs ? get("cellsize") : get("dy");
    if (!(dx > 0.0 && dy > 0.0)) {
        fail(w + ": the cell size (cellsize, or dx and dy) must be > 0; it is "
             + all_digits(dx) + " by " + all_digits(dy));
    }
    const double x0 = has("xllcorner") ? get("xllcorner") : get("xllcenter") - 0.5 * dx;
    const double y0 = has("yllcorner") ? get("yllcorner") : get("yllcenter") - 0.5 * dy;
    if (!(std::isfinite(x0 + static_cast<double>(nx) * dx)
          && std::isfinite(y0 + static_cast<double>(ny) * dy))) {
        fail(w + ": the raster's extent is not finite");
    }
    const bool has_nodata = has("nodata_value");
    const double nodata = has_nodata ? get("nodata_value") : 0.0;
    // Each value takes at least two characters (a digit and a separator):
    // a header that asks for more cannot be met, and nothing is allocated.
    const long long n = nx * ny;
    if (static_cast<double>(n) > 0.5 * static_cast<double>(buf.size()) + 1.0) {
        fail(w + ": the header says " + std::to_string(nx) + " x " + std::to_string(ny) + " = "
             + std::to_string(n) + " values, more than a file of " + std::to_string(buf.size())
             + " bytes can hold");
    }

    // The columns and rows that hold the domain, and one more each side.
    const double lo[2] = {static_cast<double>(geom.ProbLo(0)), static_cast<double>(geom.ProbLo(1))};
    const double hi[2] = {static_cast<double>(geom.ProbHi(0)), static_cast<double>(geom.ProbHi(1))};
    const double r0[2] = {x0, y0};
    const double rd[2] = {dx, dy};
    const long long rn[2] = {nx, ny};
    const char* const dir[2] = {"x", "y"};
    long long first[2] = {0, 0};
    long long last[2] = {0, 0};
    for (int d = 0; d < 2; ++d) {
        const double r1 = r0[d] + static_cast<double>(rn[d]) * rd[d];
        const double tol = 1.0e-6 * rd[d];
        if (!(lo[d] >= r0[d] - tol && hi[d] <= r1 + tol)) {
            fail(w + ": the raster covers " + all_digits(r0[d]) + " to " + all_digits(r1) + " in "
                 + dir[d] + ", which does not hold the domain, " + all_digits(lo[d]) + " to "
                 + all_digits(hi[d]) + " (flume.prob_lo, flume.prob_hi): check that the file is"
                 " on the run's coordinates and in metres");
        }
        const auto cell = [&] (double x) {
            return static_cast<long long>(std::floor((x - r0[d]) / rd[d]));
        };
        first[d] = std::clamp(cell(lo[d]) - 1, 0LL, rn[d] - 1);
        last[d] = std::clamp(cell(hi[d]) + 1, 0LL, rn[d] - 1);
    }
    const long long kx = last[0] - first[0] + 1;
    const long long ky = last[1] - first[1] + 1;
    if (kx * ky > static_cast<long long>(INT_MAX)) {
        fail(w + ": the part of the raster the domain reads has " + std::to_string(kx) + " x "
             + std::to_string(ky) + " cells, more than 2^31 - 1: use a coarser raster");
    }

    Raster r;
    r.key = key;
    r.path = path;
    r.nx = static_cast<int>(kx);
    r.ny = static_cast<int>(ky);
    r.dx = dx;
    r.dy = dy;
    r.x0 = x0 + static_cast<double>(first[0]) * dx;
    r.y0 = y0 + static_cast<double>(first[1]) * dy;
    r.value.assign(static_cast<std::size_t>(kx * ky), 0.0);
    r.nodata.assign(static_cast<std::size_t>(kx * ky), 0);

    // Values: the northern row first in the file, the southern first here.
    const double real_max = static_cast<double>(std::numeric_limits<amrex::Real>::max());
    long long count = 0;
    while (!t.empty()) {
        double v = 0.0;
        if (!to_double(t, v)) {
            fail(w + ": value " + std::to_string(count + 1) + " is not a finite number: \"" + t + "\"");
        }
        if (count < n) {
            const long long i = count % nx;
            const long long j = ny - 1 - count / nx;
            if (i >= first[0] && i <= last[0] && j >= first[1] && j <= last[1]) {
                const auto m = static_cast<std::size_t>((i - first[0]) + kx * (j - first[1]));
                if (has_nodata && v == nodata) {
                    r.nodata[m] = 1;
                } else if (!(std::abs(v) <= real_max)) {
                    fail(w + ": value " + std::to_string(count + 1) + " (column " + std::to_string(i)
                         + " from the west, row " + std::to_string(j) + " from the south) is "
                         + all_digits(v) + ", beyond the range of the build's precision");
                } else {
                    r.value[m] = v;
                }
            }
        }
        ++count;
        t = tok.next();
    }
    if (count != n) {
        fail(w + ": the file holds " + std::to_string(count) + " values; the header says "
             + std::to_string(nx) + " x " + std::to_string(ny) + " = " + std::to_string(n));
    }
    return r;
}

//! The raster cells as device arrays: values in amrex::Real, NODATA flags.
struct DeviceRaster
{
    amrex::Gpu::DeviceVector<amrex::Real> value;
    amrex::Gpu::DeviceVector<int> nodata;
};

DeviceRaster
to_device (const Raster& r)
{
    const std::size_t n = r.value.size();
    std::vector<amrex::Real> v(n);
    std::vector<int> flag(n);
    for (std::size_t m = 0; m < n; ++m) {
        v[m] = static_cast<amrex::Real>(r.value[m]);
        flag[m] = r.nodata[m];
    }
    DeviceRaster d;
    d.value.resize(n);
    d.nodata.resize(n);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, v.begin(), v.end(), d.value.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, flag.begin(), flag.end(), d.nodata.begin());
    amrex::Gpu::streamSynchronize();
    return d;
}

//! Abort, naming the key, if `n` cells (summed over the ranks) read NODATA.
void
abort_on_nodata (amrex::Long n, const Raster& r)
{
    amrex::ParallelDescriptor::ReduceLongSum(n);
    if (n > 0) {
        amrex::Abort(where(r.key, r.path) + ": " + std::to_string(n)
                     + " cell(s) of the domain read a NODATA_value cell of the raster");
    }
}

//! Where the cells of `geom` sit in the raster, in the precision of the build:
//! the offset of the domain's lower corner from the raster's west and south
//! edges, and the two cell sizes. The offset is the difference of the corner
//! (already in the precision of the build) and the raster's edge, formed in
//! double: it is at most a few raster cells, so it loses no digits.
struct Placement
{
    amrex::Real ox, oy;   //!< prob_lo - raster edge, m
    amrex::Real gdx, gdy; //!< cell size of the grid, m
    amrex::Real rdx, rdy; //!< cell size of the raster, m
    int ilo, jlo;         //!< first cell of the domain
    int nx, ny;           //!< raster columns and rows
};

Placement
placement (const Raster& r, const amrex::Geometry& geom)
{
    Placement p{};
    p.ox = static_cast<amrex::Real>(static_cast<double>(geom.ProbLo(0)) - r.x0);
    p.oy = static_cast<amrex::Real>(static_cast<double>(geom.ProbLo(1)) - r.y0);
    p.gdx = geom.CellSize(0);
    p.gdy = geom.CellSize(1);
    p.rdx = static_cast<amrex::Real>(r.dx);
    p.rdy = static_cast<amrex::Real>(r.dy);
    p.ilo = geom.Domain().smallEnd(0);
    p.jlo = geom.Domain().smallEnd(1);
    p.nx = r.nx;
    p.ny = r.ny;
    return p;
}

//! q, or the whole number nearest to it when within 16 epsilons of it (in
//! raster cells, relative to max(|q|, 1)).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
snap_whole (amrex::Real q) noexcept
{
    const amrex::Real whole = std::round(q);
    const amrex::Real tol = amrex::Real(16.0) * std::numeric_limits<amrex::Real>::epsilon()
                            * amrex::max(std::abs(q), amrex::Real(1.0));
    return (std::abs(q - whole) <= tol) ? whole : q;
}

//! The raster cell that holds position q (in raster cells from the edge):
//! a q within 16 epsilons of a whole number counts as that number, so a
//! centre on a raster edge takes the cell above it whatever the round-off.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
holding_cell (amrex::Real q, int n) noexcept
{
    return amrex::Clamp(static_cast<int>(std::floor(snap_whole(q))), 0, n - 1);
}

//! The length, in raster cells, that cell c (of n) shares with [a, b]
//! (a <= b, in raster cells from the raster's edge). The first and the last
//! cell extend without end beyond the raster's edges.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
shared_length (amrex::Real a, amrex::Real b, int c, int n) noexcept
{
    const amrex::Real lo = (c == 0) ? a : amrex::max(a, static_cast<amrex::Real>(c));
    const amrex::Real hi = (c == n - 1) ? b : amrex::min(b, static_cast<amrex::Real>(c + 1));
    return amrex::max(hi - lo, amrex::Real(0.0));
}

} // namespace

Raster
read_raster (const std::string& key, const std::string& path, const amrex::Geometry& geom)
{
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    Raster r;
    std::string error;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        try {
            r = parse(key, path, geom);
        } catch (const RasterError& e) {
            error = e.message;
        }
    }
    // The message (empty when the file is good), then the part kept.
    long long len = static_cast<long long>(error.size());
    amrex::ParallelDescriptor::Bcast(&len, 1, io);
    if (len > 0) {
        error.resize(static_cast<std::size_t>(len));
        amrex::ParallelDescriptor::Bcast(error.data(), static_cast<std::size_t>(len), io);
        amrex::Abort(error);
    }
    r.key = key;
    r.path = path;
    int ni[2] = {r.nx, r.ny};
    double geo[4] = {r.x0, r.y0, r.dx, r.dy};
    amrex::ParallelDescriptor::Bcast(ni, 2, io);
    amrex::ParallelDescriptor::Bcast(geo, 4, io);
    r.nx = ni[0];
    r.ny = ni[1];
    r.x0 = geo[0];
    r.y0 = geo[1];
    r.dx = geo[2];
    r.dy = geo[3];
    const auto n = static_cast<std::size_t>(r.nx) * static_cast<std::size_t>(r.ny);
    r.value.resize(n);
    r.nodata.resize(n);
    amrex::ParallelDescriptor::Bcast(r.value.data(), n, io);
    amrex::ParallelDescriptor::Bcast(r.nodata.data(), n, io);
    return r;
}

void
sample_bilinear (amrex::MultiFab& dst, int comp, const Raster& r, const amrex::Geometry& geom)
{
    const Placement p = placement(r, geom);
    const DeviceRaster dr = to_device(r);
    const amrex::Real* v = dr.value.data();
    const int* flag = dr.nodata.data();
    const auto& a = dst.arrays();
    const amrex::Long bad = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                             amrex::TypeList<amrex::Long>{}, dst, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            const amrex::Real zero = amrex::Real(0.0);
            const amrex::Real one = amrex::Real(1.0);
            const amrex::Real half = amrex::Real(0.5);
            // Position in raster cells from the first cell centre, held at
            // the first and last centres.
            const amrex::Real xr = p.ox + (static_cast<amrex::Real>(i - p.ilo) + half) * p.gdx;
            const amrex::Real yr = p.oy + (static_cast<amrex::Real>(j - p.jlo) + half) * p.gdy;
            const amrex::Real fx = amrex::Clamp(xr / p.rdx - half, zero,
                                                static_cast<amrex::Real>(p.nx - 1));
            const amrex::Real fy = amrex::Clamp(yr / p.rdy - half, zero,
                                                static_cast<amrex::Real>(p.ny - 1));
            const int i0 = amrex::min(static_cast<int>(std::floor(fx)), amrex::max(p.nx - 2, 0));
            const int j0 = amrex::min(static_cast<int>(std::floor(fy)), amrex::max(p.ny - 2, 0));
            const int i1 = amrex::min(i0 + 1, p.nx - 1);
            const int j1 = amrex::min(j0 + 1, p.ny - 1);
            const amrex::Real wx = fx - static_cast<amrex::Real>(i0);
            const amrex::Real wy = fy - static_cast<amrex::Real>(j0);
            const amrex::Long nx = p.nx;
            const amrex::Long m00 = i0 + nx * j0;
            const amrex::Long m10 = i1 + nx * j0;
            const amrex::Long m01 = i0 + nx * j1;
            const amrex::Long m11 = i1 + nx * j1;
            const amrex::Real w00 = (one - wx) * (one - wy);
            const amrex::Real w10 = wx * (one - wy);
            const amrex::Real w01 = (one - wx) * wy;
            const amrex::Real w11 = wx * wy;
            const bool no = (w00 > zero && flag[m00] != 0) || (w10 > zero && flag[m10] != 0)
                            || (w01 > zero && flag[m01] != 0) || (w11 > zero && flag[m11] != 0);
            a[box_no](i, j, k, comp) = w00 * v[m00] + w10 * v[m10] + w01 * v[m01] + w11 * v[m11];
            return {no ? amrex::Long(1) : amrex::Long(0)};
        });
    abort_on_nodata(bad, r);
}

void
sample_nearest (amrex::MultiFab& dst, int comp, const Raster& r, const amrex::Geometry& geom)
{
    const Placement p = placement(r, geom);
    const DeviceRaster dr = to_device(r);
    const amrex::Real* v = dr.value.data();
    const int* flag = dr.nodata.data();
    const auto& a = dst.arrays();
    const amrex::Long bad = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                             amrex::TypeList<amrex::Long>{}, dst, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            const amrex::Real half = amrex::Real(0.5);
            const amrex::Real xr = p.ox + (static_cast<amrex::Real>(i - p.ilo) + half) * p.gdx;
            const amrex::Real yr = p.oy + (static_cast<amrex::Real>(j - p.jlo) + half) * p.gdy;
            const int ci = holding_cell(xr / p.rdx, p.nx);
            const int cj = holding_cell(yr / p.rdy, p.ny);
            const amrex::Long m = ci + amrex::Long(p.nx) * cj;
            a[box_no](i, j, k, comp) = v[m];
            return {(flag[m] != 0) ? amrex::Long(1) : amrex::Long(0)};
        });
    abort_on_nodata(bad, r);
}

void
sample_conservative (amrex::MultiFab& dst, int comp, const Raster& r, const amrex::Geometry& geom)
{
    const Placement p = placement(r, geom);
    const DeviceRaster dr = to_device(r);
    const amrex::Real* v = dr.value.data();
    const int* flag = dr.nodata.data();
    const auto& a = dst.arrays();
    const amrex::Long bad = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                             amrex::TypeList<amrex::Long>{}, dst, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            const amrex::Real zero = amrex::Real(0.0);
            // The cell's edges in raster cells from the raster's edge: the
            // west edge of cell i is the east edge of cell i - 1, bit for bit.
            // An edge within round-off of a raster edge is put on it, so a
            // raster cell beyond it is not read through a sliver of round-off.
            const auto edge_x = [&] (int ii) {
                return snap_whole((p.ox + static_cast<amrex::Real>(ii - p.ilo) * p.gdx) / p.rdx);
            };
            const auto edge_y = [&] (int jj) {
                return snap_whole((p.oy + static_cast<amrex::Real>(jj - p.jlo) * p.gdy) / p.rdy);
            };
            const amrex::Real xa = edge_x(i);
            const amrex::Real xb = edge_x(i + 1);
            const amrex::Real ya = edge_y(j);
            const amrex::Real yb = edge_y(j + 1);
            const int c0 = amrex::Clamp(static_cast<int>(std::floor(xa)), 0, p.nx - 1);
            const int c1 = amrex::Clamp(static_cast<int>(std::floor(xb)), 0, p.nx - 1);
            const int r0 = amrex::Clamp(static_cast<int>(std::floor(ya)), 0, p.ny - 1);
            const int r1 = amrex::Clamp(static_cast<int>(std::floor(yb)), 0, p.ny - 1);
            amrex::Real sum = zero;
            amrex::Real area = zero;
            amrex::Long no = 0;
            for (int rr = r0; rr <= r1; ++rr) {
                const amrex::Real ly = shared_length(ya, yb, rr, p.ny);
                if (ly <= zero) { continue; }
                for (int cc = c0; cc <= c1; ++cc) {
                    const amrex::Real lx = shared_length(xa, xb, cc, p.nx);
                    if (lx <= zero) { continue; }
                    const amrex::Long m = cc + amrex::Long(p.nx) * rr;
                    if (flag[m] != 0) { no = 1; }
                    sum += lx * ly * v[m];
                    area += lx * ly;
                }
            }
            // area is the cell's area in raster cells (the shared lengths of
            // a row of raster cells sum to xb - xa); dividing by it rather
            // than by (xb - xa) (yb - ya) makes a uniform raster map to its
            // value to round-off (to the last bit when few raster cells feed
            // a cell).
            a[box_no](i, j, k, comp) = (area > zero) ? sum / area : zero;
            return {no};
        });
    abort_on_nodata(bad, r);
}

} // namespace flume
