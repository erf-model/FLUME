#include "FLUME_Inputs.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace flume {

namespace {

//! Abort unless the key holds exactly `want` values.
void
require_count (const amrex::ParmParse& pp, const std::string& key, int want)
{
    const int have = pp.countval(key);
    if (have != want) {
        std::string msg = "flume.";
        msg += key;
        msg += (want == 2) ? " needs exactly 2 values (x and y); it has "
                           : " needs exactly 1 value; it has ";
        msg += std::to_string(have);
        amrex::Abort(msg);
    }
}

//! Abort unless lo <= v <= hi. The value is read as a long integer, so a
//! value too large for an int is reported as typed and does not wrap.
int
checked_int (long long v, long long lo, long long hi, const std::string& key,
             const std::string& where)
{
    if (v < lo || v > hi) {
        std::string msg = "flume.";
        msg += key;
        msg += " must be >= ";
        msg += std::to_string(lo);
        msg += " and <= ";
        msg += std::to_string(hi);
        msg += where;
        msg += "; it is ";
        msg += std::to_string(v);
        amrex::Abort(msg);
    }
    return static_cast<int>(v);
}

//! A real as text with every digit, so that a value just outside a limit
//! does not print as the limit.
template <typename T>
std::string
all_digits (T v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<T>::max_digits10) << v;
    return s.str();
}

//! Read one finite real (float or double) if the key is given. Returns
//! whether it was.
template <typename T>
bool
query_finite (const amrex::ParmParse& pp, const std::string& key, T& value)
{
    if (!pp.contains(key)) { return false; }
    require_count(pp, key, 1);
    T v = T(0.0);
    pp.get(key, v);
    if (!std::isfinite(v)) {
        std::string msg = "flume.";
        msg += key;
        msg += " must be finite";
        amrex::Abort(msg);
    }
    value = v;
    return true;
}

//! Read two finite reals (x and y) if the key is given. Returns whether it was.
bool
query_finite_pair (const amrex::ParmParse& pp, const std::string& key,
                   amrex::Array<amrex::Real, 2>& value)
{
    if (!pp.contains(key)) { return false; }
    require_count(pp, key, 2);
    amrex::Vector<amrex::Real> v;
    pp.getarr(key, v);
    for (int d = 0; d < 2; ++d) {
        if (!std::isfinite(v[d])) {
            std::string msg = "flume.";
            msg += key;
            msg += " must be finite in x and y";
            amrex::Abort(msg);
        }
        value[d] = v[d];
    }
    return true;
}

//! Read one integer in [lo, hi] if the key is given. Returns whether it was.
bool
query_int (const amrex::ParmParse& pp, const std::string& key, long long lo, long long hi,
           int& value)
{
    if (!pp.contains(key)) { return false; }
    require_count(pp, key, 1);
    long long v = 0;
    pp.get(key, v);
    value = checked_int(v, lo, hi, key, "");
    return true;
}

} // namespace

Inputs
read_inputs ()
{
    Inputs in;
    const amrex::ParmParse pp("flume");
    const char* const dir_name[2] = {"x", "y"};

    // n_cell: required, two integers, each in [1, n_cell_max]
    {
        if (!pp.contains("n_cell")) {
            amrex::Abort("flume.n_cell is required: the number of cells in x and y");
        }
        require_count(pp, "n_cell", 2);
        std::vector<long long> v;
        pp.getarr("n_cell", v);
        for (int d = 0; d < 2; ++d) {
            in.n_cell[d] = checked_int(v[d], 1, n_cell_max, "n_cell", " in x and y");
        }
    }

    // prob_lo: optional, two finite reals
    if (pp.contains("prob_lo")) {
        require_count(pp, "prob_lo", 2);
        amrex::Vector<amrex::Real> v;
        pp.getarr("prob_lo", v);
        for (int d = 0; d < 2; ++d) {
            if (!std::isfinite(v[d])) {
                amrex::Abort("flume.prob_lo must be finite in x and y");
            }
            in.prob_lo[d] = v[d];
        }
    }

    // prob_hi: required, two finite reals, each above prob_lo
    {
        if (!pp.contains("prob_hi")) {
            amrex::Abort("flume.prob_hi is required: the upper corner of the domain in x and y, m");
        }
        require_count(pp, "prob_hi", 2);
        amrex::Vector<amrex::Real> v;
        pp.getarr("prob_hi", v);
        for (int d = 0; d < 2; ++d) {
            if (!std::isfinite(v[d])) {
                amrex::Abort("flume.prob_hi must be finite in x and y");
            }
            if (!(v[d] > in.prob_lo[d])) {
                amrex::Abort("flume.prob_hi must be greater than flume.prob_lo in x and y");
            }
            in.prob_hi[d] = v[d];
        }
    }

    // The cell size (prob_hi - prob_lo) / n_cell must be one the precision of
    // the build can tell apart from the coordinates: finite, a normal number,
    // and at least machine epsilon times the largest coordinate. Below that,
    // neighbouring cell edges fall on the same floating-point number.
    for (int d = 0; d < 2; ++d) {
        const amrex::Real extent = in.prob_hi[d] - in.prob_lo[d];
        const amrex::Real dx = extent / static_cast<amrex::Real>(in.n_cell[d]);
        const amrex::Real big = std::max(std::abs(in.prob_lo[d]), std::abs(in.prob_hi[d]));
        const amrex::Real dx_min = std::max(std::numeric_limits<amrex::Real>::min(),
                                            std::numeric_limits<amrex::Real>::epsilon() * big);
        // Room for prob_hi - prob_lo and for AMReX's own sums of coordinates.
        const amrex::Real big_max = std::numeric_limits<amrex::Real>::max() / amrex::Real(8.0);
        if (big > big_max) {
            std::ostringstream msg;
            msg << "flume.prob_lo and flume.prob_hi must be at most " << big_max
                << " in magnitude in x and y";
            amrex::Abort(msg.str());
        }
        if (!(dx >= dx_min)) {
            std::ostringstream msg;
            msg << "flume.prob_lo, flume.prob_hi and flume.n_cell give a cell size in "
                << dir_name[d] << " that this build's floating-point precision cannot resolve: "
                << dx << " m; the smallest it can resolve at these coordinates is " << dx_min
                << " m";
            amrex::Abort(msg.str());
        }
    }

    const long long int_max = std::numeric_limits<int>::max();

    // max_grid_size: optional, one integer in [1, INT_MAX]
    query_int(pp, "max_grid_size", 1, int_max, in.max_grid_size);

    // plot_file: optional, one string, not empty
    if (pp.contains("plot_file")) {
        require_count(pp, "plot_file", 1);
        pp.get("plot_file", in.plot_file);
        if (in.plot_file.empty()) {
            amrex::Abort("flume.plot_file must not be empty");
        }
    }

    // plot_int: optional, one integer in [0, INT_MAX]
    query_int(pp, "plot_int", 0, int_max, in.plot_int);

    // stop_time: optional, one finite real >= 0
    {
        double v = 0.0;
        if (query_finite(pp, "stop_time", v)) {
            if (v < 0.0) {
                amrex::Abort("flume.stop_time must be >= 0; it is " + all_digits(v));
            }
            in.stop_time = v;
        }
    }

    // max_step: optional, one integer in [0, INT_MAX]
    {
        int v = 0;
        if (query_int(pp, "max_step", 0, int_max, v)) { in.max_step = v; }
    }

    // cfl: optional, one real in (0, 1]
    query_finite(pp, "cfl", in.cfl);
    if (!(in.cfl > amrex::Real(0.0) && in.cfl <= amrex::Real(1.0))) {
        amrex::Abort("flume.cfl must be > 0 and <= 1; it is " + all_digits(in.cfl));
    }

    // gravity: optional, one finite real > 0
    query_finite(pp, "gravity", in.gravity);
    if (!(in.gravity > amrex::Real(0.0))) {
        amrex::Abort("flume.gravity must be > 0; it is " + all_digits(in.gravity));
    }

    // bed.elevation, bed.slope: optional, finite
    query_finite(pp, "bed.elevation", in.bed_elevation);
    query_finite_pair(pp, "bed.slope", in.bed_slope);

    // init.surface_elevation: optional, finite. Without it the domain is dry,
    // and a velocity would have nothing to move.
    {
        amrex::Real v = amrex::Real(0.0);
        if (query_finite(pp, "init.surface_elevation", v)) { in.init_surface_elevation = v; }
    }
    if (query_finite_pair(pp, "init.velocity", in.init_velocity)
        && !in.init_surface_elevation.has_value()) {
        amrex::Abort("flume.init.velocity needs flume.init.surface_elevation: without a water"
                     " surface the domain starts dry");
    }

    // A key in the flume namespace that nothing above read is an error.
    if (amrex::ParmParse::hasUnusedInputs("flume")) {
        std::string msg = "Unknown key(s) in the flume namespace:";
        for (const auto& s : amrex::ParmParse::getUnusedInputs("flume")) {
            msg += "\n  ";
            msg += s;
        }
        amrex::Abort(msg);
    }

    return in;
}

} // namespace flume
