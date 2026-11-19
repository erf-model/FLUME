#include "FLUME_Inputs.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
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

    // max_grid_size: optional, one integer in [1, INT_MAX]
    if (pp.contains("max_grid_size")) {
        require_count(pp, "max_grid_size", 1);
        long long v = 0;
        pp.get("max_grid_size", v);
        in.max_grid_size = checked_int(v, 1, std::numeric_limits<int>::max(), "max_grid_size", "");
    }

    // plot_file: optional, one string, not empty
    if (pp.contains("plot_file")) {
        require_count(pp, "plot_file", 1);
        pp.get("plot_file", in.plot_file);
        if (in.plot_file.empty()) {
            amrex::Abort("flume.plot_file must not be empty");
        }
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
