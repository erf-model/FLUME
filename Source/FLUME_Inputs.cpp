#include "FLUME_Inputs.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
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

//! Read one file name, not empty, if the key is given. Returns whether it
//! was. `what` names the kind of file in the message.
bool
query_file (const amrex::ParmParse& pp, const std::string& key, std::string& value,
            const char* what = "an ESRI ASCII grid file")
{
    if (!pp.contains(key)) { return false; }
    require_count(pp, key, 1);
    pp.get(key, value);
    if (value.empty()) {
        amrex::Abort("flume." + key + " must not be empty: the name of " + what);
    }
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

    // periodic: optional, two integers, each 0 or 1
    if (pp.contains("periodic")) {
        require_count(pp, "periodic", 2);
        std::vector<long long> v;
        pp.getarr("periodic", v);
        for (int d = 0; d < 2; ++d) {
            in.periodic[d] = checked_int(v[d], 0, 1, "periodic", " in x and y");
        }
    }

    // A key that needs another one is reported after the check for unknown
    // keys below: the other key may be there, misspelled. `needs` holds the
    // first such message.
    std::string needs;
    const auto defer = [&needs] (const std::string& msg) { if (needs.empty()) { needs = msg; } };

    // bc.segments and bc.<name>.side, .type, .range: segments of the sides
    // that are not periodic, each with its own type. Checked here against
    // the grid: on a side that is not periodic, inside it, holding at least
    // one cell centre, not overlapping another segment.
    if (pp.contains("bc.segments")) {
        // AMReX's parser stops a key with no value, so there is at least one name.
        std::vector<std::string> names;
        pp.getarr("bc.segments", names);
        for (const auto& name : names) {
            const bool ok_name = !name.empty()
                && std::all_of(name.begin(), name.end(), [] (char c) {
                       return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_';
                   });
            if (!ok_name) {
                amrex::Abort("flume.bc.segments: a name may hold only letters, digits and _; "
                             "it is \"" + name + "\"");
            }
            for (const auto& s : in.bc_segments) {
                if (s.name == name) { amrex::Abort("flume.bc.segments lists " + name + " twice"); }
            }
            BoundarySegment seg;
            seg.name = name;
            const std::string key = "bc." + name + ".";

            // A segment without its side or type is reported after the check
            // for unknown keys (the key may be there, misspelled); its other
            // keys are marked as read so that they are not reported too.
            const bool has_side = pp.contains(key + "side");
            const bool has_type = pp.contains(key + "type");
            if (!has_side || !has_type) {
                defer(!has_side ? "flume." + key + "side is required: xlo, xhi, ylo or yhi"
                                : "flume." + key + "type is required: wall or outflow");
                static_cast<void>(pp.contains(key + "range"));
                continue;
            }
            require_count(pp, key + "side", 1);
            std::string side;
            pp.get(key + "side", side);
            if (side == "xlo") {
                seg.direction = 0; seg.high = 0;
            } else if (side == "xhi") {
                seg.direction = 0; seg.high = 1;
            } else if (side == "ylo") {
                seg.direction = 1; seg.high = 0;
            } else if (side == "yhi") {
                seg.direction = 1; seg.high = 1;
            } else {
                amrex::Abort("flume." + key + "side must be xlo, xhi, ylo or yhi; it is " + side);
            }

            require_count(pp, key + "type", 1);
            std::string type;
            pp.get(key + "type", type);
            if (type == "wall") {
                seg.type = BoundaryType::wall;
            } else if (type == "outflow") {
                seg.type = BoundaryType::outflow;
            } else if (type == "discharge" || type == "stage") {
                amrex::Abort("flume." + key + "type = " + type + " is reserved for a later version"
                             " (a prescribed " + (type == "discharge" ? "inflow discharge" : "water level")
                             + " from a time series) and is not available yet: use wall or outflow");
            } else {
                amrex::Abort("flume." + key + "type must be wall or outflow; it is " + type);
            }

            const int d = seg.direction;
            const int t = 1 - d;
            if (in.periodic[d] == 1) {
                amrex::Abort("flume." + key + "side is " + side + ", but flume.periodic is 1 in "
                             + dir_name[d] + ": a periodic side has no boundary to set");
            }
            const amrex::Real lo = in.prob_lo[t];
            const amrex::Real hi = in.prob_hi[t];
            amrex::Real r_lo = lo;
            amrex::Real r_hi = hi;
            if (pp.contains(key + "range")) {
                const int have = pp.countval(key + "range");
                if (have != 2) {
                    amrex::Abort("flume." + key + "range needs exactly 2 values (from and to, in "
                                 + dir_name[t] + "); it has " + std::to_string(have));
                }
                amrex::Vector<amrex::Real> r;
                pp.getarr(key + "range", r);
                r_lo = r[0];
                r_hi = r[1];
                if (!std::isfinite(r_lo) || !std::isfinite(r_hi)) {
                    std::string msg = "flume." + key + "range must be finite; it is ";
                    msg += all_digits(r_lo);
                    msg += " ";
                    msg += all_digits(r_hi);
                    amrex::Abort(msg);
                }
                if (!(r_lo < r_hi && r_lo >= lo && r_hi <= hi)) {
                    std::string msg = "flume." + key + "range must be two coordinates along the side,"
                                      " in ";
                    msg += dir_name[t];
                    msg += ", the first below the second, both from ";
                    msg += all_digits(lo);
                    msg += " to ";
                    msg += all_digits(hi);
                    msg += " (flume.prob_lo and flume.prob_hi); it is ";
                    msg += all_digits(r_lo);
                    msg += " ";
                    msg += all_digits(r_hi);
                    amrex::Abort(msg);
                }
            }
            // The cells whose centre is at r_lo or above and below r_hi, with
            // the centres computed as the plotfile and the kernels place them.
            const int n = in.n_cell[t];
            const amrex::Real dx = (hi - lo) / static_cast<amrex::Real>(n);
            auto centre = [&] (long long j) {
                return lo + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx;
            };
            // The first cell whose centre is at x or above (n if none): a
            // guess from the division, moved to the exact answer.
            auto first_at_or_above = [&] (amrex::Real x) {
                long long j = static_cast<long long>(std::ceil((x - lo) / dx - amrex::Real(0.5)));
                j = std::clamp(j, 0LL, static_cast<long long>(n));
                while (j > 0 && centre(j - 1) >= x) { --j; }
                while (j < n && centre(j) < x) { ++j; }
                return static_cast<int>(j);
            };
            // An end at the side's own end is its end cell, exactly: in single
            // precision the last centre of a long side can round up to hi.
            const int first = (r_lo == lo) ? 0 : first_at_or_above(r_lo);
            const int last = (r_hi == hi) ? n - 1 : first_at_or_above(r_hi) - 1;
            if (last < first) {
                std::string msg = "flume." + key + "range holds no cell centre: a segment takes the"
                                  " cells whose centre is at the first coordinate or above and below"
                                  " the second; it is ";
                msg += all_digits(r_lo);
                msg += " ";
                msg += all_digits(r_hi);
                amrex::Abort(msg);
            }
            seg.cell_lo = first;
            seg.cell_hi = last;
            for (const auto& s : in.bc_segments) {
                if (s.direction == d && s.high == seg.high && s.cell_lo <= last
                    && first <= s.cell_hi) {
                    amrex::Abort("flume.bc." + s.name + " and flume.bc." + name
                                 + " overlap on side " + side + ": a boundary cell can belong to"
                                 " one segment only");
                }
            }
            in.bc_segments.push_back(seg);
        }
    }

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

    // order: optional, one integer, 1 or 2
    query_int(pp, "order", 1, 2, in.order);
    // At second order a side of a cell holds up to 1.5 times its depth, and
    // the scheme loses its order above half the first-order bound.
    if (in.order == 2 && !(in.cfl <= amrex::Real(0.5))) {
        amrex::Abort("flume.cfl must be <= 0.5 with flume.order = 2; it is " + all_digits(in.cfl));
    }

    // gravity: optional, one finite real > 0
    query_finite(pp, "gravity", in.gravity);
    if (!(in.gravity > amrex::Real(0.0))) {
        amrex::Abort("flume.gravity must be > 0; it is " + all_digits(in.gravity));
    }

    // dry_depth: optional, one finite real >= 0
    query_finite(pp, "dry_depth", in.dry_depth);
    if (!(in.dry_depth >= amrex::Real(0.0))) {
        amrex::Abort("flume.dry_depth must be >= 0; it is " + all_digits(in.dry_depth));
    }

    // friction.manning: optional, one finite real >= 0
    query_finite(pp, "friction.manning", in.friction_manning);
    if (!(in.friction_manning >= amrex::Real(0.0))) {
        amrex::Abort("flume.friction.manning must be >= 0; it is "
                     + all_digits(in.friction_manning));
    }

    // friction.manning_file: optional, one file name, not empty; not with
    // friction.manning, which would be read and not used
    const bool manning_value = pp.contains("friction.manning");
    if (query_file(pp, "friction.manning_file", in.friction_manning_file) && manning_value) {
        defer("flume.friction.manning_file cannot be used with flume.friction.manning: the file"
              " gives the coefficient of every cell");
    }

    // rain.type: none, constant or hyetograph; rain.rate (constant only),
    // rain.file and rain.interpolation (hyetograph only)
    if (pp.contains("rain.type")) {
        require_count(pp, "rain.type", 1);
        std::string v;
        pp.get("rain.type", v);
        if (v == "none") {
            in.rain_type = RainType::none;
        } else if (v == "constant") {
            in.rain_type = RainType::constant;
        } else if (v == "hyetograph") {
            in.rain_type = RainType::hyetograph;
        } else if (v == "gridded" || v == "moving_cell" || v == "microphysics") {
            amrex::Abort("flume.rain.type = " + v + " is reserved for a later version and is not"
                         " available yet: use none, constant or hyetograph");
        } else {
            amrex::Abort("flume.rain.type must be none, constant or hyetograph; it is " + v);
        }
    }
    {
        const bool rate = query_finite(pp, "rain.rate", in.rain_rate);
        if (!(in.rain_rate >= 0.0)) {
            amrex::Abort("flume.rain.rate must be >= 0 (mm/h); it is " + all_digits(in.rain_rate));
        }
        if (in.rain_type == RainType::constant && !rate) {
            defer("flume.rain.type = constant needs flume.rain.rate: the rain rate in mm/h");
        }
        if (in.rain_type != RainType::constant && rate) {
            defer("flume.rain.rate is read only with flume.rain.type = constant");
        }
        const bool file = query_file(pp, "rain.file", in.rain_file,
                                     "a hyetograph file (a time in s and a rain rate in mm/h"
                                     " per line)");
        if (in.rain_type == RainType::hyetograph && !file) {
            defer("flume.rain.type = hyetograph needs flume.rain.file: lines of a time in s and"
                  " a rain rate in mm/h");
        }
        if (in.rain_type != RainType::hyetograph && file) {
            defer("flume.rain.file is read only with flume.rain.type = hyetograph");
        }
        if (pp.contains("rain.interpolation")) {
            require_count(pp, "rain.interpolation", 1);
            std::string v;
            pp.get("rain.interpolation", v);
            if (v == "step") {
                in.rain_interpolation = RainInterpolation::step;
            } else if (v == "linear") {
                in.rain_interpolation = RainInterpolation::linear;
            } else {
                amrex::Abort("flume.rain.interpolation must be step or linear; it is " + v);
            }
            if (in.rain_type != RainType::hyetograph) {
                defer("flume.rain.interpolation is read only with flume.rain.type = hyetograph");
            }
        }
    }

    // bed.elevation, bed.slope: optional, finite
    const bool bed_elevation = query_finite(pp, "bed.elevation", in.bed_elevation);
    const bool bed_slope = query_finite_pair(pp, "bed.slope", in.bed_slope);

    // bed.continue_slope: optional, one integer, 0 or 1
    const bool continue_slope = query_int(pp, "bed.continue_slope", 0, 1,
                                          in.bed_continue_slope);

    // bed.noise: optional, one finite real >= 0
    const bool bed_noise = query_finite(pp, "bed.noise", in.bed_noise);
    if (!(in.bed_noise >= amrex::Real(0.0))) {
        amrex::Abort("flume.bed.noise must be >= 0; it is " + all_digits(in.bed_noise));
    }

    // init.surface_elevation: optional, finite. Without it the domain is dry,
    // and a velocity would have nothing to move.
    {
        amrex::Real v = amrex::Real(0.0);
        if (query_finite(pp, "init.surface_elevation", v)) { in.init_surface_elevation = v; }
    }
    // bed.continue_slope = 1 acts on a periodic side across which the bed
    // slopes; anywhere else it would be read and not used.
    if (continue_slope && in.bed_continue_slope == 1
        && !((in.periodic[0] == 1 && in.bed_slope[0] != amrex::Real(0.0))
             || (in.periodic[1] == 1 && in.bed_slope[1] != amrex::Real(0.0)))) {
        defer("flume.bed.continue_slope = 1 needs flume.periodic = 1 in a direction where"
              " flume.bed.slope is not 0");
    }

    // bed.bowl_depth: optional, finite; bed.bowl_radius: finite, > 0, and
    // each needs the other (a depth of 0 is no bowl, and reads no radius).
    const bool bowl_depth = query_finite(pp, "bed.bowl_depth", in.bed_bowl_depth);
    const bool bowl_radius = query_finite(pp, "bed.bowl_radius", in.bed_bowl_radius);
    if (bowl_radius && !(in.bed_bowl_radius > amrex::Real(0.0))) {
        amrex::Abort("flume.bed.bowl_radius must be > 0; it is " + all_digits(in.bed_bowl_radius));
    }
    // The bed uses 1 / radius^2.
    if (bowl_radius && !(in.bed_bowl_radius * in.bed_bowl_radius
                         >= std::numeric_limits<amrex::Real>::min())) {
        amrex::Abort("flume.bed.bowl_radius must have a square that is a normal number; it is "
                     + all_digits(in.bed_bowl_radius));
    }
    const bool bowl = bowl_depth && in.bed_bowl_depth != amrex::Real(0.0);
    if (bowl && !bowl_radius) {
        defer("flume.bed.bowl_depth needs flume.bed.bowl_radius: the distance from the centre"
              " of the domain at which the bowl crosses the plane bed");
    }
    if (bowl_radius && !bowl) {
        defer("flume.bed.bowl_radius needs flume.bed.bowl_depth (not 0)");
    }

    // bed.file: optional, one file name, not empty. It is the whole bed, so
    // the keys of the analytic bed would be read and not used.
    if (query_file(pp, "bed.file", in.bed_file)) {
        const std::pair<bool, const char*> others[] = {
            {bed_elevation, "bed.elevation"}, {bed_slope, "bed.slope"},
            {continue_slope, "bed.continue_slope"}, {bowl_depth, "bed.bowl_depth"},
            {bowl_radius, "bed.bowl_radius"}, {bed_noise, "bed.noise"}};
        for (const auto& [given, name] : others) {
            if (given) {
                defer(std::string("flume.bed.file cannot be used with flume.") + name
                      + ": the file gives the bed of every cell");
            }
        }
    }

    // buildings.file: optional, one file name, not empty; buildings.height:
    // optional, finite, > 0, needs buildings.file
    const bool buildings = query_file(pp, "buildings.file", in.buildings_file);
    {
        amrex::Real v = amrex::Real(0.0);
        if (query_finite(pp, "buildings.height", v)) {
            if (!(v > amrex::Real(0.0))) {
                amrex::Abort("flume.buildings.height must be > 0; it is " + all_digits(v));
            }
            in.buildings_height = v;
            if (!buildings) {
                defer("flume.buildings.height needs flume.buildings.file: the footprints the"
                      " height is given to");
            }
        }
    }

    // init.surface_slope, init.surface_curvature: optional, finite
    const bool surface_slope = query_finite_pair(pp, "init.surface_slope", in.init_surface_slope);
    const bool surface_curvature = query_finite(pp, "init.surface_curvature",
                                                in.init_surface_curvature);
    if (surface_slope && !in.init_surface_elevation.has_value()) {
        defer("flume.init.surface_slope needs flume.init.surface_elevation: the elevation of"
              " the surface at the centre of the domain");
    }
    if (surface_curvature && !in.init_surface_elevation.has_value()) {
        defer("flume.init.surface_curvature needs flume.init.surface_elevation: the elevation"
              " of the surface at the centre of the domain");
    }
    // init.wave_amplitude: optional, finite; init.wave_direction: x or y
    const bool wave_amplitude = query_finite(pp, "init.wave_amplitude", in.init_wave_amplitude);
    const bool wave = wave_amplitude && in.init_wave_amplitude != amrex::Real(0.0);
    if (wave_amplitude && !in.init_surface_elevation.has_value()) {
        defer("flume.init.wave_amplitude needs flume.init.surface_elevation: the elevation"
              " of the level surface the wave is added to");
    }
    if (pp.contains("init.wave_direction")) {
        require_count(pp, "init.wave_direction", 1);
        std::string v;
        pp.get("init.wave_direction", v);
        if (v == "x") {
            in.init_wave_direction = 0;
        } else if (v == "y") {
            in.init_wave_direction = 1;
        } else {
            amrex::Abort("flume.init.wave_direction must be x or y; it is " + v);
        }
        if (!wave) { defer("flume.init.wave_direction needs flume.init.wave_amplitude (not 0)"); }
    }
    if (wave && (surface_slope || surface_curvature)) {
        defer("flume.init.wave_amplitude cannot be used with flume.init.surface_slope or"
              " flume.init.surface_curvature: the wave is added to a level surface");
    }
    if (wave && in.periodic[in.init_wave_direction] != 1) {
        defer("flume.init.wave_amplitude needs flume.periodic = 1 along"
              " flume.init.wave_direction: one wavelength spans the domain");
    }
    if (query_finite_pair(pp, "init.velocity", in.init_velocity)
        && !in.init_surface_elevation.has_value()) {
        defer("flume.init.velocity needs flume.init.surface_elevation: without a water"
              " surface the domain starts dry");
    }

    // init.dam_position, init.dam_surface_elevation, init.dam_direction: a dam
    // needs the surface on both of its sides, and the other two keys need a dam.
    {
        amrex::Real v = amrex::Real(0.0);
        if (query_finite(pp, "init.dam_position", v)) { in.init_dam_position = v; }
    }
    const bool dam = in.init_dam_position.has_value();
    const bool dam_surface = query_finite(pp, "init.dam_surface_elevation",
                                          in.init_dam_surface_elevation);
    if (pp.contains("init.dam_direction")) {
        require_count(pp, "init.dam_direction", 1);
        std::string v;
        pp.get("init.dam_direction", v);
        if (v == "x") {
            in.init_dam_direction = 0;
        } else if (v == "y") {
            in.init_dam_direction = 1;
        } else {
            amrex::Abort("flume.init.dam_direction must be x or y; it is " + v);
        }
        if (!dam) { defer("flume.init.dam_direction needs flume.init.dam_position"); }
    }
    // The velocity beyond the dam is that of the near side unless given.
    in.init_dam_velocity = in.init_velocity;
    if (query_finite_pair(pp, "init.dam_velocity", in.init_dam_velocity) && !dam) {
        defer("flume.init.dam_velocity needs flume.init.dam_position");
    }
    if (dam_surface && !dam) {
        defer("flume.init.dam_surface_elevation needs flume.init.dam_position");
    }
    if (dam && !in.init_surface_elevation.has_value()) {
        defer("flume.init.dam_position needs flume.init.surface_elevation: the water"
              " surface on the low-coordinate side of the dam");
    }
    if (dam && !dam_surface) {
        defer("flume.init.dam_position needs flume.init.dam_surface_elevation: the"
              " water surface on the high-coordinate side of the dam");
    }
    if (dam && surface_slope) {
        defer("flume.init.surface_slope cannot be used with flume.init.dam_position: the"
              " surface on each side of a dam is level");
    }
    if (dam && wave) {
        defer("flume.init.wave_amplitude cannot be used with flume.init.dam_position: the"
              " surface on each side of a dam is level");
    }
    if (dam && surface_curvature) {
        defer("flume.init.surface_curvature cannot be used with flume.init.dam_position: the"
              " surface on each side of a dam is level");
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
    if (!needs.empty()) { amrex::Abort(needs); }

    if (dam) {
        // A cell is on the low side when its centre is below the dam (as
        // fill_state_dam computes it). Each side must hold at least one cell:
        // otherwise the keys of the other side would be read and not used.
        const int d = in.init_dam_direction;
        const amrex::Real x = *in.init_dam_position;
        const amrex::Real dx = (in.prob_hi[d] - in.prob_lo[d]) / static_cast<amrex::Real>(in.n_cell[d]);
        const amrex::Real first = in.prob_lo[d] + amrex::Real(0.5) * dx;
        const amrex::Real last = in.prob_lo[d]
                                 + (static_cast<amrex::Real>(in.n_cell[d] - 1) + amrex::Real(0.5)) * dx;
        if (!(x > first && x <= last)) {
            std::string msg = "flume.init.dam_position must leave a cell centre on each side of the"
                              " dam: it must be > ";
            msg += all_digits(first);
            msg += " and <= ";
            msg += all_digits(last);
            msg += " in ";
            msg += dir_name[d];
            msg += " (flume.prob_lo, flume.prob_hi and flume.n_cell set these); it is ";
            msg += all_digits(x);
            amrex::Abort(msg);
        }
    }

    return in;
}

} // namespace flume
