#include "FLUME_Drains.H"
#include "FLUME_Gauges.H"
#include "FLUME_RemoveWater.H"
#include "FLUME_State.H"
#include "FLUME_WaterBalance.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFIter.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <system_error>
#include <utility>
#include <vector>

namespace flume {

namespace {

//! The whole token as a finite decimal double, or false. A hexadecimal
//! number (0x10) is not one; a number too small for a double is its nearest
//! (0 or a denormal), one too large is not finite.
bool
to_double (const std::string& t, double& v)
{
    if (t.empty() || t.find_first_of("xX") != std::string::npos) { return false; }
    char* end = nullptr;
    v = std::strtod(t.c_str(), &end);
    return end == t.c_str() + t.size() && std::isfinite(v);
}

//! A double as text with every digit.
std::string
all_digits (double v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
    return s.str();
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

//! Parse a list of inlet locations on the I/O rank. Returns the message of
//! the first failed check, or an empty string.
std::string
parse_inlets (const std::string& key, const std::string& path, const amrex::Geometry& geom,
              std::vector<double>& x, std::vector<double>& y, std::vector<int>& line_of)
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
    const char* const dir_name[2] = {"x", "y"};
    std::string line;
    int line_no = 0;
    while (std::getline(f, line)) {
        if (line_no == std::numeric_limits<int>::max()) { return w + ": the file holds too many lines"; }
        ++line_no;
        // The byte-order mark some editors put before the first character.
        if (line_no == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) { line.erase(0, 3); }
        const std::string where = w + ", line " + std::to_string(line_no);
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) { line.erase(hash); }
        std::istringstream in(line);
        std::vector<std::string> tok;
        std::string t;
        while (in >> t) { tok.push_back(t); }
        if (tok.empty()) { continue; }
        if (tok.size() != 2) {
            return where + ": a line holds the x and the y of an inlet in m; it holds "
                   + std::to_string(tok.size()) + " item(s)";
        }
        double v[2] = {0.0, 0.0};
        for (int d = 0; d < 2; ++d) {
            if (!to_double(tok[d], v[d])) {
                return where + ": " + dir_name[d] + " must be a finite number; it is \"" + tok[d]
                       + "\"";
            }
            const double lo = static_cast<double>(geom.ProbLo(d));
            const double hi = static_cast<double>(geom.ProbHi(d));
            if (!(v[d] >= lo && v[d] <= hi)) {
                return where + ": the inlet must be inside the domain: from " + all_digits(lo)
                       + " to " + all_digits(hi) + " in " + dir_name[d]
                       + " (flume.prob_lo and flume.prob_hi); it is " + all_digits(v[d]);
            }
        }
        x.push_back(v[0]);
        y.push_back(v[1]);
        line_of.push_back(line_no);
    }
    if (f.bad()) { return w + ": the file cannot be read"; }
    if (x.empty()) {
        return w + ": the file holds no inlet: it needs at least one line with an x and a y in m";
    }
    return {};
}

//! One step of the 64-bit FNV-1a hash over the 4 bytes of `v`.
std::uint64_t
fnv1a (std::uint64_t hash, int v)
{
    auto u = static_cast<std::uint32_t>(v);
    for (int b = 0; b < 4; ++b) {
        hash ^= static_cast<std::uint64_t>(u & 0xffU);
        hash *= 1099511628211ULL;
        u >>= 8;
    }
    return hash;
}

//! An inlet cell: the entry of in.drains it belongs to, the line of the
//! file that placed it first (0 for an inlet of a location), and the number
//! of inlets it holds.
struct InletCell
{
    amrex::IntVect cell;
    std::size_t group = 0;
    int line = 0;
    int count = 0;
};

//! The rates of an inlet cell (i, j): a and b of the weir and orifice law,
//! or r and alpha of the approach-flow law; the others are 0.
struct InletRate
{
    int i;
    int j;
    amrex::Real weir;
    amrex::Real orifice;
    amrex::Real rate;
    amrex::Real exponent;
};

//! Write the rates of the inlets of one box into its cells. A free function:
//! nvcc does not accept a device lambda in a constructor.
void
set_inlet_rates (const std::vector<InletRate>& host, const amrex::Array4<amrex::Real>& wa,
                 const amrex::Array4<amrex::Real>& oa, const amrex::Array4<amrex::Real>& ra,
                 const amrex::Array4<amrex::Real>& ea)
{
    amrex::Gpu::DeviceVector<InletRate> device(host.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
    const InletRate* const p = device.data();
    amrex::ParallelFor(static_cast<long>(host.size()),
        [=] AMREX_GPU_DEVICE (long m) noexcept
        {
            const InletRate& v = p[m];
            wa(v.i, v.j, 0) = v.weir;
            oa(v.i, v.j, 0) = v.orifice;
            ra(v.i, v.j, 0) = v.rate;
            ea(v.i, v.j, 0) = v.exponent;
        });
    // Before `device` goes out of scope.
    amrex::Gpu::streamSynchronize();
}

} // namespace

void
read_inlet_file (const std::string& key, const std::string& path, const amrex::Geometry& geom,
                 std::vector<double>& x, std::vector<double>& y, std::vector<int>& line)
{
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    x.clear();
    y.clear();
    line.clear();
    std::string error;
    if (amrex::ParallelDescriptor::IOProcessor()) { error = parse_inlets(key, path, geom, x, y, line); }
    bcast_string(error, io);
    if (!error.empty()) { amrex::Abort(error); }
    long long n = static_cast<long long>(x.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, io);
    const auto un = static_cast<std::size_t>(n);
    x.resize(un);
    y.resize(un);
    line.resize(un);
    amrex::ParallelDescriptor::Bcast(x.data(), un, io);
    amrex::ParallelDescriptor::Bcast(y.data(), un, io);
    amrex::ParallelDescriptor::Bcast(line.data(), un, io);
}

Drains::Drains (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings)
    : m_weir(ba, dm, 1, 0), m_orifice(ba, dm, 1, 0), m_rate(ba, dm, 1, 0),
      m_exponent(ba, dm, 1, 0), m_drained(ba, dm, 1, 0), m_drained_carry(ba, dm, 1, 0),
      m_dry_depth(in.dry_depth), m_cell_area(cell_area(geom))
{
    m_weir.setVal(amrex::Real(0.0));
    m_orifice.setVal(amrex::Real(0.0));
    m_rate.setVal(amrex::Real(0.0));
    m_exponent.setVal(amrex::Real(0.0));
    m_drained.setVal(amrex::Real(0.0));
    m_drained_carry.setVal(amrex::Real(0.0));
    const double gravity = static_cast<double>(in.gravity);
    const double root_2g = std::sqrt(2.0 * gravity);
    // The inlet cells in the order they are placed, and where each is.
    std::vector<InletCell> cells;
    std::map<std::pair<int, int>, std::size_t> at;
    // The inputs are amrex::Real: at most its own digits10 (and the 10 of
    // the default stream), so a single-precision build prints 1.66 for an
    // input of 1.66 rather than its float, 1.659999967.
    const int digits = std::min(10, std::numeric_limits<amrex::Real>::digits10);
    for (std::size_t n = 0; n < in.drains.size(); ++n) {
        const DrainInlet& d = in.drains[n];
        const bool approach = (d.type == DrainType::approach_flow);
        (approach ? m_any_approach : m_any_weir) = true;
        // The parameters as the checkpoint keeps them, and as they are printed.
        std::ostringstream par;
        par << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
        std::ostringstream text;
        text << std::setprecision(digits);
        if (approach) {
            par << "approach_flow:" << d.area << ':' << d.coefficient << ':' << d.froude_exponent;
            text << "approach flow, area " << d.area << " m2, coefficient " << d.coefficient
                 << ", Froude exponent " << d.froude_exponent;
        } else {
            // The depth at which the weir and the orifice rates are equal
            // (in a single-precision build the kernel's, from a and b
            // rounded to amrex::Real, differs from this one by about 1e-7).
            const double ds = static_cast<double>(d.orifice_coefficient)
                              * static_cast<double>(d.area) * root_2g
                              / (static_cast<double>(d.weir_coefficient)
                                 * static_cast<double>(d.perimeter));
            par << d.perimeter << ':' << d.area << ':' << d.weir_coefficient << ':'
                << d.orifice_coefficient;
            text << "perimeter " << d.perimeter << " m, open area " << d.area
                 << " m2, weir coefficient " << d.weir_coefficient << ", orifice coefficient "
                 << d.orifice_coefficient << "; orifice above a depth of " << ds << " m";
        }
        if (d.file.empty()) {
            const std::string key = "flume.drains." + d.name + ".location";
            const amrex::IntVect iv = gauge_cell(d.location, geom);
            const auto it = at.find({iv[0], iv[1]});
            if (it != at.end()) {
                amrex::Abort(key + " puts inlet " + d.name + " in cell ("
                             + std::to_string(iv[0]) + ", " + std::to_string(iv[1])
                             + "), which holds inlet " + in.drains[cells[it->second].group].name
                             + " already: a cell holds the inlets of one name (merge two"
                               " inlets into one, or list them in one file)");
            }
            at[{iv[0], iv[1]}] = cells.size();
            cells.push_back({iv, n, 0, 1});
            std::ostringstream w;
            w << d.name << ':' << iv[0] << ':' << iv[1] << ':' << par.str();
            m_words.push_back(w.str());
            amrex::Print().SetPrecision(digits)
                << "FLUME: drain inlet " << d.name << " at (" << d.location[0] << ", "
                << d.location[1] << ") m in cell (" << iv[0] << ", " << iv[1]
                << "): " << text.str() << "\n";
            continue;
        }
        // The inlets of a list file: the cells they share hold their number.
        const std::string key = "flume.drains." + d.name + ".file";
        std::vector<double> x;
        std::vector<double> y;
        std::vector<int> line;
        read_inlet_file(key, d.file, geom, x, y, line);
        const std::size_t first = cells.size();
        for (std::size_t m = 0; m < x.size(); ++m) {
            const amrex::IntVect iv = gauge_cell({static_cast<amrex::Real>(x[m]),
                                                  static_cast<amrex::Real>(y[m])}, geom);
            const auto it = at.find({iv[0], iv[1]});
            if (it == at.end()) {
                at[{iv[0], iv[1]}] = cells.size();
                cells.push_back({iv, n, line[m], 1});
            } else if (cells[it->second].group == n) {
                ++cells[it->second].count;
            } else {
                amrex::Abort(key + " (" + d.file + "), line " + std::to_string(line[m])
                             + " puts an inlet in cell (" + std::to_string(iv[0]) + ", "
                             + std::to_string(iv[1]) + "), which holds inlet "
                             + in.drains[cells[it->second].group].name
                             + " already: only the inlets of one file share a cell");
            }
        }
        // The inlets of a cell fit in it, as one inlet does.
        for (std::size_t m = first; m < cells.size(); ++m) {
            const double sum = static_cast<double>(cells[m].count) * static_cast<double>(d.area);
            if (sum > m_cell_area) {
                amrex::Abort(key + " (" + d.file + "), line " + std::to_string(cells[m].line)
                             + ": cell (" + std::to_string(cells[m].cell[0]) + ", "
                             + std::to_string(cells[m].cell[1]) + ") holds "
                             + std::to_string(cells[m].count) + " inlets, whose areas (flume.drains."
                             + d.name + ".area) sum to " + all_digits(sum)
                             + " m2: more than the area of a cell (" + all_digits(m_cell_area)
                             + " m2)");
            }
        }
        // The cells and their counts, in the order of the cells.
        std::uint64_t hash = 14695981039346656037ULL;
        long long shared = 0;
        for (const auto& [ij, m] : at) {
            if (cells[m].group != n) { continue; }
            hash = fnv1a(fnv1a(fnv1a(hash, ij.first), ij.second), cells[m].count);
            if (cells[m].count > 1) { ++shared; }
        }
        std::ostringstream w;
        w << d.name << ":file:" << x.size() << ':' << (cells.size() - first) << ':' << std::hex
          << hash << std::dec << ':' << par.str();
        m_words.push_back(w.str());
        amrex::Print() << "FLUME: drain inlets " << d.name << ": " << x.size() << " from "
                       << d.file << " in " << (cells.size() - first) << " cell(s), " << shared
                       << " with more than one inlet: " << text.str() << "\n";
    }
    // The rates of each inlet cell: those of one inlet times the number.
    // They are binned by the box that holds the cell, and each box of this
    // rank is filled in one pass.
    std::map<int, std::vector<InletRate>> bins;
    std::vector<std::pair<int, amrex::Box>> hits;
    for (const InletCell& c : cells) {
        const DrainInlet& d = in.drains[c.group];
        const auto count = static_cast<double>(c.count);
        InletRate v{c.cell[0], c.cell[1], amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                    amrex::Real(0.0)};
        bool positive = false;
        if (d.type == DrainType::approach_flow) {
            const double alpha = static_cast<double>(d.froude_exponent);
            const double r = count * static_cast<double>(d.coefficient)
                             * static_cast<double>(d.area) * std::pow(gravity, -0.5 * alpha)
                             / m_cell_area;
            v.rate = static_cast<amrex::Real>(r);
            v.exponent = d.froude_exponent;
            positive = v.rate > amrex::Real(0.0);
        } else {
            const double a = count * (static_cast<double>(d.weir_coefficient)
                                      * static_cast<double>(d.perimeter) / m_cell_area);
            const double b = count * (static_cast<double>(d.orifice_coefficient)
                                      * static_cast<double>(d.area) * root_2g / m_cell_area);
            v.weir = static_cast<amrex::Real>(a);
            v.orifice = static_cast<amrex::Real>(b);
            positive = v.weir > amrex::Real(0.0) && v.orifice > amrex::Real(0.0);
        }
        // A rate of 0 would read as a cell without an inlet.
        if (!positive) {
            amrex::Abort("flume.drains." + d.name + ": the rate of the inlet in cell ("
                         + std::to_string(c.cell[0]) + ", " + std::to_string(c.cell[1])
                         + ") is 0 in the precision of this build: its area and coefficients"
                           " are too small for the cell");
        }
        ba.intersections(amrex::Box(c.cell, c.cell), hits);
        for (const auto& hit : hits) {
            if (dm[hit.first] == amrex::ParallelDescriptor::MyProc()) {
                bins[hit.first].push_back(v);
            }
        }
    }
    for (amrex::MFIter mfi(m_rate); mfi.isValid(); ++mfi) {
        const auto bin = bins.find(mfi.index());
        if (bin == bins.end()) { continue; }
        set_inlet_rates(bin->second, m_weir.array(mfi), m_orifice.array(mfi),
                        m_rate.array(mfi), m_exponent.array(mfi));
    }
    // No inlet in a building cell. One pass over the grid for all the
    // inlets (the weir rate or the approach rate is > 0 in the inlet cells
    // only); the cells are searched one by one only to name the inlet
    // before stopping.
    if (buildings != nullptr && !cells.empty()) {
        amrex::MultiFab both(m_weir.boxArray(), m_weir.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(both, m_weir, 0, 0, 1, 0);
        amrex::MultiFab::Add(both, m_rate, 0, 0, 1, 0);
        amrex::MultiFab::Multiply(both, *buildings, 0, 0, 1, 0);
        if (both.max(0, 0) > amrex::Real(0.0)) {
            for (const InletCell& c : cells) {
                if (buildings->max(amrex::Box(c.cell, c.cell), 0, 0) > amrex::Real(0.0)) {
                    const DrainInlet& d = in.drains[c.group];
                    if (d.file.empty()) {
                        amrex::Abort("flume.drains." + d.name + ".location puts inlet " + d.name
                                     + " in a building cell (flume.buildings.file)");
                    }
                    amrex::Abort("flume.drains." + d.name + ".file (" + d.file + "), line "
                                 + std::to_string(c.line) + " puts an inlet in a building cell"
                                   " (flume.buildings.file)");
                }
            }
        }
    }
}

void
Drains::apply (amrex::MultiFab& state, const SourceContext& /*ctx*/, amrex::Real dt)
{
    // A cell holds inlets of one law, so the two passes touch other cells.
    if (m_any_weir) {
        apply_drains(state, m_weir, m_orifice, m_drained, m_drained_carry, m_dry_depth, dt);
    }
    if (m_any_approach) {
        apply_approach_drains(state, m_rate, m_exponent, m_drained, m_drained_carry, m_dry_depth,
                              dt);
    }
}

double
Drains::balance_volume () const
{
    // As for infiltration: with the carry the volume follows the state to
    // round-off of a double, also in a single-precision build.
    return -(volume(m_cell_area, local_sum(m_drained))
             + volume(m_cell_area, local_sum(m_drained_carry)));
}

void
Drains::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_words(name() + "_inlets", m_words);
    chk.put_multifab(name() + "_depth", m_drained);
    chk.put_multifab(name() + "_depth_carry", m_drained_carry);
}

void
Drains::read_checkpoint (const CheckpointReader& chk)
{
    // The drained depth of a cell belongs to its inlet: the inlets must be
    // those the checkpoint was run with.
    if (chk.get_words(name() + "_inlets") != m_words) {
        std::string was;
        std::string now;
        for (const auto& w : chk.get_words(name() + "_inlets")) { was += (was.empty() ? "" : " ") + w; }
        for (const auto& w : m_words) { now += (now.empty() ? "" : " ") + w; }
        chk.fail("its drain inlets (" + was + ") differ from those of the inputs (" + now
                 + "): for an inlet of a location name:i:j, for the inlets of a file"
                   " name:file:inlets:cells:hash of the cells, then"
                   " perimeter:area:weir_coefficient:orifice_coefficient or"
                   " approach_flow:area:coefficient:froude_exponent, from flume.drains.*");
    }
    chk.get_multifab(name() + "_depth", m_drained);
    chk.get_multifab(name() + "_depth_carry", m_drained_carry);
}

void
apply_drains (amrex::MultiFab& state, const amrex::MultiFab& weir, const amrex::MultiFab& orifice,
              amrex::MultiFab& drained, amrex::MultiFab& drained_carry, amrex::Real dry_depth,
              amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& aa = weir.const_arrays();
    const auto& bb = orifice.const_arrays();
    const auto& ff = drained.arrays();
    const auto& fc = drained_carry.arrays();
    // Each cell reads and writes only itself. In double precision.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const double a = static_cast<double>(aa[box_no](i, j, k));
            if (!(a > 0.0)) { return; }  // no inlet
            const auto& c = s[box_no];
            const amrex::Real h = c(i, j, k, State::h);
            // No water, or a NaN depth (left for the check of the state; an
            // infinite one gives a NaN loss below and is left too).
            if (!(h > amrex::Real(0.0))) { return; }
            const double b = static_cast<double>(bb[box_no](i, j, k));
            const double w = static_cast<double>(h);
            const double h1 = drain_depth_after(w, a, b, static_cast<double>(dt));
            // drain_depth_after returns a depth in [0, w], so w - h1 is in
            // [0, w], and so is its rounding to amrex::Real, since w = h is
            // an amrex::Real: the loss never exceeds the water of the cell.
            const amrex::Real loss = static_cast<amrex::Real>(w - h1);
            if (!(loss > amrex::Real(0.0))) { return; }
            remove_water(c, ff[box_no](i, j, k), fc[box_no](i, j, k), i, j, k, h, loss, dry_depth);
        });
    amrex::Gpu::streamSynchronize();
}

void
apply_approach_drains (amrex::MultiFab& state, const amrex::MultiFab& rate,
                       const amrex::MultiFab& exponent, amrex::MultiFab& drained,
                       amrex::MultiFab& drained_carry, amrex::Real dry_depth, amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& rr = rate.const_arrays();
    const auto& ee = exponent.const_arrays();
    const auto& ff = drained.arrays();
    const auto& fc = drained_carry.arrays();
    // Each cell reads and writes only itself. In double precision.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const double r = static_cast<double>(rr[box_no](i, j, k));
            if (!(r > 0.0)) { return; }  // no inlet
            const auto& c = s[box_no];
            const amrex::Real h = c(i, j, k, State::h);
            // A dry cell holds no discharge; a NaN depth is left for the
            // check of the state.
            if (!(h > dry_depth) || !(h > amrex::Real(0.0))) { return; }
            const double w = static_cast<double>(h);
            const double u = std::hypot(static_cast<double>(c(i, j, k, State::hu)),
                                        static_cast<double>(c(i, j, k, State::hv))) / w;
            const double alpha = static_cast<double>(ee[box_no](i, j, k));
            // dh/dt = -cc h^(-alpha/2) at the speed u. With alpha = -1 the
            // speed drops out (u^0 = 1, also at rest).
            const double cc = r * std::pow(u, 1.0 + alpha);
            // Still water (alpha above -1), or a speed that is not a finite
            // number (left for the check of the state).
            if (!(cc > 0.0) || !(cc <= std::numeric_limits<double>::max())) { return; }
            double h1 = approach_depth_after(w, cc, -0.5 * alpha, static_cast<double>(dt));
            // A film at or below the dry depth holds no discharge, so this
            // law could never take it later: the inlet takes it now.
            if (!(h1 > static_cast<double>(dry_depth))) { h1 = 0.0; }
            // As in apply_drains: h1 is in [0, w], so the loss never exceeds
            // the water of the cell.
            const amrex::Real loss = static_cast<amrex::Real>(w - h1);
            if (!(loss > amrex::Real(0.0))) { return; }
            remove_water(c, ff[box_no](i, j, k), fc[box_no](i, j, k), i, j, k, h, loss, dry_depth);
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
