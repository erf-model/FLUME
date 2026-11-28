#include "FLUME_InflowSource.H"
#include "FLUME_Gauges.H"
#include "FLUME_State.H"
#include "FLUME_WaterBalance.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace flume {

namespace {

//! The cell that holds x along direction d (the coordinate of a piece's
//! middle, inside the domain or on its edge). A coordinate within rounding
//! of a face is on the face and goes to the cell above it, by the rule of
//! gauge_cell (FLUME_Gauges.cpp), so a line along a face takes the cell a
//! point on that face takes.
int
cell_of (double x, int d, const amrex::Geometry& geom)
{
    const double lo = static_cast<double>(geom.ProbLo(d));
    const double hi = static_cast<double>(geom.ProbHi(d));
    const int n = geom.Domain().length(d);
    const double f = (x - lo) / (hi - lo) * static_cast<double>(n);
    const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    const double big = std::max({std::abs(x), std::abs(lo), std::abs(hi)});
    const double tol = std::min(4.0 * eps * big / (hi - lo) * static_cast<double>(n), 1e-3);
    const double r = std::round(f);
    const int c = (std::abs(f - r) <= tol) ? static_cast<int>(r) : static_cast<int>(std::floor(f));
    return std::clamp(c, 0, n - 1);
}

//! A number with every digit of amrex::Real, for the checkpoint words.
std::string
digits (amrex::Real v)
{
    std::ostringstream s;
    s << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << v;
    return s.str();
}

//! hit[m] = 1 where cell e[m] is a building cell (height above 0), else 0.
void
mark_building_cells (const InflowCell* e, int n, const amrex::Array4<const amrex::Real>& b, int* hit)
{
    amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (int m) noexcept
    {
        hit[m] = (b(e[m].i, e[m].j, 0) > amrex::Real(0.0)) ? 1 : 0;
    });
}

//! Add the depth e[m].fraction * scale, m, to each of the n cells e (all
//! different) with a compensated sum, as add_rain_depth: the depth a step
//! rounds away is kept in the carry c and added with the next step's, and a
//! negative carry waits while the cell holds less. Sums the depth asked for
//! into `data`.
void
add_inflow_depth (amrex::ReduceOps<amrex::ReduceOpSum>& op, amrex::ReduceData<double>& data,
                  const InflowCell* e, int n, double scale, const amrex::Array4<amrex::Real>& s,
                  const amrex::Array4<amrex::Real>& c)
{
    op.eval(n, data, [=] AMREX_GPU_DEVICE (int m) noexcept -> amrex::GpuTuple<double>
    {
        const int i = e[m].i;
        const int j = e[m].j;
        const auto dr = static_cast<amrex::Real>(e[m].fraction * scale);
        const amrex::Real h = s(i, j, 0, State::h);
        const amrex::Real y = dr + c(i, j, 0);
        if (!(h + y >= amrex::Real(0.0))) {
            c(i, j, 0) = y;
            return {static_cast<double>(dr)};
        }
        const amrex::Real h_new = h + y;
        c(i, j, 0) = y - (h_new - h);
        s(i, j, 0, State::h) = h_new;
        return {static_cast<double>(dr)};
    });
}

} // namespace

//! The shortest piece of a line, as a part of its length, that line_cells
//! gives to a cell.
constexpr double kSliver = 1e-12;

std::vector<InflowCell>
line_cells (double x0, double y0, double x1, double y1, const amrex::Geometry& geom)
{
    const double p0[2] = {x0, y0};
    const double p1[2] = {x1, y1};
    // The parameters t in (0, 1) where the line crosses a line of the grid,
    // with the ends.
    std::vector<double> t{0.0, 1.0};
    for (int d = 0; d < 2; ++d) {
        const double span = p1[d] - p0[d];
        if (span == 0.0) { continue; }
        const double lo = static_cast<double>(geom.ProbLo(d));
        const double hi = static_cast<double>(geom.ProbHi(d));
        const int n = geom.Domain().length(d);
        for (int f = 1; f < n; ++f) {
            const double x = lo + (hi - lo) * static_cast<double>(f) / static_cast<double>(n);
            const double s = (x - p0[d]) / span;
            if (s > 0.0 && s < 1.0) { t.push_back(s); }
        }
    }
    std::sort(t.begin(), t.end());
    std::map<std::pair<int, int>, double> parts;
    double total = 0.0;
    for (std::size_t m = 0; m + 1 < t.size(); ++m) {
        const double w = t[m + 1] - t[m];
        // A piece shorter than kSliver of the line is rounding where the
        // line passes through a corner of the grid, or ends on a face: the
        // crossings of x and y differ by an ulp, and the middle of the sliver
        // lies in a cell the line only touches.
        if (!(w > kSliver)) { continue; }
        const double mid = 0.5 * (t[m] + t[m + 1]);
        const int i = cell_of(p0[0] + mid * (p1[0] - p0[0]), 0, geom);
        const int j = cell_of(p0[1] + mid * (p1[1] - p0[1]), 1, geom);
        parts[{i, j}] += w;
        total += w;
    }
    std::vector<InflowCell> cells;
    cells.reserve(parts.size());
    for (const auto& [ij, w] : parts) {
        cells.push_back(InflowCell{ij.first, ij.second, w / total});
    }
    return cells;
}

InflowSource::InflowSource (const InflowSpec& spec, const amrex::Geometry& geom,
                            const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                            const amrex::MultiFab* buildings)
    : m_name(spec.name),
      m_q([&spec] {
          std::vector<double> time;
          std::vector<double> q;
          read_time_table("flume.inflows." + spec.name + ".file", spec.file, "discharge", "m3/s",
                          time, q);
          return TimeTable(std::move(time), std::move(q), spec.interpolation);
      }()),
      m_carry(ba, dm, 1, 0), m_cell_area(cell_area(geom))
{
    m_carry.setVal(amrex::Real(0.0));
    const std::string key = "flume.inflows." + m_name + (spec.line ? ".line" : ".location");
    std::ostringstream w;
    w << m_name;
    if (spec.line) {
        m_cells = line_cells(static_cast<double>(spec.coords[0]), static_cast<double>(spec.coords[1]),
                             static_cast<double>(spec.coords[2]), static_cast<double>(spec.coords[3]),
                             geom);
        w << ":line";
        for (int c = 0; c < 4; ++c) { w << ':' << digits(spec.coords[c]); }
    } else {
        const amrex::IntVect iv = gauge_cell({spec.coords[0], spec.coords[1]}, geom);
        m_cells.push_back(InflowCell{iv[0], iv[1], 1.0});
        w << ":point:" << iv[0] << ':' << iv[1];
    }
    m_words.push_back(w.str());
    for (const auto& c : m_cells) { m_fraction_max = std::max(m_fraction_max, c.fraction); }

    // The cells of each local box, on the device.
    std::vector<int> owner(m_cells.size(), -1);  // local index of the box, or -1
    m_local.resize(m_carry.local_size());
    // The host lists live until the copies are done (streamSynchronize).
    std::vector<std::vector<InflowCell>> host(m_local.size());
    for (amrex::MFIter mfi(m_carry); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto& mine = host[mfi.LocalIndex()];
        for (std::size_t n = 0; n < m_cells.size(); ++n) {
            if (bx.contains(amrex::IntVect(m_cells[n].i, m_cells[n].j, 0))) {
                mine.push_back(m_cells[n]);
                owner[n] = mfi.LocalIndex();
            }
        }
        auto& dv = m_local[mfi.LocalIndex()];
        dv.resize(mine.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, mine.begin(), mine.end(), dv.begin());
    }
    amrex::Gpu::streamSynchronize();

    // No cell in a building. The building heights of this rank's cells of
    // the inflow, then the first building cell over the ranks.
    if (buildings != nullptr) {
        int first = static_cast<int>(m_cells.size());
        for (amrex::MFIter mfi(*buildings); mfi.isValid(); ++mfi) {
            const auto& dv = m_local[mfi.LocalIndex()];
            const int n = static_cast<int>(dv.size());
            if (n == 0) { continue; }
            amrex::Gpu::DeviceVector<int> hit(n);
            mark_building_cells(dv.data(), n, buildings->const_array(mfi), hit.data());
            std::vector<int> h(n);
            amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, hit.begin(), hit.end(), h.begin());
            amrex::Gpu::streamSynchronize();
            // The entries of m_local follow m_cells: the index in m_cells of
            // each building cell of this box.
            int k = 0;
            for (std::size_t q = 0; q < m_cells.size(); ++q) {
                if (owner[q] != mfi.LocalIndex()) { continue; }
                if (h[k] != 0) { first = std::min(first, static_cast<int>(q)); }
                ++k;
            }
        }
        amrex::ParallelDescriptor::ReduceIntMin(first);
        if (first < static_cast<int>(m_cells.size())) {
            amrex::Abort(key + " puts inflow " + m_name + " in a building cell ("
                         + std::to_string(m_cells[first].i) + ", "
                         + std::to_string(m_cells[first].j) + ") (flume.buildings.file)");
        }
    }

    const auto& times = m_q.times();
    const auto& values = m_q.values();
    if (times.front() > 0.0 && values.front() > 0.0) {
        amrex::Print() << "FLUME: inflow " << m_name << ": note: the first time of " << spec.file
                       << " is " << times.front() << " s; its discharge, " << values.front()
                       << " m3/s, enters from the start of the run (a line \"0 0\" first"
                          " starts the hydrograph at " << times.front() << " s)\n";
    }
    const bool linear = spec.interpolation == RainInterpolation::linear;
    amrex::Print() << "FLUME: inflow " << m_name << ": hydrograph from " << spec.file << ", "
                   << times.size() << " entries, " << (linear ? "linear" : "step")
                   << " in time, up to " << *std::max_element(values.begin(), values.end())
                   << " m3/s into " << m_cells.size() << " cell(s)";
    // The coordinates to the digits amrex::Real holds, so that 1.3 reads 1.3
    // in a single-precision build too.
    std::ostringstream where;
    where << std::setprecision(std::numeric_limits<amrex::Real>::digits10);
    if (spec.line) {
        where << " along the line from (" << spec.coords[0] << ", " << spec.coords[1]
              << ") to (" << spec.coords[2] << ", " << spec.coords[3] << ") m\n";
    } else {
        where << ", (" << m_cells[0].i << ", " << m_cells[0].j << "), at ("
              << spec.coords[0] << ", " << spec.coords[1] << ") m\n";
    }
    amrex::Print() << where.str();
}

void
InflowSource::apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real /*dt*/)
{
    // The ends of the step in double precision, as for the rain. Never below
    // 0: the integral can only grow, but its two values carry their own
    // rounding.
    const double dv = std::max(m_q.integral(ctx.time_end) - m_q.integral(ctx.time), 0.0);
    if (!(dv > 0.0)) { return; }
    const double scale = dv / m_cell_area;  // m
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<double> data(op);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const auto& dvec = m_local[mfi.LocalIndex()];
        const int n = static_cast<int>(dvec.size());
        if (n == 0) { continue; }
        add_inflow_depth(op, data, dvec.data(), n, scale, state.array(mfi), m_carry.array(mfi));
    }
    compensated_add(m_volume, m_volume_carry,
                    volume(m_cell_area, amrex::get<0>(data.value(op))));
}

double
InflowSource::balance_volume () const
{
    return m_volume - volume(m_cell_area, local_sum(m_carry));
}

double
InflowSource::max_depth_rate (double t0, double t1) const
{
    return m_fraction_max * m_q.max_value(t0, t1) / m_cell_area;
}

void
InflowSource::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_words(name() + "_cells", m_words);
    chk.put_multifab(name() + "_carry", m_carry);
    chk.put_per_rank(name() + "/volume", {m_volume, m_volume_carry});
}

void
InflowSource::read_checkpoint (const CheckpointReader& chk)
{
    // The carry and the volume belong to the cells of the inflow.
    const std::vector<std::string> was = chk.get_words(name() + "_cells");
    if (was != m_words) {
        chk.fail("its inflow " + (was.empty() ? std::string("(none)") : was[0])
                 + " differs from that of the inputs (" + m_words[0]
                 + "): name:point:i:j or name:line:x0:y0:x1:y1, from flume.inflows." + m_name
                 + ".location or .line");
    }
    chk.get_multifab(name() + "_carry", m_carry);
    const std::vector<double> v = chk.get_per_rank(name() + "/volume", 2);
    m_volume = v[0];
    m_volume_carry = v[1];
}

} // namespace flume
