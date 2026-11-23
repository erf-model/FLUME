#include "FLUME_WaterBalance.H"
#include "FLUME_State.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace flume {

double
cell_area (const amrex::Geometry& geom)
{
    return static_cast<double>(geom.CellSize(0)) * static_cast<double>(geom.CellSize(1));
}

double
local_sum (const amrex::MultiFab& a)
{
    const auto& x = a.const_arrays();
    return amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{}, amrex::TypeList<double>{},
                            a, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<double>
        {
            return {static_cast<double>(x[box_no](i, j, k, 0))};
        });
}

WaterBalance::WaterBalance (const Inputs& in, const amrex::Geometry& geom)
    : m_geom(geom), m_cell_area(cell_area(geom)), m_on(in.log_int > 0)
{
    for (const auto& s : in.bc_segments) {
        if (s.type != BoundaryType::outflow) { continue; }
        m_segments.push_back(Segment{s.name, s.direction, s.high, s.cell_lo, s.cell_hi});
    }
    m_out.assign(m_segments.size(), 0.0);
    m_out_carry.assign(m_segments.size(), 0.0);
}

void
WaterBalance::add_outflow (const amrex::Array<amrex::MultiFab, 2>& flux, double dt)
{
    if (!m_on) { return; }
    const amrex::Box& domain = m_geom.Domain();
    for (std::size_t n = 0; n < m_segments.size(); ++n) {
        const Segment& s = m_segments[n];
        const int t = 1 - s.d;  // along the side
        // The faces of the segment: across direction d, at the low edge of
        // the domain or one past its last cell.
        const int face = (s.high == 1) ? domain.bigEnd(s.d) + 1 : domain.smallEnd(s.d);
        amrex::IntVect lo(0, 0, 0);
        amrex::IntVect hi(0, 0, 0);
        lo[s.d] = face;
        hi[s.d] = face;
        lo[t] = s.lo;
        hi[t] = s.hi;
        const amrex::Box faces(lo, hi, amrex::IntVect::TheDimensionVector(s.d));
        // Each boundary face belongs to one box only.
        amrex::ReduceOps<amrex::ReduceOpSum> op;
        amrex::ReduceData<double> data(op);
        for (amrex::MFIter mfi(flux[s.d]); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox() & faces;
            if (!b.ok()) { continue; }
            const auto& f = flux[s.d].const_array(mfi);
            op.eval(b, data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<double>
                {
                    return {static_cast<double>(f(i, j, k, State::h))};
                });
        }
        const double sum = amrex::get<0>(data.value(op));
        // A positive mass flux runs towards higher coordinates: out of the
        // high side, into the low side.
        const double sign = (s.high == 1) ? 1.0 : -1.0;
        compensated_add(m_out[n], m_out_carry[n],
                        sign * sum * static_cast<double>(m_geom.CellSize(t)) * dt);
    }
}

std::vector<double>
WaterBalance::totals (const amrex::MultiFab& state,
                      const std::vector<std::unique_ptr<Source>>& sources) const
{
    std::vector<double> t;
    t.push_back(volume(m_cell_area, local_sum(state)));  // State::h is component 0
    for (const auto& s : sources) {
        if (!s->balance_name().empty()) { t.push_back(s->balance_volume()); }
    }
    t.insert(t.end(), m_out.begin(), m_out.end());
    amrex::ParallelAllReduce::Sum(t.data(), static_cast<int>(t.size()),
                                  amrex::ParallelDescriptor::Communicator());
    return t;
}

void
WaterBalance::start (const amrex::MultiFab& state,
                     const std::vector<std::unique_ptr<Source>>& sources)
{
    if (!m_on) { return; }
    static_assert(State::h == 0, "local_sum reads component 0 of the state");
    m_source_names.clear();
    m_source_removes.clear();
    for (const auto& s : sources) {
        if (s->balance_name().empty()) { continue; }
        m_source_names.push_back(s->balance_name());
        m_source_removes.push_back(s->balance_removes());
    }
    const std::vector<double> t = totals(state, sources);
    m_v0 = t[0];
    if (amrex::ParallelDescriptor::IOProcessor()) {
        m_file.open("water_balance.txt");
        if (!m_file) { amrex::Abort("FLUME: cannot open water_balance.txt for writing"); }
        m_file << "# FLUME water balance: volumes in m3, cumulative from step 0\n"
               << "# epsilon " << std::setprecision(17)
               << std::numeric_limits<amrex::Real>::epsilon() << " (machine epsilon of the"
               << (sizeof(amrex::Real) == sizeof(double) ? " double" : " single")
               << "-precision state)\n"
               << "# residual = stored - stored at step 0 - (added) + (removed) + (out)\n"
               << "# added:";
        for (std::size_t n = 0; n < m_source_names.size(); ++n) {
            if (!m_source_removes[n]) { m_file << " " << m_source_names[n]; }
        }
        m_file << "\n# removed:";
        for (std::size_t n = 0; n < m_source_names.size(); ++n) {
            if (m_source_removes[n]) { m_file << " " << m_source_names[n]; }
        }
        m_file << "\n# out:";
        for (const auto& s : m_segments) { m_file << " out_" << s.name; }
        m_file << "\n# step time stored";
        for (std::size_t n = 0; n < m_source_names.size(); ++n) {
            m_file << " " << m_source_names[n];
        }
        for (const auto& s : m_segments) { m_file << " out_" << s.name; }
        m_file << " residual\n";
    }
    write_line(0, 0.0, t);
}

void
WaterBalance::log (const amrex::MultiFab& state,
                   const std::vector<std::unique_ptr<Source>>& sources, int step, double time)
{
    if (!m_on) { return; }
    write_line(step, time, totals(state, sources));
}

void
WaterBalance::write_line (int step, double time, const std::vector<double>& t)
{
    const std::size_t ns = m_source_names.size();
    // An infinite total (a state that overflows; the time step check aborts
    // such runs at step 0 with its own message) would make inf - inf, which
    // traps with amrex.fpe_trap_invalid = 1 before that message. The
    // residual is then NaN without the subtraction.
    const bool finite = std::isfinite(m_v0)
        && std::all_of(t.begin(), t.end(), [] (double v) { return std::isfinite(v); });
    double r = std::numeric_limits<double>::quiet_NaN();
    if (finite) {
        r = t[0] - m_v0;
        for (std::size_t n = 0; n < ns; ++n) { r -= t[1 + n]; }
        for (std::size_t n = 0; n < m_segments.size(); ++n) { r += t[1 + ns + n]; }
    }
    m_residual = r;
    m_last = t;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        m_file << step << " " << std::setprecision(17) << time << " " << t[0];
        for (std::size_t n = 0; n < ns; ++n) {
            m_file << " " << (m_source_removes[n] ? -t[1 + n] : t[1 + n]);
        }
        for (std::size_t n = 0; n < m_segments.size(); ++n) { m_file << " " << t[1 + ns + n]; }
        m_file << " " << r << "\n";
        m_file.flush();
    }
}

void
WaterBalance::print_summary () const
{
    if (!m_on || m_last.empty()) { return; }
    const std::size_t ns = m_source_names.size();
    std::ostringstream s;
    s << std::setprecision(10) << "FLUME: water balance (m3): stored " << m_last[0]
      << " (at step 0 " << m_v0 << ")";
    for (std::size_t n = 0; n < ns; ++n) {
        s << ", " << m_source_names[n] << " "
          << (m_source_removes[n] ? -m_last[1 + n] : m_last[1 + n]);
    }
    for (std::size_t n = 0; n < m_segments.size(); ++n) {
        s << ", out " << m_segments[n].name << " " << m_last[1 + ns + n];
    }
    s << std::setprecision(3) << ", residual " << m_residual << "\n";
    amrex::Print() << s.str();
}

} // namespace flume
