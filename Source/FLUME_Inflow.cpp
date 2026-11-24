#include "FLUME_Inflow.H"
#include "FLUME_State.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <limits>
#include <utility>

namespace flume {

namespace {

/**
 * The speed q / h_b at the face. h_b >= h_c > 0 when q > 0; with q = 0 the
 * face may be dry (h_b = 0), and the speed is 0. The denominator is kept at
 * or above `tiny` (> 0) rather than the division skipped: compilers evaluate
 * both sides of a select, and 0 / 0 would raise the invalid-operation flag.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
face_speed (amrex::Real q, amrex::Real hb, amrex::Real tiny) noexcept
{
    return q / amrex::max(hb, tiny);
}

//! The smallest normal number of the build, for face_speed; taken on the
//! host and passed to the kernels by value.
constexpr amrex::Real tiny_depth = std::numeric_limits<amrex::Real>::min();

//! The boundary faces of the cells lo to hi along the side, on the low
//! (high = 0) or high side across direction d.
amrex::Box
segment_faces (const amrex::Box& domain, int d, int high, int lo, int hi)
{
    const int t = 1 - d;
    const int face = (high == 1) ? domain.bigEnd(d) + 1 : domain.smallEnd(d);
    amrex::IntVect a(0, 0, 0);
    amrex::IntVect b(0, 0, 0);
    a[d] = face;
    b[d] = face;
    a[t] = lo;
    b[t] = hi;
    return amrex::Box(a, b, amrex::IntVect::TheDimensionVector(d));
}

/**
 * The fluxes of the boundary faces `faces` (across direction d, on the
 * high side when high = 1) in the fluxes `flux` of every box: mass flux q
 * into the domain, momentum flux q u_b + g h_b^2 / 2 along the normal, zero
 * along the side. The cell inside a face is in the same box (the face is on
 * the edge of the domain).
 *
 * The face takes the cell average h, not a reconstructed depth. That keeps
 * still water at rest over any bed only because the ghost cells copy the
 * depth and the bed: the reconstruction then has no slope in the boundary
 * cell, and with q = 0 the pressure flux g h^2 / 2 cancels the bed source
 * of the cell. A ghost bed that continues a slope would need the face depth
 * to change with it.
 */
void
set_face_fluxes (amrex::MultiFab& flux, const amrex::MultiFab& state, const amrex::Box& faces,
                 int d, int high, amrex::Real q, amrex::Real g, amrex::Real dry)
{
    // From a face to the cell inside it, and the sign of the normal into the
    // domain along +d.
    const amrex::IntVect to_cell = (high == 1) ? -amrex::IntVect::TheDimensionVector(d)
                                               : amrex::IntVect(0, 0, 0);
    const amrex::Real into = (high == 1) ? amrex::Real(-1.0) : amrex::Real(1.0);
    const int normal = State::hu + d;
    const int along = State::hu + (1 - d);
    const amrex::Real tiny = tiny_depth;
    for (amrex::MFIter mfi(flux); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox() & faces;
        if (!b.ok()) { continue; }
        const auto& f = flux.array(mfi);
        const auto& s = state.const_array(mfi);
        amrex::ParallelFor(b,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::IntVect c = amrex::IntVect(i, j, k) + to_cell;
                amrex::Real h = s(c, State::h);
                amrex::Real un = amrex::Real(0.0);
                if (h > dry) {
                    un = into * s(c, normal) / h;
                } else {
                    h = amrex::Real(0.0);
                }
                const amrex::Real hb = inflow_depth(q, un, h, g);
                const amrex::Real ub = face_speed(q, hb, tiny);
                f(i, j, k, State::h) = into * q;
                f(i, j, k, normal) = q * ub + amrex::Real(0.5) * g * hb * hb;
                f(i, j, k, along) = amrex::Real(0.0);
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! The largest (u_b + c_b) * dni + c_b * dti over the faces `faces`, on this
//! rank, with the discharge q per unit width (as set_face_fluxes).
amrex::Real
face_rate (const amrex::MultiFab& state, const amrex::Box& faces, int d, int high, amrex::Real q,
           amrex::Real g, amrex::Real dry, amrex::Real dni, amrex::Real dti)
{
    const amrex::IntVect to_cell = (high == 1) ? -amrex::IntVect::TheDimensionVector(d)
                                               : amrex::IntVect(0, 0, 0);
    const amrex::Real into = (high == 1) ? amrex::Real(-1.0) : amrex::Real(1.0);
    const int normal = State::hu + d;
    // The cells inside the faces: the same indices, shifted by one below a
    // face on the high side.
    const amrex::Box cells = amrex::Box(faces.smallEnd(), faces.bigEnd()).shift(to_cell);
    const amrex::Real tiny = tiny_depth;
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real> data(op);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox() & cells;
        if (!b.ok()) { continue; }
        const auto& s = state.const_array(mfi);
        op.eval(b, data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real>
            {
                amrex::Real h = s(i, j, k, State::h);
                amrex::Real un = amrex::Real(0.0);
                if (h > dry) {
                    un = into * s(i, j, k, normal) / h;
                } else {
                    h = amrex::Real(0.0);
                }
                const amrex::Real hb = inflow_depth(q, un, h, g);
                const amrex::Real ub = face_speed(q, hb, tiny);
                const amrex::Real cb = std::sqrt(g * hb);
                return {(ub + cb) * dni + cb * dti};
            });
    }
    return amrex::get<0>(data.value(op));
}

const char*
side_name (int d, int high)
{
    const char* const n[2][2] = {{"xlo", "xhi"}, {"ylo", "yhi"}};
    return n[d][high];
}

} // namespace

Inflow::Inflow (const Inputs& in, const amrex::Geometry& geom)
    : m_geom(geom), m_gravity(in.gravity), m_dry_depth(in.dry_depth)
{
    for (const auto& s : in.bc_segments) {
        if (s.type != BoundaryType::discharge) { continue; }
        std::vector<double> time;
        std::vector<double> q;
        read_time_table("flume.bc." + s.name + ".file", s.file, "discharge", "m3/s", time, q);
        const double q_max = *std::max_element(q.begin(), q.end());
        if (time.front() > 0.0 && q.front() > 0.0) {
            amrex::Print() << "FLUME: inflow " << s.name << ": note: the first time of " << s.file
                           << " is " << time.front() << " s; its discharge, " << q.front()
                           << " m3/s, enters from the start of the run (a line \"0 0\" first"
                              " starts the hydrograph at " << time.front() << " s)\n";
        }
        const int t = 1 - s.direction;
        const double width = static_cast<double>(s.cell_hi - s.cell_lo + 1)
                             * static_cast<double>(geom.CellSize(t));
        const bool linear = s.interpolation == RainInterpolation::linear;
        amrex::Print() << "FLUME: inflow " << s.name << ": hydrograph from " << s.file << ", "
                       << time.size() << " entries, " << (linear ? "linear" : "step")
                       << " in time, up to " << q_max << " m3/s over " << width << " m of "
                       << side_name(s.direction, s.high) << "\n";
        m_segments.push_back(Segment{s.name, s.direction, s.high, s.cell_lo, s.cell_hi, width,
                                     TimeTable(std::move(time), std::move(q), s.interpolation)});
    }
}

void
Inflow::set_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                    double t0, double t1, double dt) const
{
    for (const auto& s : m_segments) {
        // Never below 0: the integral can only grow, but its two values carry
        // their own rounding.
        const double volume = std::max(s.q.integral(t1) - s.q.integral(t0), 0.0);
        const auto q = static_cast<amrex::Real>(volume / (dt * s.width));
        set_face_fluxes(flux[s.d], state, segment_faces(m_geom.Domain(), s.d, s.high, s.lo, s.hi),
                        s.d, s.high, q, m_gravity, m_dry_depth);
    }
}

amrex::Real
Inflow::max_rate (const amrex::MultiFab& state, double t0, double t1) const
{
    amrex::Real rate = amrex::Real(0.0);
    for (const auto& s : m_segments) {
        const auto q = static_cast<amrex::Real>(s.q.max_value(t0, t1) / s.width);
        rate = std::max(rate, face_rate(state, segment_faces(m_geom.Domain(), s.d, s.high, s.lo,
                                                             s.hi),
                                        s.d, s.high, q, m_gravity, m_dry_depth,
                                        m_geom.InvCellSize(s.d), m_geom.InvCellSize(1 - s.d)));
    }
    if (!m_segments.empty()) {
        amrex::ParallelAllReduce::Max(rate, amrex::ParallelDescriptor::Communicator());
    }
    return rate;
}

} // namespace flume
