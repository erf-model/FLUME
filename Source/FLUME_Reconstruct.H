#ifndef FLUME_RECONSTRUCT_H_
#define FLUME_RECONSTRUCT_H_

#include "FLUME_State.H"

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * The minmod limiter: zero where `a` and `b` differ in sign (or one is zero),
 * otherwise the one of smaller magnitude.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
minmod (amrex::Real a, amrex::Real b) noexcept
{
    const amrex::Real zero = amrex::Real(0.0);
    if (a > zero && b > zero) { return (a < b) ? a : b; }
    if (a < zero && b < zero) { return (a > b) ? a : b; }
    return zero;
}

//! What a cell presents on its low and on its high side along one direction.
struct CellSides
{
    amrex::Real h_lo;  //!< depth on the low side, m
    amrex::Real h_hi;  //!< depth on the high side, m
    amrex::Real z_lo;  //!< bed on the low side minus the bed of the cell, m
    amrex::Real z_hi;  //!< bed on the high side minus the bed of the cell, m
    amrex::Real un_lo; //!< velocity along the direction on the low side, m/s
    amrex::Real un_hi; //!< and on the high side
    amrex::Real ut_lo; //!< velocity across the direction on the low side, m/s
    amrex::Real ut_hi; //!< and on the high side
    amrex::Real d_eta; //!< limited difference of the surface elevation over the cell, m
};

/**
 * Piecewise-linear reconstruction of cell `c` along the direction of the unit
 * vector `e` (Docs/sphinx_doc/theory/Discretization.rst). `a` is the state, `z` the bed;
 * `i_norm` and `i_tang` are the components of the discharge along and across
 * the direction. Reads the cells c - e, c and c + e.
 *
 * The depth, the surface elevation and the two velocities each get a minmod
 * difference from the two one-sided differences. The depth on a side is the
 * depth of the cell plus or minus half its difference, and is >= 0 when the
 * three depths are. The bed on a side follows from the surface and the depth,
 * so under a level surface the surface on both sides is level whatever the
 * bed does. Only differences of bed elevations enter, so no digits are lost
 * when the elevations are large.
 *
 * The velocity on a side is weighted with the depth on the other side, so
 * that the mean of depth times velocity over the two sides is the discharge
 * of the cell.
 *
 * A cell whose depth is at or below `dry_depth` is not reconstructed: both
 * sides hold its own depth and bed and no velocity. The velocity of a dry
 * neighbour counts as zero.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE CellSides
reconstruct (const amrex::Array4<const amrex::Real>& a, const amrex::Array4<const amrex::Real>& z,
             const amrex::IntVect& c, const amrex::IntVect& e, int i_norm, int i_tang,
             amrex::Real dry_depth) noexcept
{
    const amrex::Real zero = amrex::Real(0.0);
    const amrex::Real half = amrex::Real(0.5);
    const amrex::Real h0 = a(c, State::h);
    if (!(h0 > dry_depth)) { return {h0, h0, zero, zero, zero, zero, zero, zero, zero}; }

    const amrex::IntVect m = c - e;
    const amrex::IntVect p = c + e;
    const amrex::Real hm = a(m, State::h);
    const amrex::Real hp = a(p, State::h);
    const amrex::Real dh_m = h0 - hm;
    const amrex::Real dh_p = hp - h0;
    const amrex::Real dh = minmod(dh_m, dh_p);
    const amrex::Real d_eta = minmod(dh_m + (z(c) - z(m)), dh_p + (z(p) - z(c)));
    const amrex::Real dz = d_eta - dh;

    const amrex::Real un0 = a(c, i_norm) / h0;
    const amrex::Real ut0 = a(c, i_tang) / h0;
    const amrex::Real un_m = (hm > dry_depth) ? a(m, i_norm) / hm : zero;
    const amrex::Real ut_m = (hm > dry_depth) ? a(m, i_tang) / hm : zero;
    const amrex::Real un_p = (hp > dry_depth) ? a(p, i_norm) / hp : zero;
    const amrex::Real ut_p = (hp > dry_depth) ? a(p, i_tang) / hp : zero;
    const amrex::Real dun = minmod(un0 - un_m, un_p - un0);
    const amrex::Real dut = minmod(ut0 - ut_m, ut_p - ut0);

    CellSides s{};
    s.h_lo = h0 - half * dh;
    s.h_hi = h0 + half * dh;
    s.z_lo = -half * dz;
    s.z_hi = half * dz;
    // |dh| <= h0, so the two weights lie between 1/4 and 3/4.
    const amrex::Real w_lo = half * s.h_hi / h0;
    const amrex::Real w_hi = half * s.h_lo / h0;
    s.un_lo = un0 - w_lo * dun;
    s.un_hi = un0 + w_hi * dun;
    s.ut_lo = ut0 - w_lo * dut;
    s.ut_hi = ut0 + w_hi * dut;
    s.d_eta = d_eta;
    return s;
}

} // namespace flume

#endif
