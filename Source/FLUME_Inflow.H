#ifndef FLUME_INFLOW_H_
#define FLUME_INFLOW_H_

#include "FLUME_Inputs.H"
#include "FLUME_TimeTable.H"

#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>
#include <vector>

namespace flume {

/**
 * The depth of the water at the boundary face of an inflow segment, m, for
 * the discharge q >= 0 per unit width that enters (m2/s), given the cell
 * inside: depth h (m) and velocity un (m/s) along the normal that points
 * into the domain (Docs/sphinx_doc/theory/Boundaries.rst).
 *
 * The characteristic of speed un - c leaves the domain when the flow is
 * subcritical and carries the invariant R = un - 2 c, c = sqrt(g h), from
 * inside. The face state (h_b, u_b = q / h_b) keeps it:
 *
 *   q / h_b - 2 sqrt(g h_b) = R .
 *
 * The left side falls from +infinity to -infinity as h_b grows, so there is
 * one root. It is subcritical (u_b < c_b) exactly when it lies above the
 * critical depth h_c = (q / sqrt(g))^(2/3), that is when R < -sqrt(g h_c).
 * Otherwise both characteristics enter the domain, the state at the face
 * comes from outside alone, and it is taken at the critical depth: the
 * result is max(root, h_c). With s = sqrt(h_b) the root is found by Newton's
 * method from sqrt(h_c): the function is convex and decreasing in s, so the
 * iterates rise to the root; the iteration stops when one no longer rises.
 *
 * With q = 0 nothing enters and the face keeps R with u_b = 0:
 * c_b = max(-R / 2, 0), h_b = c_b^2 / g; for still water (un = 0) that is
 * h itself, which the copy returns exactly. A cell at or below the dry depth
 * is passed as h = 0, un = 0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
inflow_depth (amrex::Real q, amrex::Real un, amrex::Real h, amrex::Real g) noexcept
{
    const amrex::Real c = std::sqrt(g * h);
    const amrex::Real r = un - amrex::Real(2.0) * c;
    if (!(q > amrex::Real(0.0))) {
        if (un == amrex::Real(0.0)) { return h; }
        const amrex::Real cb = amrex::max(amrex::Real(-0.5) * r, amrex::Real(0.0));
        return cb * cb / g;
    }
    const amrex::Real sg = std::sqrt(g);
    const amrex::Real a = std::cbrt(q / sg);
    const amrex::Real hc = a * a;
    if (!(r < -sg * a)) { return hc; }  // sqrt(g h_c) = sqrt(g) a
    amrex::Real s = a;
    for (int it = 0; it < 100; ++it) {
        const amrex::Real s2 = s * s;
        const amrex::Real f = q / s2 - amrex::Real(2.0) * sg * s - r;
        const amrex::Real fp = amrex::Real(-2.0) * q / (s2 * s) - amrex::Real(2.0) * sg;
        const amrex::Real next = s - f / fp;
        if (!(next > s)) { break; }
        s = next;
    }
    return amrex::max(s * s, hc);
}

/**
 * Inflow hydrographs on boundary segments (flume.bc.<name>.type =
 * discharge): the discharge Q(t), m3/s, of each segment, from its file
 * (flume.bc.<name>.file), enters through the boundary faces of the segment,
 * spread evenly over its width W: q = Q / W per unit width.
 *
 * Each face carries the mass flux q exactly, so the water that enters over a
 * step is the integral of the hydrograph over the step (TimeTable). The
 * momentum flux along the normal is q u_b + g h_b^2 / 2, with h_b from the
 * cell inside the face (inflow_depth) and u_b = q / h_b; the flux of the
 * momentum along the side is zero (the water enters along the normal). The
 * ghost cells of the segment are filled as for an outflow segment (a copy of
 * the boundary cell), which the face fluxes then replace.
 *
 * Nothing changes from step to step but the time, so a checkpoint holds
 * nothing of its own; the water balance counts the inflow of each segment
 * (FLUME_WaterBalance.H).
 */
class Inflow
{
public:
    //! The discharge segments of `in` on the grid of `geom`: reads each file
    //! (collective) and prints the hydrograph.
    Inflow (const Inputs& in, const amrex::Geometry& geom);

    //! Whether there is no discharge segment.
    [[nodiscard]] bool empty () const { return m_segments.empty(); }

    /**
     * Replace the fluxes of `flux` on the boundary faces of each discharge
     * segment: the mean discharge of the step from `t0` to `t1` (s, t1 > t0)
     * over the step `dt` (s, the step of the update), (I(t1) - I(t0)) / dt,
     * per unit width, and its momentum flux from `state` (valid cells read
     * only). Both stages of the two-stage step take the same discharge, so
     * the step adds the integral of the hydrograph.
     */
    void set_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                     double t0, double t1, double dt) const;

    /**
     * The largest wave rate, 1/s, of the face states of the discharge
     * segments with the largest discharge from `t0` to `t1` (t1 may be
     * +infinity): (u_b + c_b) / dn + c_b / dt, dn the cell size along the
     * normal and dt along the side, as max_wave_rate of FLUME_TimeStep.H
     * counts a cell. The face depth grows with q, so this bounds the rates of
     * the step. Zero without a discharge segment. Collective.
     */
    [[nodiscard]] amrex::Real max_rate (const amrex::MultiFab& state, double t0, double t1) const;

private:
    struct Segment
    {
        std::string name;
        int d = 0;
        int high = 0;
        int lo = 0;
        int hi = -1;
        double width = 0.0;  //!< m
        TimeTable q;         //!< discharge, m3/s
    };

    amrex::Geometry m_geom;
    amrex::Real m_gravity = amrex::Real(9.81);
    amrex::Real m_dry_depth = amrex::Real(0.0);
    std::vector<Segment> m_segments;
};

} // namespace flume

#endif
