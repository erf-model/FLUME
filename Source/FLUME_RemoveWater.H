#ifndef FLUME_REMOVEWATER_H_
#define FLUME_REMOVEWATER_H_

#include "FLUME_State.H"

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * Remove `loss` (0 <= loss <= h) from the cell and add the depth removed,
 * h - h_new as stored, to the running total `f`, a compensated sum (rounding
 * carried in `carry`, as the depth in add_rain_depth). So the depth plus the
 * total of the cell stays the water it got: adding `loss`, which h - loss
 * rounds, or rounding the sum would let them drift apart with the steps. The
 * discharge keeps the velocity, and a cell left at or below `dry_depth`
 * holds none. Used by the sinks: infiltration (FLUME_Infiltration.H) and the
 * drain inlets (FLUME_Drains.H).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
remove_water (amrex::Array4<amrex::Real> const& a, amrex::Real& f, amrex::Real& carry, int i,
              int j, int k, amrex::Real h, amrex::Real loss, amrex::Real dry_depth) noexcept
{
    const amrex::Real h_new = h - loss;
    const amrex::Real y = (h - h_new) + carry;
    const amrex::Real f_new = f + y;
    carry = y - (f_new - f);
    f = f_new;
    a(i, j, k, State::h) = h_new;
    if (h_new > dry_depth) {
        const amrex::Real r = h_new / h;
        a(i, j, k, State::hu) *= r;
        a(i, j, k, State::hv) *= r;
    } else {
        a(i, j, k, State::hu) = amrex::Real(0.0);
        a(i, j, k, State::hv) = amrex::Real(0.0);
    }
}

} // namespace flume

#endif
