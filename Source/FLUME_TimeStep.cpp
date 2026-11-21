#include "FLUME_TimeStep.H"
#include "FLUME_State.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeList.H>

#include <AMReX.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace flume {

WaveRate
max_wave_rate (const amrex::MultiFab& state, const amrex::Geometry& geom, amrex::Real gravity,
               amrex::Real dry_depth)
{
    const amrex::Real dxi = geom.InvCellSize(0);
    const amrex::Real dyi = geom.InvCellSize(1);
    const auto& s = state.const_arrays();
    const auto r = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpMax, amrex::ReduceOpSum>{},
                                    amrex::TypeList<amrex::Real, amrex::Long>{}, state,
                                    amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Real, amrex::Long>
        {
            const amrex::Real h = s[box_no](i, j, k, State::h);
            // A NaN or infinite depth is counted here, before any comparison:
            // comparing a NaN raises the invalid-operation exception, which
            // the tests trap. A negative depth is counted too. Neither
            // reaches the square root or the division.
            if (!amrex::Math::isfinite(h) || h < amrex::Real(0.0)) {
                return {amrex::Real(0.0), amrex::Long(1)};
            }
            if (h <= dry_depth) {
                // A dry cell has no wave speed, but its discharges must still
                // be numbers: they are written to the plotfile.
                const bool finite = amrex::Math::isfinite(s[box_no](i, j, k, State::hu))
                                    && amrex::Math::isfinite(s[box_no](i, j, k, State::hv));
                return {amrex::Real(0.0), amrex::Long(finite ? 0 : 1)};
            }
            const amrex::Real c = std::sqrt(gravity * h);
            const amrex::Real u = s[box_no](i, j, k, State::hu) / h;
            const amrex::Real v = s[box_no](i, j, k, State::hv) / h;
            const amrex::Real rate = (std::abs(u) + c) * dxi + (std::abs(v) + c) * dyi;
            if (!amrex::Math::isfinite(rate)) {
                return {amrex::Real(0.0), amrex::Long(1)};
            }
            return {rate, amrex::Long(0)};
        });
    WaveRate w;
    w.rate = amrex::get<0>(r);
    w.bad_cells = amrex::get<1>(r);
    amrex::ParallelDescriptor::ReduceRealMax(w.rate);
    amrex::ParallelDescriptor::ReduceLongSum(w.bad_cells);
    return w;
}

double
rain_limited_step (double rate, double r_max, double gravity, double dry_depth, double dxi,
                   double dyi, double cfl)
{
    const double inv = dxi + dyi;
    const double a = std::sqrt(gravity * r_max) * inv;
    const double b = std::max(rate, std::sqrt(gravity * dry_depth) * inv);
    // f(y) = a y^3 + b y^2 - cfl. Start at or above the root: each term alone
    // reaches cfl at its own root, and f >= 0 at both, since the other term
    // is >= 0. The smaller of the two is the closer start.
    double y = std::cbrt(cfl / a);
    if (b > 0.0) { y = std::min(y, std::sqrt(cfl / b)); }
    for (int it = 0; it < 100; ++it) {
        const double f = (a * y + b) * y * y - cfl;
        const double df = (3.0 * a * y + 2.0 * b) * y;
        const double y_new = y - f / df;
        // From above the iterates fall; stop when they stop falling.
        if (!(y_new < y)) { break; }
        y = y_new;
    }
    // The iteration stops at the root to round-off; anything else is a bug.
    const double f = (a * y + b) * y * y - cfl;
    if (!(std::isfinite(y) && y > 0.0 && std::abs(f) <= 1.0e-10 * cfl)) {
        std::ostringstream msg;
        msg << std::setprecision(17) << "The rain-limited time step did not converge: rate " << rate
            << " 1/s, rain rate " << r_max << " m/s, sqrt(dt) " << y << " s^1/2, residual " << f;
        amrex::Abort(msg.str());
    }
    return y * y;
}

} // namespace flume
