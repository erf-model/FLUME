#include "FLUME_TimeStep.H"
#include "FLUME_State.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeList.H>

#include <cmath>

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

} // namespace flume
