#include "FLUME_Friction.H"
#include "FLUME_State.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>

#include <cmath>

namespace flume {

Friction::Friction (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                    amrex::Real manning, amrex::Real gravity, amrex::Real dry_depth)
    : m_manning(ba, dm, 1, 0), m_gravity(gravity), m_dry_depth(dry_depth)
{
    AMREX_ALWAYS_ASSERT(manning > amrex::Real(0.0));
    m_manning.setVal(manning);
    amrex::Print() << "FLUME: Manning friction, n = " << manning << " s/m^(1/3)\n";
}

void
Friction::apply (amrex::MultiFab& state, const SourceContext& /*ctx*/, amrex::Real dt)
{
    apply_manning_friction(state, m_manning, m_gravity, m_dry_depth, dt);
}

void
apply_manning_friction (amrex::MultiFab& state, const amrex::MultiFab& manning,
                        amrex::Real gravity, amrex::Real dry_depth, amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& mn = manning.const_arrays();
    const amrex::Real one = amrex::Real(1.0);
    // Each cell reads and writes only itself.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            // A dry cell holds no discharge, and a depth that is not finite
            // is left for the check of the state.
            if (!(h > dry_depth)) { return; }
            const amrex::Real qx = a(i, j, k, State::hu);
            const amrex::Real qy = a(i, j, k, State::hv);
            const amrex::Real q = std::sqrt(qx * qx + qy * qy);
            const amrex::Real n = mn[box_no](i, j, k);
            const amrex::Real h73 = h * h * std::cbrt(h);
            // The root in (0, 1] of dt g n^2 |q| / h^(7/3) r^2 + r - 1 = 0.
            // Written without h^(14/3), which can overflow; a ratio that
            // overflows gives r = 0, its limit.
            amrex::Real r = amrex::Real(0.0);
            if (h73 > amrex::Real(0.0)) {
                const amrex::Real stiff = amrex::Real(4.0) * dt * gravity * n * n * q / h73;
                r = amrex::Real(2.0) / (one + std::sqrt(one + stiff));
            }
            a(i, j, k, State::hu) = r * qx;
            a(i, j, k, State::hv) = r * qy;
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
