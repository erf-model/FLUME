#include "FLUME_FloodMaps.H"
#include "FLUME_State.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Math.H>

#include <cmath>

namespace flume {

namespace {
enum : int { max_depth = 0, max_speed, max_depth_speed, arrival_time, n_maps };
} // namespace

FloodMaps::FloodMaps (const Inputs& in, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm)
    : m_maps(ba, dm, n_maps, 0), m_dry_depth(in.dry_depth), m_arrival_depth(in.arrival_depth)
{
    m_maps.setVal(amrex::Real(0.0), 0, arrival_time);
    m_maps.setVal(amrex::Real(-1.0), arrival_time, 1);
    m_alias.reserve(n_maps);
    for (int n = 0; n < n_maps; ++n) {
        m_alias.emplace_back(m_maps, amrex::make_alias, n, 1);
    }
}

void
FloodMaps::update (const amrex::MultiFab& state, double time)
{
    update_flood_maps(m_maps, state, m_dry_depth, m_arrival_depth,
                      static_cast<amrex::Real>(time));
}

std::vector<PlotField>
FloodMaps::plot_fields () const
{
    return {{"max_depth", &m_alias[max_depth]},
            {"max_speed", &m_alias[max_speed]},
            {"max_depth_speed", &m_alias[max_depth_speed]},
            {"arrival_time", &m_alias[arrival_time]}};
}

void
update_flood_maps (amrex::MultiFab& maps, const amrex::MultiFab& state, amrex::Real dry_depth,
                   amrex::Real arrival_depth, amrex::Real time)
{
    const auto& m = maps.arrays();
    const auto& s = state.const_arrays();
    // Each cell reads and writes only itself.
    amrex::ParallelFor(maps,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const auto& b = m[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            const amrex::Real zero = amrex::Real(0.0);
            amrex::Real q = zero;
            amrex::Real speed = zero;
            if (h > dry_depth) {
                const amrex::Real hu = a(i, j, k, State::hu);
                const amrex::Real hv = a(i, j, k, State::hv);
                q = std::sqrt(hu * hu + hv * hv);
                speed = q / h;
            }
            b(i, j, k, max_depth) = amrex::max(b(i, j, k, max_depth), h);
            b(i, j, k, max_speed) = amrex::max(b(i, j, k, max_speed), speed);
            b(i, j, k, max_depth_speed) = amrex::max(b(i, j, k, max_depth_speed), q);
            if (b(i, j, k, arrival_time) < zero && h > arrival_depth) {
                b(i, j, k, arrival_time) = time;
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
