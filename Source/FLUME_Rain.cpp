#include "FLUME_Rain.H"
#include "FLUME_State.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>

namespace flume {

Rain::Rain (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
            const Forcing& forcing)
    : m_forcing(forcing), m_depth(ba, dm, 1, 0), m_carry(ba, dm, 1, 0)
{
    m_carry.setVal(amrex::Real(0.0));
}

void
Rain::apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real /*dt*/)
{
    // The ends of the step in double precision, as the run's time is: in
    // single precision the step dt differs from time_end - time by rounding.
    m_forcing.rain_depth(m_depth, ctx.time, ctx.time_end);
    add_rain_depth(state, m_depth, m_carry);
}

void
add_rain_depth (amrex::MultiFab& state, const amrex::MultiFab& depth, amrex::MultiFab& carry)
{
    const auto& s = state.arrays();
    const auto& d = depth.const_arrays();
    const auto& c = carry.arrays();
    // Each cell reads and writes only itself. Compensated (Kahan) sum: the
    // part of the rain that the addition rounds away is kept in `carry` and
    // added with the next step's rain, so that it is not lost when the rain
    // of a step is below the rounding of the depth (single precision on
    // deep water).
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real h = s[box_no](i, j, k, State::h);
            const amrex::Real y = d[box_no](i, j, k) + c[box_no](i, j, k);
            const amrex::Real h_new = h + y;
            c[box_no](i, j, k) = y - (h_new - h);
            s[box_no](i, j, k, State::h) = h_new;
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
