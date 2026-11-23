#include "FLUME_Rain.H"
#include "FLUME_State.H"
#include "FLUME_WaterBalance.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>

namespace flume {

Rain::Rain (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
            const amrex::Geometry& geom, const Forcing& forcing)
    : m_forcing(forcing), m_depth(ba, dm, 1, 0), m_carry(ba, dm, 1, 0),
      m_cell_area(cell_area(geom))
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
    compensated_add(m_volume, m_volume_carry, volume(m_cell_area, local_sum(m_depth)));
}

double
Rain::balance_volume () const
{
    return m_volume - volume(m_cell_area, local_sum(m_carry));
}

void
Rain::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_multifab(name() + "_carry", m_carry);
    chk.put_per_rank(name() + "/volume", {m_volume, m_volume_carry});
}

void
Rain::read_checkpoint (const CheckpointReader& chk)
{
    chk.get_multifab(name() + "_carry", m_carry);
    const std::vector<double> v = chk.get_per_rank(name() + "/volume", 2);
    m_volume = v[0];
    m_volume_carry = v[1];
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
            // A negative carry (the rain a step rounded up) waits while the
            // cell holds less: rain never takes water away, and a cell that
            // another source emptied stays at 0 instead of below it.
            if (!(h + y >= amrex::Real(0.0))) {
                c[box_no](i, j, k) = y;
                return;
            }
            const amrex::Real h_new = h + y;
            c[box_no](i, j, k) = y - (h_new - h);
            s[box_no](i, j, k, State::h) = h_new;
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
