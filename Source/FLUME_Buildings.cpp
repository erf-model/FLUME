#include "FLUME_Buildings.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeList.H>

namespace flume {

namespace {

//! Number of valid cells, over all ranks, whose first component is above 0.
amrex::Long
count_above_zero (const amrex::MultiFab& mf)
{
    const auto& a = mf.const_arrays();
    amrex::Long n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                     amrex::TypeList<amrex::Long>{}, mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            return {(a[box_no](i, j, k) > amrex::Real(0.0)) ? amrex::Long(1) : amrex::Long(0)};
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace

amrex::Long
add_building_heights (amrex::MultiFab& bed, const amrex::MultiFab& height)
{
    const auto& zb = bed.arrays();
    const auto& hb = height.const_arrays();
    amrex::ParallelFor(bed,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            zb[box_no](i, j, k) += hb[box_no](i, j, k);
        });
    amrex::Gpu::streamSynchronize();
    return count_above_zero(height);
}

amrex::Long
add_building_footprints (amrex::MultiFab& bed, const amrex::MultiFab& height, amrex::Real fixed)
{
    const auto& zb = bed.arrays();
    const auto& hb = height.const_arrays();
    amrex::ParallelFor(bed,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            if (hb[box_no](i, j, k) > amrex::Real(0.0)) { zb[box_no](i, j, k) += fixed; }
        });
    amrex::Gpu::streamSynchronize();
    return count_above_zero(height);
}

amrex::Long
count_below_zero (const amrex::MultiFab& mf)
{
    const auto& a = mf.const_arrays();
    amrex::Long n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                     amrex::TypeList<amrex::Long>{}, mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            return {(a[box_no](i, j, k) < amrex::Real(0.0)) ? amrex::Long(1) : amrex::Long(0)};
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace flume
