#include "FLUME_Infiltration.H"
#include "FLUME_Forcing.H"
#include "FLUME_Raster.H"
#include "FLUME_State.H"
#include "FLUME_WaterBalance.H"

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <limits>
#include <sstream>

namespace flume {

namespace {

//! Number of valid cells, over all ranks, whose value is outside the range
//! from `lo` to `hi` (each end allowed or not), or not finite.
amrex::Long
count_outside (const amrex::MultiFab& mf, amrex::Real lo, bool lo_in, amrex::Real hi, bool hi_in)
{
    const auto& a = mf.const_arrays();
    amrex::Long n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                     amrex::TypeList<amrex::Long>{}, mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            const amrex::Real v = a[box_no](i, j, k);
            const bool ok = (lo_in ? v >= lo : v > lo) && (hi_in ? v <= hi : v < hi);
            return {ok ? amrex::Long(0) : amrex::Long(1)};
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

/**
 * Fill `dst` with the soil parameter `p` (flume.infiltration.<key>): its
 * value, or its file sampled at each cell centre. A file whose cells take a
 * value outside the range aborts, naming the key and the file. Returns a
 * description of what was used, in the units of the input, for the
 * start-up message.
 */
std::string
fill_soil (amrex::MultiFab& dst, const SoilParameter& p, const std::string& key,
           const amrex::Geometry& geom, amrex::Real lo, bool lo_in, amrex::Real hi, bool hi_in,
           const std::string& range)
{
    std::ostringstream what;
    if (p.file.empty()) {
        dst.setVal(p.value);
        what << p.value;
        return what.str();
    }
    const Raster r = read_raster("infiltration." + key + "_file", p.file, geom);
    sample_nearest(dst, 0, r, geom);
    const amrex::Long bad = count_outside(dst, lo, lo_in, hi, hi_in);
    if (bad > 0) {
        amrex::Abort("flume.infiltration." + key + "_file (" + p.file + "): "
                     + std::to_string(bad) + " cell(s) of the domain take a value that is not "
                     + range);
    }
    what << "from " << dst.min(0) << " to " << dst.max(0) << " (" << p.file << ")";
    return what.str();
}

//! Set `field` to `value` in the cells where `buildings` (may be null) is
//! above 0. Returns the number of such cells over all ranks.
amrex::Long
zero_in_buildings (amrex::MultiFab& field, const amrex::MultiFab* buildings, amrex::Real value)
{
    if (buildings == nullptr) { return 0; }
    const auto& f = field.arrays();
    const auto& b = buildings->const_arrays();
    amrex::Long n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                     amrex::TypeList<amrex::Long>{}, field, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Long>
        {
            if (b[box_no](i, j, k) > amrex::Real(0.0)) {
                f[box_no](i, j, k) = value;
                return {amrex::Long(1)};
            }
            return {amrex::Long(0)};
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace

Infiltration::Infiltration (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                            const amrex::Geometry& geom, amrex::Real dry_depth)
    : m_infiltrated(ba, dm, 1, 0), m_infiltrated_carry(ba, dm, 1, 0), m_dry_depth(dry_depth),
      m_cell_area(cell_area(geom))
{
    m_infiltrated.setVal(amrex::Real(0.0));
    m_infiltrated_carry.setVal(amrex::Real(0.0));
}

double
Infiltration::balance_volume () const
{
    // Each cell has lost its total plus the carry of remove_water: with the
    // carry the volume follows the state to round-off of a double, also in
    // a single-precision build.
    return -(volume(m_cell_area, local_sum(m_infiltrated))
             + volume(m_cell_area, local_sum(m_infiltrated_carry)));
}

void
Infiltration::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_multifab(name() + "_depth", m_infiltrated);
    chk.put_multifab(name() + "_depth_carry", m_infiltrated_carry);
}

void
Infiltration::read_checkpoint (const CheckpointReader& chk)
{
    chk.get_multifab(name() + "_depth", m_infiltrated);
    chk.get_multifab(name() + "_depth_carry", m_infiltrated_carry);
}

void
CurveNumber::write_checkpoint (CheckpointWriter& chk) const
{
    Infiltration::write_checkpoint(chk);
    chk.put_multifab(name() + "_rain", m_rain_total);
    chk.put_multifab(name() + "_rain_carry", m_rain_carry);
}

void
CurveNumber::read_checkpoint (const CheckpointReader& chk)
{
    Infiltration::read_checkpoint(chk);
    chk.get_multifab(name() + "_rain", m_rain_total);
    chk.get_multifab(name() + "_rain_carry", m_rain_carry);
}

GreenAmpt::GreenAmpt (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings)
    : Infiltration(ba, dm, geom, in.dry_depth), m_conductivity(ba, dm, 1, 0), m_storage(ba, dm, 1, 0)
{
    const amrex::Real zero = amrex::Real(0.0);
    const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    const std::string k = fill_soil(m_conductivity, in.ga_conductivity,
                                    "green_ampt.conductivity", geom, zero, true, inf, false,
                                    ">= 0 (mm/h)");
    amrex::MultiFab deficit(ba, dm, 1, 0);
    const std::string psi = fill_soil(m_storage, in.ga_suction, "green_ampt.suction", geom, zero,
                                      true, inf, false, ">= 0 (mm)");
    const std::string theta = fill_soil(deficit, in.ga_deficit, "green_ampt.deficit", geom, zero,
                                        true, amrex::Real(1.0), true, ">= 0 and <= 1");
    // mm/h to m/s; M = psi theta, mm to m.
    m_conductivity.mult(static_cast<amrex::Real>(mm_per_h));
    amrex::MultiFab::Multiply(m_storage, deficit, 0, 0, 1, 0);
    m_storage.mult(amrex::Real(1.0e-3));
    const amrex::Long n_building = zero_in_buildings(m_conductivity, buildings, zero);
    amrex::Print() << "FLUME: Green-Ampt infiltration: conductivity " << k << " mm/h, suction "
                   << psi << " mm, deficit " << theta;
    if (n_building > 0) { amrex::Print() << "; " << n_building << " building cell(s) impervious"; }
    amrex::Print() << "\n";
}

void
GreenAmpt::apply (amrex::MultiFab& state, const SourceContext& /*ctx*/, amrex::Real dt)
{
    apply_green_ampt(state, m_conductivity, m_storage, m_infiltrated, m_infiltrated_carry,
                     m_dry_depth, dt);
}

CurveNumber::CurveNumber (const Inputs& in, const amrex::Geometry& geom,
                          const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                          const amrex::MultiFab* buildings, const Rain& rain)
    : Infiltration(ba, dm, geom, in.dry_depth), m_rain(rain), m_retention(ba, dm, 1, 0),
      m_rain_total(ba, dm, 1, 0), m_rain_carry(ba, dm, 1, 0), m_ia_ratio(in.cn_ia_ratio)
{
    const std::string cn = fill_soil(m_retention, in.cn_number, "cn.number", geom,
                                     amrex::Real(0.0), false, amrex::Real(100.0), true,
                                     "> 0 and <= 100");
    // S = (1000 / CN - 10) inches, in m.
    const auto& s = m_retention.arrays();
    amrex::ParallelFor(m_retention,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const amrex::Real c = s[box_no](i, j, k);
            const amrex::Real inches = amrex::Real(1000.0) / c - amrex::Real(10.0);
            s[box_no](i, j, k) = amrex::max(inches, amrex::Real(0.0)) * amrex::Real(0.0254);
        });
    amrex::Gpu::streamSynchronize();
    m_rain_total.setVal(amrex::Real(0.0));
    m_rain_carry.setVal(amrex::Real(0.0));
    const amrex::Long n_building = zero_in_buildings(m_retention, buildings, amrex::Real(0.0));
    amrex::Print() << "FLUME: curve-number infiltration: CN " << cn << ", Ia / S = "
                   << m_ia_ratio;
    if (n_building > 0) { amrex::Print() << "; " << n_building << " building cell(s) impervious"; }
    amrex::Print() << "\n";
}

void
CurveNumber::apply (amrex::MultiFab& state, const SourceContext& /*ctx*/, amrex::Real /*dt*/)
{
    apply_curve_number(state, m_rain.step_depth(), m_retention, m_rain_total, m_rain_carry,
                       m_infiltrated, m_infiltrated_carry, m_ia_ratio, m_dry_depth);
}

namespace {

//! Remove `loss` (0 <= loss <= h) from the cell and add the depth removed,
//! h - h_new as stored, to the infiltrated depth `f`, a compensated sum
//! (rounding carried in `carry`, as the depth in add_rain_depth). So the
//! depth plus the infiltrated depth of the cell stays the water it got:
//! adding `loss`, which h - loss rounds, or rounding the sum would let them
//! drift apart with the steps. The discharge keeps the velocity, and a cell
//! left at or below `dry_depth` holds none.
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

} // namespace

void
apply_green_ampt (amrex::MultiFab& state, const amrex::MultiFab& conductivity,
                  const amrex::MultiFab& storage, amrex::MultiFab& infiltrated,
                  amrex::MultiFab& infiltrated_carry, amrex::Real dry_depth, amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& kk = conductivity.const_arrays();
    const auto& mm = storage.const_arrays();
    const auto& ff = infiltrated.arrays();
    const auto& fc = infiltrated_carry.arrays();
    // Each cell reads and writes only itself. In double precision.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            const double w = static_cast<double>(h);
            const double kc = static_cast<double>(kk[box_no](i, j, k));
            // No water, a depth that is not finite (left for the check of
            // the state), or an impervious cell.
            if (!(w > 0.0) || !(kc > 0.0)) { return; }
            const double m = static_cast<double>(mm[box_no](i, j, k));
            const double f = static_cast<double>(ff[box_no](i, j, k));
            const double tau = static_cast<double>(dt);
            const double rate = w / tau;  // the water available, as a rate
            double df = w;
            if (rate > kc) {
                const double fs = kc * m / (rate - kc);  // eq. 4-29
                if (f >= fs) {
                    df = green_ampt_ponded_step(f, kc, m, tau);
                } else if (f + w >= fs) {
                    // Saturated after (F_s - F) / i of the step.
                    const double t_sat = (fs - f) / rate;
                    df = (fs - f) + green_ampt_ponded_step(fs, kc, m, tau - t_sat);
                }
            }
            // Never more than the water the cell holds. When i <= K the cell
            // loses W; when ponded the rate K (1 + M / F) is at most i (F >= F_s),
            // so the loss is at most i dt = W in exact arithmetic: this guards
            // against rounding only.
            const amrex::Real loss = (df < w) ? static_cast<amrex::Real>(df) : h;
            remove_water(a, ff[box_no](i, j, k), fc[box_no](i, j, k), i, j, k, h, loss, dry_depth);
        });
    amrex::Gpu::streamSynchronize();
}

void
apply_curve_number (amrex::MultiFab& state, const amrex::MultiFab& rain,
                    const amrex::MultiFab& retention, amrex::MultiFab& rain_total,
                    amrex::MultiFab& rain_carry, amrex::MultiFab& infiltrated,
                    amrex::MultiFab& infiltrated_carry, amrex::Real ia_ratio,
                    amrex::Real dry_depth)
{
    const auto& s = state.arrays();
    const auto& r = rain.const_arrays();
    const auto& ss = retention.const_arrays();
    const auto& pp = rain_total.arrays();
    const auto& pc = rain_carry.arrays();
    const auto& ff = infiltrated.arrays();
    const auto& fc = infiltrated_carry.arrays();
    // Each cell reads and writes only itself.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            // P is a compensated sum, as the depth in add_rain_depth: the
            // rain of a step can be below the rounding of P (single
            // precision, many steps).
            const amrex::Real p0 = pp[box_no](i, j, k);
            const amrex::Real y = r[box_no](i, j, k) + pc[box_no](i, j, k);
            const amrex::Real p1 = p0 + y;
            pc[box_no](i, j, k) = y - (p1 - p0);
            pp[box_no](i, j, k) = p1;
            const amrex::Real sr = ss[box_no](i, j, k);
            const amrex::Real ia = ia_ratio * sr;
            const amrex::Real want = curve_number_loss(p1, sr, ia) - curve_number_loss(p0, sr, ia);
            const auto& a = s[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            if (!(want > amrex::Real(0.0)) || !(h > amrex::Real(0.0))) { return; }
            // The loss of a step is at most its rain (dL/dP <= 1), which the
            // cell holds after the rain source; the cap guards against rounding.
            const amrex::Real loss = (want < h) ? want : h;
            remove_water(a, ff[box_no](i, j, k), fc[box_no](i, j, k), i, j, k, h, loss, dry_depth);
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
