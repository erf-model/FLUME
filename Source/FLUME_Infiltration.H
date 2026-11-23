#ifndef FLUME_INFILTRATION_H_
#define FLUME_INFILTRATION_H_

#include "FLUME_Inputs.H"
#include "FLUME_Rain.H"
#include "FLUME_Source.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
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
 * Infiltration: water that leaves the surface into the soil
 * (flume.infiltration.type, Docs/sphinx_doc/theory/SourceTerms.rst).
 *
 * Each step a loss model gives a depth for each cell; the cell loses that
 * depth, and never more than the water it holds, so the depth stays at 0 or
 * above. The water that leaves takes its share of the discharge with it:
 * the velocity does not change, and a cell left at or below the dry depth
 * holds no discharge. The depth each cell has lost since the start of the run
 * is kept and written to the plotfiles as `infiltrated_depth`, so that the
 * water of the run balances: the water on the surface plus the infiltrated
 * water is the initial water plus the rain, less what leaves through the
 * sides.
 *
 * The soil parameters are per-cell fields, filled from one value or from an
 * ESRI ASCII grid (each cell takes the raster cell that holds its centre).
 * Building cells (flume.buildings.file) infiltrate nothing. The source acts
 * after the rain and before the friction (make_sources), so it sees the
 * rain of the step.
 */
class Infiltration : public Source
{
public:
    [[nodiscard]] std::vector<PlotField> plot_fields () const override
    {
        return {{"infiltrated_depth", &m_infiltrated}};
    }

    //! Depth each cell has lost to the soil since the start, m: one
    //! component, no ghost cells.
    [[nodiscard]] const amrex::MultiFab& infiltrated () const { return m_infiltrated; }

    [[nodiscard]] std::string balance_name () const override { return "infiltrated"; }
    [[nodiscard]] bool balance_removes () const override { return true; }

    //! Minus the infiltrated volume of this rank's cells, m3: the cell area
    //! times the sum of the infiltrated depth and the rounding carried in it.
    [[nodiscard]] double balance_volume () const override;

protected:
    Infiltration (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  const amrex::Geometry& geom, amrex::Real dry_depth);

    amrex::MultiFab m_infiltrated;
    amrex::MultiFab m_infiltrated_carry;  //!< rounding carried in the sum m_infiltrated
    amrex::Real m_dry_depth = amrex::Real(0.0);
    double m_cell_area = 0.0;  //!< m2
};

/**
 * Green-Ampt infiltration for one storm, after Rossman and Huber (2016,
 * section 4.4): EPA SWMM Reference Manual, Volume I, EPA/600/R-15/162A.
 *
 * With K the saturated hydraulic conductivity, psi the suction head at the
 * wetting front, theta the moisture deficit, M = psi theta, F the cumulative
 * infiltration of the cell and W the water the cell holds at the end of the
 * step (rain of the step and ponded water), the step of `dt` is (i = W / dt):
 *
 *   - i <= K: the cell loses all of W;
 *   - otherwise the surface is saturated once F >= F_s = K M / (i - K)
 *     (eq. 4-29). Below F_s the cell loses W if F + W < F_s; else F rises to
 *     F_s at the rate i, and from F_s on follows the integrated equation for
 *     the rest of the step (eq. 4-32):
 *
 *       K tau = (F_2 - F_1) - M ln((F_2 + M) / (F_1 + M)) ,
 *
 *     never losing more than W.
 *
 * The integrated equation is exact for any step, so the ponded infiltration
 * does not depend on the step. The recovery of the soil between storms
 * (section 4.4.2) is not modelled: a run is one storm.
 */
class GreenAmpt : public Infiltration
{
public:
    //! Fill K (m/s), M (m) from flume.infiltration.green_ampt.*; K = 0 in
    //! the cells where `buildings` (may be null) is above 0. Aborts, naming
    //! the key and the file, when a cell takes a value out of range.
    GreenAmpt (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
               const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings);

    [[nodiscard]] std::string name () const override { return "green_ampt_infiltration"; }
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

private:
    amrex::MultiFab m_conductivity;  //!< K, m/s
    amrex::MultiFab m_storage;       //!< M = psi theta, m
};

/**
 * Runoff curve number as a loss on the rain (NRCS National Engineering
 * Handbook Part 630, chapter 10; Woodward et al., 2003).
 *
 * With P the rain fallen on the cell since the start, S the potential
 * retention, S = 25.4 (1000 / CN - 10) mm, and Ia = lambda S the initial
 * abstraction (lambda = flume.infiltration.cn.ia_ratio), the runoff is
 * Q = (P - Ia)^2 / (P - Ia + S) for P > Ia and 0 below, and the cumulative
 * loss is L(P) = P - Q. Each step the cell loses L(P_end) - L(P_start),
 * never more than the water it holds. Water that runs onto a cell from its
 * neighbours is not infiltrated. CN = 100 (S = 0) is an impervious cell.
 */
class CurveNumber : public Infiltration
{
public:
    //! Fill S (m) from flume.infiltration.cn.*; S = 0 in the cells where
    //! `buildings` (may be null) is above 0. `rain` gives the rain of each
    //! step and must act before this source and outlive it.
    CurveNumber (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                 const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings,
                 const Rain& rain);

    [[nodiscard]] std::string name () const override { return "curve_number_infiltration"; }
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

private:
    const Rain& m_rain;
    amrex::MultiFab m_retention;   //!< S, m
    amrex::MultiFab m_rain_total;  //!< P, m
    amrex::MultiFab m_rain_carry;  //!< rounding carried in the sum P, m
    amrex::Real m_ia_ratio = amrex::Real(0.2);
};

/**
 * The depth x = F_2 - F_1 >= 0, m, that Green-Ampt infiltration under
 * ponding adds in the time `tau` (s) from the cumulative infiltration `f1`
 * (m), with K = `k` (m/s) and M = `m` (m): the root of
 *
 *   x - M ln(1 + x / (F_1 + M)) = K tau .
 *
 * Newton's method on the convex, increasing left side, from the bound
 * K tau + sqrt(2 M K tau), so every iterate stays above the root; the
 * left side is written as F_1 u + M (u - ln(1 + u)), u = x / (F_1 + M), with
 * a series for small u, without cancellation. In double precision.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
green_ampt_ponded_step (double f1, double k, double m, double tau) noexcept
{
    const double kt = k * tau;
    if (!(kt > 0.0)) { return 0.0; }
    if (!(m > 0.0)) { return kt; }  // no suction: the rate is K
    const double a = f1 + m;
    // u - ln(1 + u), and its derivative u / (1 + u).
    auto phi = [] (double u) {
        if (u < 0.05) {
            // u^2/2 - u^3/3 + ... - u^9/9: error below 4e-12 u^2 at u = 0.05.
            double term = -u;
            double sum = 0.0;
            for (int n = 2; n <= 9; ++n) {
                term *= -u;
                sum += term / static_cast<double>(n);
            }
            return sum;
        }
        return u - std::log1p(u);
    };
    double u = (kt + std::sqrt(2.0 * m * kt)) / a;
    for (int it = 0; it < 100; ++it) {
        const double g = f1 * u + m * phi(u) - kt;
        const double dg = f1 + m * u / (1.0 + u);
        if (!(g > 0.0) || !(dg > 0.0)) { break; }
        const double u_new = u - g / dg;
        // From above the iterates fall to the root; stop when they no longer do.
        if (!(u_new < u)) { break; }
        u = (u_new > 0.0) ? u_new : 0.0;
    }
    return u * a;
}

/**
 * The cumulative curve-number loss L(P) = P - Q(P), m, for the rain `p` (m),
 * the retention `s` (m) and the initial abstraction `ia` = lambda S (m),
 * written as Ia + (P - Ia) S / (P - Ia + S) above Ia, without cancellation.
 * A curve number so small that S overflows to infinity makes Ia infinite
 * (lambda > 0) or NaN (lambda = 0, 0 x inf); either way the first test
 * returns P: all the rain is lost.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
curve_number_loss (amrex::Real p, amrex::Real s, amrex::Real ia) noexcept
{
    if (!(p > ia)) { return p; }
    const amrex::Real e = p - ia;
    const amrex::Real d = e + s;
    return (d > amrex::Real(0.0)) ? ia + e * s / d : ia;
}

//! The step of GreenAmpt::apply on the valid cells of `state`, with K
//! (`conductivity`, m/s) and M (`storage`, m) per cell; adds the depth
//! removed to `infiltrated` (rounding carried in `infiltrated_carry`). The
//! ghost cells of `state` are not read.
void apply_green_ampt (amrex::MultiFab& state, const amrex::MultiFab& conductivity,
                       const amrex::MultiFab& storage, amrex::MultiFab& infiltrated,
                       amrex::MultiFab& infiltrated_carry, amrex::Real dry_depth,
                       amrex::Real dt);

//! The step of CurveNumber::apply on the valid cells of `state`: adds the
//! step's rain `rain` (m) to `rain_total` (rounding carried in `rain_carry`)
//! and the depth removed to `infiltrated` (rounding carried in
//! `infiltrated_carry`).
void apply_curve_number (amrex::MultiFab& state, const amrex::MultiFab& rain,
                         const amrex::MultiFab& retention, amrex::MultiFab& rain_total,
                         amrex::MultiFab& rain_carry, amrex::MultiFab& infiltrated,
                         amrex::MultiFab& infiltrated_carry, amrex::Real ia_ratio,
                         amrex::Real dry_depth);

} // namespace flume

#endif
