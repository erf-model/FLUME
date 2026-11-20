#ifndef FLUME_RIEMANN_H_
#define FLUME_RIEMANN_H_

#include <AMReX_Algorithm.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace flume {

//! Flux across a face: mass, momentum normal to the face, momentum along it.
struct FaceFlux
{
    amrex::Real mass;    //!< m2/s
    amrex::Real normal;  //!< m3/s2
    amrex::Real tangent; //!< m3/s2
};

//! The depths the two cells of a face present to it.
struct FaceDepths
{
    amrex::Real lo; //!< of the cell on the low side
    amrex::Real hi; //!< of the cell on the high side
};

/**
 * Depths the two cells of a face present to it (hydrostatic reconstruction,
 * Audusse et al. 2004). `hL`, `zL` are the depth and the bed elevation of the
 * cell on the low side, `hR`, `zR` those of the cell on the high side. The
 * face bed is the higher of the two beds, and the depth of a cell is what its
 * surface leaves above it: h - max(0, z_other - z). A depth at or below
 * `dry_depth` (>= 0) counts as none: it is returned as zero.
 *
 * Under a level surface the two depths are equal up to rounding (they are
 * computed from different numbers). Where one is above `dry_depth`, the other
 * is not and they differ by rounding only (4 machine epsilons of the larger
 * cell depth), both are returned as zero: the cut is the same on the two
 * sides, and water at rest stays at rest for any `dry_depth`. The flux and
 * the bed-slope term call this with the same arguments, so they see the same
 * depths.
 *
 * Only the difference of the two bed elevations enters, so the result does
 * not lose digits when the elevations are large.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE FaceDepths
face_depths (amrex::Real hL, amrex::Real zL, amrex::Real hR, amrex::Real zR,
             amrex::Real dry_depth) noexcept
{
    const amrex::Real zero = amrex::Real(0.0);
    const amrex::Real dz = zR - zL;
    const amrex::Real hsL = hL - amrex::max(dz, zero);
    const amrex::Real hsR = hR - amrex::max(-dz, zero);
    bool wetL = hsL > dry_depth;
    bool wetR = hsR > dry_depth;
    if (wetL != wetR) {
        const amrex::Real tol = amrex::Real(4.0) * std::numeric_limits<amrex::Real>::epsilon()
                                * amrex::max(hL, hR);
        if (amrex::Math::abs(hsL - hsR) <= tol) {
            wetL = false;
            wetR = false;
        }
    }
    return {wetL ? hsL : zero, wetR ? hsR : zero};
}

/**
 * HLLC flux of the shallow-water equations across a face, from the state on
 * its low side (L) and on its high side (R): depth h (>= 0), velocity normal
 * to the face un and along it ut. `g` is the acceleration of gravity.
 *
 * The wave speeds are estimated from the middle state that two rarefactions
 * would leave (Docs/sphinx_doc/Model.rst), with the speeds of a front on a
 * dry bed when one side has no water. Mass and
 * normal momentum take the HLL flux. The tangential momentum is carried by
 * the mass flux from its upwind side, which restores the contact wave.
 * A face that is dry on both sides has no flux.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE FaceFlux
hllc_flux (amrex::Real hL, amrex::Real unL, amrex::Real utL, amrex::Real hR, amrex::Real unR,
           amrex::Real utR, amrex::Real g) noexcept
{
    const amrex::Real zero = amrex::Real(0.0);
    const amrex::Real half = amrex::Real(0.5);
    if (!(hL > zero) && !(hR > zero)) { return {zero, zero, zero}; }

    const amrex::Real cL = std::sqrt(g * hL);
    const amrex::Real cR = std::sqrt(g * hR);
    amrex::Real sL = zero;
    amrex::Real sR = zero;
    if (!(hR > zero)) {
        // Dry on the high side: the front runs at uL + 2 cL.
        sL = unL - cL;
        sR = unL + amrex::Real(2.0) * cL;
    } else if (!(hL > zero)) {
        sL = unR - amrex::Real(2.0) * cR;
        sR = unR + cR;
    } else {
        const amrex::Real us = half * (unL + unR) + cL - cR;
        // sqrt(g h*) of the two-rarefaction star state; not below zero.
        const amrex::Real cs = amrex::max(half * (cL + cR) + amrex::Real(0.25) * (unL - unR), zero);
        sL = amrex::min(unL - cL, us - cs);
        sR = amrex::max(unR + cR, us + cs);
    }

    const amrex::Real fL_mass = hL * unL;
    const amrex::Real fR_mass = hR * unR;
    const amrex::Real fL_norm = fL_mass * unL + half * g * hL * hL;
    const amrex::Real fR_norm = fR_mass * unR + half * g * hR * hR;

    FaceFlux f{zero, zero, zero};
    if (sL >= zero) {
        f.mass = fL_mass;
        f.normal = fL_norm;
    } else if (sR <= zero) {
        f.mass = fR_mass;
        f.normal = fR_norm;
    } else {
        const amrex::Real inv = amrex::Real(1.0) / (sR - sL);
        f.mass = (sR * fL_mass - sL * fR_mass + sL * sR * (hR - hL)) * inv;
        f.normal = (sR * fL_norm - sL * fR_norm + sL * sR * (fR_mass - fL_mass)) * inv;
    }
    // The contact wave moves with the sign of the mass flux.
    f.tangent = f.mass * ((f.mass >= zero) ? utL : utR);
    return f;
}

} // namespace flume

#endif
