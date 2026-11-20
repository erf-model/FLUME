#include "FLUME_Advance.H"
#include "FLUME_Reconstruct.H"
#include "FLUME_Riemann.H"
#include "FLUME_State.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFParallelFor.H>

#include <limits>

namespace flume {

void
compute_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                const amrex::MultiFab& bed, amrex::Real gravity, amrex::Real dry_depth)
{
    const auto& s = state.const_arrays();
    const auto& zb = bed.const_arrays();
    for (int d = 0; d < 2; ++d) {
        // Momentum normal to the face and along it.
        const int i_norm = State::hu + d;
        const int i_tang = State::hu + 1 - d;
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        const auto& f = flux[d].arrays();
        amrex::ParallelFor(flux[d],
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                // Face (i, j) across direction d lies between cell (i, j) - e
                // and cell (i, j).
                const amrex::IntVect R(i, j, k);
                const amrex::IntVect L = R - e;
                const auto& a = s[box_no];
                const amrex::Real hL = a(L, State::h);
                const amrex::Real hR = a(R, State::h);
                const amrex::Real zL = zb[box_no](L);
                const amrex::Real zR = zb[box_no](R);
                // The velocity of a dry cell is not defined: zero.
                const amrex::Real zero = amrex::Real(0.0);
                const amrex::Real unL = (hL > dry_depth) ? a(L, i_norm) / hL : zero;
                const amrex::Real utL = (hL > dry_depth) ? a(L, i_tang) / hL : zero;
                const amrex::Real unR = (hR > dry_depth) ? a(R, i_norm) / hR : zero;
                const amrex::Real utR = (hR > dry_depth) ? a(R, i_tang) / hR : zero;
                const FaceDepths hs = face_depths(hL, zL, hR, zR, dry_depth);
                const FaceFlux F = hllc_flux(hs.lo, unL, utL, hs.hi, unR, utR, gravity);
                f[box_no](i, j, k, State::h) = F.mass;
                f[box_no](i, j, k, i_norm) = F.normal;
                f[box_no](i, j, k, i_tang) = F.tangent;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
compute_fluxes_second_order (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                             const amrex::MultiFab& bed, amrex::Real gravity,
                             amrex::Real dry_depth)
{
    const auto& s = state.const_arrays();
    const auto& zb = bed.const_arrays();
    for (int d = 0; d < 2; ++d) {
        const int i_norm = State::hu + d;
        const int i_tang = State::hu + 1 - d;
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        const auto& f = flux[d].arrays();
        amrex::ParallelFor(flux[d],
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::IntVect R(i, j, k);
                const amrex::IntVect L = R - e;
                const auto& a = s[box_no];
                const auto& z = zb[box_no];
                // The high side of the cell below the face and the low side
                // of the cell above it.
                const CellSides lo = reconstruct(a, z, L, e, i_norm, i_tang, dry_depth);
                const CellSides hi = reconstruct(a, z, R, e, i_norm, i_tang, dry_depth);
                // Step of the bed across the face, from differences only.
                const amrex::Real dz = (z(R) - z(L)) + (hi.z_lo - lo.z_hi);
                const FaceDepths hs = face_depths(lo.h_hi, amrex::Real(0.0), hi.h_lo, dz, dry_depth);
                const FaceFlux F = hllc_flux(hs.lo, lo.un_hi, lo.ut_hi, hs.hi, hi.un_lo, hi.ut_lo,
                                             gravity);
                f[box_no](i, j, k, State::h) = F.mass;
                f[box_no](i, j, k, i_norm) = F.normal;
                f[box_no](i, j, k, i_tang) = F.tangent;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
limit_outflow (amrex::Array<amrex::MultiFab, 2>& flux, amrex::MultiFab& factor,
               const amrex::MultiFab& state, const amrex::Geometry& geom, amrex::Real dt)
{
    const amrex::Real dtdx = dt / geom.CellSize(0);
    const amrex::Real dtdy = dt / geom.CellSize(1);
    // Ghost cells outside the domain keep 1. Through an outflow piece, water
    // leaving is limited by the factor of the cell inside, and water entering
    // takes the ghost cell's 1 (no limit). The mass flux through a wall is
    // zero, or the rounding of one product where the compiler contracts the
    // flux into a fused multiply-add: mass then leaves the cell inside, and
    // the factor of the ghost cell is not the one used.
    factor.setVal(amrex::Real(1.0));
    {
        const auto& s = state.const_arrays();
        const auto& fx = flux[0].const_arrays();
        const auto& fy = flux[1].const_arrays();
        const auto& p = factor.arrays();
        amrex::ParallelFor(state,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::Real zero = amrex::Real(0.0);
                const auto& f = fx[box_no];
                const auto& g = fy[box_no];
                const amrex::Real out =
                    dtdx * (amrex::max(f(i + 1, j, k, State::h), zero)
                            + amrex::max(-f(i, j, k, State::h), zero))
                    + dtdy * (amrex::max(g(i, j + 1, k, State::h), zero)
                              + amrex::max(-g(i, j, k, State::h), zero));
                const amrex::Real h = s[box_no](i, j, k, State::h);
                if (out > h && out > zero) { p[box_no](i, j, k) = amrex::max(h, zero) / out; }
            });
        amrex::Gpu::streamSynchronize();
    }
    factor.FillBoundary(geom.periodicity());
    const auto& p = factor.const_arrays();
    for (int d = 0; d < 2; ++d) {
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        const auto& f = flux[d].arrays();
        amrex::ParallelFor(flux[d],
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::IntVect R(i, j, k);
                const amrex::Real mass = f[box_no](R, State::h);
                amrex::Real scale = amrex::Real(1.0);
                if (mass > amrex::Real(0.0)) {
                    scale = p[box_no](R - e);
                } else if (mass < amrex::Real(0.0)) {
                    scale = p[box_no](R);
                }
                for (int n = 0; n < State::ncomp; ++n) { f[box_no](R, n) *= scale; }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
compute_bed_source (amrex::MultiFab& source, const amrex::MultiFab& state,
                    const amrex::MultiFab& bed, amrex::Real gravity, amrex::Real dry_depth)
{
    const auto& s = state.const_arrays();
    const auto& zb = bed.const_arrays();
    const auto& src = source.arrays();
    const amrex::Real half_g = amrex::Real(0.5) * gravity;
    amrex::ParallelFor(source,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const auto& z = zb[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            const amrex::Real zc = z(i, j, k);
            // The same calls as in compute_fluxes: low cell first.
            const amrex::Real h_xlo =
                face_depths(a(i - 1, j, k, State::h), z(i - 1, j, k), h, zc, dry_depth).hi;
            const amrex::Real h_xhi =
                face_depths(h, zc, a(i + 1, j, k, State::h), z(i + 1, j, k), dry_depth).lo;
            const amrex::Real h_ylo =
                face_depths(a(i, j - 1, k, State::h), z(i, j - 1, k), h, zc, dry_depth).hi;
            const amrex::Real h_yhi =
                face_depths(h, zc, a(i, j + 1, k, State::h), z(i, j + 1, k), dry_depth).lo;
            // h_hi^2 - h_lo^2 as a product: exactly zero for equal depths,
            // also where the compiler contracts a * a - b * b into a fused
            // multiply-add (which would leave the rounding error of b * b).
            src[box_no](i, j, k, 0) = half_g * ((h_xhi - h_xlo) * (h_xhi + h_xlo));
            src[box_no](i, j, k, 1) = half_g * ((h_yhi - h_ylo) * (h_yhi + h_ylo));
        });
    amrex::Gpu::streamSynchronize();
}

void
compute_bed_source_second_order (amrex::MultiFab& source, const amrex::MultiFab& state,
                                 const amrex::MultiFab& bed, amrex::Real gravity,
                                 amrex::Real dry_depth)
{
    const auto& s = state.const_arrays();
    const auto& zb = bed.const_arrays();
    const auto& src = source.arrays();
    const amrex::Real half_g = amrex::Real(0.5) * gravity;
    for (int d = 0; d < 2; ++d) {
        const int i_norm = State::hu + d;
        const int i_tang = State::hu + 1 - d;
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        amrex::ParallelFor(source,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                const amrex::IntVect c(i, j, k);
                const auto& a = s[box_no];
                const auto& z = zb[box_no];
                const CellSides below = reconstruct(a, z, c - e, e, i_norm, i_tang, dry_depth);
                const CellSides here = reconstruct(a, z, c, e, i_norm, i_tang, dry_depth);
                const CellSides above = reconstruct(a, z, c + e, e, i_norm, i_tang, dry_depth);
                // The same calls as in compute_fluxes_second_order.
                const amrex::Real dz_lo = (z(c) - z(c - e)) + (here.z_lo - below.z_hi);
                const amrex::Real dz_hi = (z(c + e) - z(c)) + (above.z_lo - here.z_hi);
                const amrex::Real zero = amrex::Real(0.0);
                const amrex::Real h_lo = face_depths(below.h_hi, zero, here.h_lo, dz_lo, dry_depth).hi;
                const amrex::Real h_hi = face_depths(here.h_hi, zero, above.h_lo, dz_hi, dry_depth).lo;
                // The steps of the bed at the two faces and the slope of the
                // bed inside the cell add up to the first-order term with the
                // depths of the sides, minus g h d_eta
                // (Docs/sphinx_doc/Model.rst).
                src[box_no](c, d) = half_g * ((h_hi - h_lo) * (h_hi + h_lo))
                                    - gravity * a(c, State::h) * here.d_eta;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
update_state (amrex::MultiFab& state, const amrex::Array<amrex::MultiFab, 2>& flux,
              const amrex::MultiFab& source, const amrex::Geometry& geom,
              amrex::Real dry_depth, amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& sb = source.const_arrays();
    const auto& fx = flux[0].const_arrays();
    const auto& fy = flux[1].const_arrays();
    const amrex::Real dtdx = dt / geom.CellSize(0);
    const amrex::Real dtdy = dt / geom.CellSize(1);
    const amrex::Real round_off = amrex::Real(16.0) * std::numeric_limits<amrex::Real>::epsilon();
    // A cell reads its own state, its bed-slope term (computed from the
    // state before the step) and the fluxes on its four faces, and writes
    // its own state: the update is in place and needs no atomics.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const auto& f = fx[box_no];
            const auto& g = fy[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            for (int n = 0; n < State::ncomp; ++n) {
                a(i, j, k, n) -= dtdx * (f(i + 1, j, k, n) - f(i, j, k, n))
                                 + dtdy * (g(i, j + 1, k, n) - g(i, j, k, n));
            }
            a(i, j, k, State::hu) += dtdx * sb[box_no](i, j, k, 0);
            a(i, j, k, State::hv) += dtdy * sb[box_no](i, j, k, 1);
            // The rounding of a cell that limit_outflow emptied. A depth
            // further below zero is left for the check of the state to stop.
            const amrex::Real h_new = a(i, j, k, State::h);
            if (h_new < amrex::Real(0.0) && h_new >= -round_off * h) {
                a(i, j, k, State::h) = amrex::Real(0.0);
            }
            // A dry cell holds no momentum. A depth that is not finite is
            // kept for the check of the state (or stops here, where invalid
            // operations trap).
            if (a(i, j, k, State::h) <= dry_depth) {
                a(i, j, k, State::hu) = amrex::Real(0.0);
                a(i, j, k, State::hv) = amrex::Real(0.0);
            }
        });
    amrex::Gpu::streamSynchronize();
}

void
average_with (amrex::MultiFab& state, const amrex::MultiFab& old, amrex::Real dry_depth)
{
    const auto& s = state.arrays();
    const auto& o = old.const_arrays();
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            for (int n = 0; n < State::ncomp; ++n) {
                a(i, j, k, n) = amrex::Real(0.5) * (o[box_no](i, j, k, n) + a(i, j, k, n));
            }
            if (a(i, j, k, State::h) <= dry_depth) {
                a(i, j, k, State::hu) = amrex::Real(0.0);
                a(i, j, k, State::hv) = amrex::Real(0.0);
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
