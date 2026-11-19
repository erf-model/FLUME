#include "FLUME_Advance.H"
#include "FLUME_Riemann.H"
#include "FLUME_State.H"

#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFParallelFor.H>

namespace flume {

void
compute_fluxes (amrex::Array<amrex::MultiFab, 2>& flux, const amrex::MultiFab& state,
                const amrex::MultiFab& bed, amrex::Real gravity)
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
                // The velocity of a cell with no water is not defined: zero.
                const amrex::Real zero = amrex::Real(0.0);
                const amrex::Real unL = (hL > zero) ? a(L, i_norm) / hL : zero;
                const amrex::Real utL = (hL > zero) ? a(L, i_tang) / hL : zero;
                const amrex::Real unR = (hR > zero) ? a(R, i_norm) / hR : zero;
                const amrex::Real utR = (hR > zero) ? a(R, i_tang) / hR : zero;
                const FaceFlux F = hllc_flux(face_depth(hL, zL, zR), unL, utL,
                                             face_depth(hR, zR, zL), unR, utR, gravity);
                f[box_no](i, j, k, State::h) = F.mass;
                f[box_no](i, j, k, i_norm) = F.normal;
                f[box_no](i, j, k, i_tang) = F.tangent;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void
update_state (amrex::MultiFab& state, const amrex::Array<amrex::MultiFab, 2>& flux,
              const amrex::MultiFab& bed, const amrex::Geometry& geom, amrex::Real gravity,
              amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& zb = bed.const_arrays();
    const auto& fx = flux[0].const_arrays();
    const auto& fy = flux[1].const_arrays();
    const amrex::Real half_g = amrex::Real(0.5) * gravity;
    const amrex::Real dtdx = dt / geom.CellSize(0);
    const amrex::Real dtdy = dt / geom.CellSize(1);
    // Both directions in one kernel, so the bed-slope term of each uses the
    // depth from before the step. A cell reads its own state, the bed of its
    // four neighbours and the fluxes on its four faces, and writes its own
    // state: the update is in place and needs no atomics.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const auto& a = s[box_no];
            const auto& z = zb[box_no];
            const auto& f = fx[box_no];
            const auto& g = fy[box_no];
            const amrex::Real h = a(i, j, k, State::h);
            const amrex::Real zc = z(i, j, k);
            const amrex::Real h_xlo = face_depth(h, zc, z(i - 1, j, k));
            const amrex::Real h_xhi = face_depth(h, zc, z(i + 1, j, k));
            const amrex::Real h_ylo = face_depth(h, zc, z(i, j - 1, k));
            const amrex::Real h_yhi = face_depth(h, zc, z(i, j + 1, k));
            for (int n = 0; n < State::ncomp; ++n) {
                a(i, j, k, n) -= dtdx * (f(i + 1, j, k, n) - f(i, j, k, n))
                                 + dtdy * (g(i, j + 1, k, n) - g(i, j, k, n));
            }
            // h_hi^2 - h_lo^2 as a product: exactly zero for equal depths,
            // also where the compiler contracts a * a - b * b into a fused
            // multiply-add (which would leave the rounding error of b * b).
            a(i, j, k, State::hu) += dtdx * half_g * ((h_xhi - h_xlo) * (h_xhi + h_xlo));
            a(i, j, k, State::hv) += dtdy * half_g * ((h_yhi - h_ylo) * (h_yhi + h_ylo));
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
