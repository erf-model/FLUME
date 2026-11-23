#include "FLUME_Gauges.H"
#include "FLUME_Checkpoint.H"
#include "FLUME_State.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

namespace flume {

namespace {
//! Digits that print an amrex::Real exactly.
constexpr int digits = std::numeric_limits<amrex::Real>::max_digits10;
} // namespace

amrex::IntVect
gauge_cell (const amrex::Array<amrex::Real, 2>& x, const amrex::Geometry& geom)
{
    amrex::IntVect iv(0, 0, 0);
    const amrex::Box& domain = geom.Domain();
    for (int d = 0; d < 2; ++d) {
        // In double: the point and the corner are within the domain, and
        // the quotient is at most the number of cells.
        const double lo = static_cast<double>(geom.ProbLo(d));
        const double hi = static_cast<double>(geom.ProbHi(d));
        const int n = domain.length(d);
        const double f = (static_cast<double>(x[d]) - lo) / (hi - lo) * static_cast<double>(n);
        // A point within rounding of a face is on the face: 4 machine
        // epsilons of the build of the largest coordinate involved, in
        // cells, and never more than 1e-3 of a cell (with coordinates far
        // from 0 in single precision the rounding is a large part of a
        // cell, and a point in the upper part of a cell would be moved up).
        const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
        const double big = std::max({std::abs(static_cast<double>(x[d])), std::abs(lo), std::abs(hi)});
        const double tol = std::min(4.0 * eps * big / (hi - lo) * static_cast<double>(n), 1e-3);
        const double r = std::round(f);
        const int c = (std::abs(f - r) <= tol) ? static_cast<int>(r)
                                               : static_cast<int>(std::floor(f));
        iv[d] = std::clamp(c, 0, n - 1);
    }
    return iv;
}

Gauges::Gauges (const Inputs& in, const amrex::Geometry& geom) : m_dry_depth(in.dry_depth)
{
    for (const auto& g : in.gauges) {
        m_names.push_back(g.name);
        m_locations.push_back(g.location);
        m_cells.push_back(gauge_cell(g.location, geom));
        amrex::Print() << "FLUME: gauge " << g.name << " at (" << g.location[0] << ", "
                       << g.location[1] << ") m reads cell (" << m_cells.back()[0] << ", "
                       << m_cells.back()[1] << ")\n";
    }
}

void
Gauges::start ()
{
    open(-1, false);
}

void
Gauges::resume (int step, bool keep_step)
{
    open(step, keep_step);
}

std::vector<std::string>
Gauges::words () const
{
    std::vector<std::string> w;
    for (std::size_t n = 0; n < m_names.size(); ++n) {
        w.push_back(m_names[n] + ":" + std::to_string(m_cells[n][0]) + ":"
                    + std::to_string(m_cells[n][1]));
    }
    return w;
}

void
Gauges::open (int resume_step, bool keep_step)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    m_files.clear();
    for (std::size_t n = 0; n < m_names.size(); ++n) {
        const std::string file = "gauge_" + m_names[n] + ".txt";
        if (resume_step >= 0 && rewind_log(file, resume_step, keep_step)) {
            auto f = std::make_unique<std::ofstream>(file, std::ios::app);
            if (!*f) { amrex::Abort("FLUME: cannot open " + file + " for writing"); }
            m_files.push_back(std::move(f));
            continue;
        }
        auto f = std::make_unique<std::ofstream>(file);
        if (!*f) { amrex::Abort("FLUME: cannot open " + file + " for writing"); }
        *f << "# FLUME gauge " << m_names[n] << " at x = " << std::setprecision(digits)
           << m_locations[n][0] << " m, y = " << m_locations[n][1] << " m, cell ("
           << m_cells[n][0] << ", " << m_cells[n][1] << ")\n"
           << "# step time depth surface_elevation velocity_x velocity_y\n";
        f->flush();
        m_files.push_back(std::move(f));
    }
}

void
Gauges::log (const amrex::MultiFab& state, const amrex::MultiFab& bed, int step, double time)
{
    if (m_names.empty()) { return; }
    constexpr int nv = 4;  // depth, bed, hu, hv
    const int ng = static_cast<int>(m_names.size());
    // Only the rank that owns a gauge's cell writes its values; the others
    // hold 0, so the sum over the ranks is the value.
    amrex::Gpu::DeviceVector<amrex::Real> d_vals(static_cast<std::size_t>(nv * ng),
                                                 amrex::Real(0.0));
    amrex::Real* p = d_vals.data();
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const amrex::Box& vb = mfi.validbox();
        const auto& a = state.const_array(mfi);
        const auto& z = bed.const_array(mfi);
        for (int g = 0; g < ng; ++g) {
            const amrex::IntVect iv = m_cells[static_cast<std::size_t>(g)];
            if (!vb.contains(iv)) { continue; }
            amrex::ParallelFor(amrex::Box(iv, iv),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    p[nv * g + 0] = a(i, j, k, State::h);
                    p[nv * g + 1] = z(i, j, k);
                    p[nv * g + 2] = a(i, j, k, State::hu);
                    p[nv * g + 3] = a(i, j, k, State::hv);
                });
        }
    }
    amrex::Vector<amrex::Real> vals(static_cast<std::size_t>(nv * ng));
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_vals.begin(), d_vals.end(), vals.begin());
    amrex::Gpu::streamSynchronize();
    const int io = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::ReduceRealSum(vals.data(), nv * ng, io);
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    for (int g = 0; g < ng; ++g) {
        const amrex::Real h = vals[nv * g + 0];
        const amrex::Real zb = vals[nv * g + 1];
        const bool wet = h > m_dry_depth;
        const amrex::Real u = wet ? vals[nv * g + 2] / h : amrex::Real(0.0);
        const amrex::Real v = wet ? vals[nv * g + 3] / h : amrex::Real(0.0);
        auto& f = *m_files[static_cast<std::size_t>(g)];
        f << step << " " << std::setprecision(17) << time << " " << std::setprecision(digits) << h
          << " " << zb + h << " " << u << " " << v << "\n";
        f.flush();
    }
}

} // namespace flume
