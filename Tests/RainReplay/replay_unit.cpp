// The writer of flume.rain.replay_file (flume::RainReplay) on its own, with
// rain that a host model can hand over but that the coupled test does not
// reach: negative depths, a step of no rain, a step that starts an
// instant before the end of the last one (as a single-precision host gives
// it), and a restart while a negative depth is held. See
// Tests/RainReplay/README.md.
//
// Usage: replay_unit run | restart2 | restart3 | resume3 | early
//
//   run       start at 0 s, four host steps of 0.5 s (the fourth from
//             1.5 - 1e-8 s), the checkpoints chk2 and chk3 after the second
//             and the third step;
//   restart2  continue from chk2 with the third and fourth step (the third
//             starts after a step of no rain);
//   restart3  continue from chk3 with the fourth step (which takes the place
//             of the closing line of the third);
//   resume3   continue from chk3 and stop without a step;
//   early     start at 0 s, the first step, then a step from 0.25 s (before
//             the end of the first): this must stop.
//
// The grid is 4 x 2 cells of 1 m in boxes of 2 x 2 cells; the index is
// rain/replay.txt.

#include "FLUME_Checkpoint.H"
#include "FLUME_RainReplay.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <string>

namespace {

// The depth of rain, m, of each cell (i + 4 j) in each host step.
constexpr std::array<std::array<double, 8>, 4> depth_of_step = {{
    {1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, 5.0e-4, 0.0, 0.0},  // step 1: cell 1 holds 2e-3 back
    {0.0, 1.0e-3, 0.0, 0.0, 0.0, 0.0, -3.0e-4, 0.0},     // step 2: no cell has rain
    {0.0, 3.0e-3, 0.0, 2.0e-3, 0.0, 0.0, 1.0e-4, 0.0},   // step 3: cell 1 gives 2e-3, cell 6 holds
    {1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0e-4, 0.0},      // step 4: cell 6 gives 2e-4
}};
constexpr std::array<std::array<double, 2>, 4> time_of_step = {{
    {0.0, 0.5}, {0.5, 1.0}, {1.0, 1.5}, {1.5 - 1.0e-8, 2.0},
}};

void
step (flume::RainReplay& replay, amrex::MultiFab& depth, int n, double t0, double t1)
{
    const auto& d = depth.arrays();
    amrex::GpuArray<double, 8> row{};
    for (int c = 0; c < 8; ++c) {
        row[c] = depth_of_step[static_cast<std::size_t>(n)][static_cast<std::size_t>(c)];
    }
    amrex::ParallelFor(depth, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        d[b](i, j, k) = amrex::Real(row[i + 4 * j]);
    });
    amrex::Gpu::streamSynchronize();
    replay.record(depth, t0, t1);
}

} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode != "run" && mode != "restart2" && mode != "restart3" && mode != "resume3"
            && mode != "early") {
            amrex::Abort("usage: replay_unit run | restart2 | restart3 | resume3 | early");
        }
        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(AMREX_D_DECL(3, 1, 0)));
        const amrex::RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(4.0, 2.0, 1.0)});
        const amrex::Array<int, AMREX_SPACEDIM> periodic{AMREX_D_DECL(0, 0, 0)};
        const amrex::Geometry geom(domain, rb, 0, periodic);
        amrex::BoxArray ba(domain);
        ba.maxSize(2);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab depth(ba, dm, 1, 0);

        flume::RainReplay replay("rain/replay.txt", geom, ba, dm);
        int first = 0;
        if (mode == "restart2" || mode == "restart3" || mode == "resume3") {
            first = mode == "restart2" ? 2 : 3;
            const flume::CheckpointReader chk("chk" + std::to_string(first));
            replay.resume(chk);
            if (mode == "resume3") { first = 4; }
        } else {
            replay.start(0.0);
        }
        const int last = mode == "early" ? 1 : 4;
        for (int n = first; n < last; ++n) {
            const auto& t = time_of_step[static_cast<std::size_t>(n)];
            step(replay, depth, n, t[0], t[1]);
            if (mode == "run" && (n == 1 || n == 2)) {
                flume::CheckpointWriter chk("chk" + std::to_string(n + 1));
                replay.write_checkpoint(chk);
                chk.finish();
            }
        }
        if (mode == "early") { step(replay, depth, 2, 0.25, 0.75); }
        amrex::Print() << "replay_unit: " << mode << " done\n";
    }
    amrex::Finalize();
    return 0;
}
