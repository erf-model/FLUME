#include "FLUME_Run.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Utility.H>

namespace flume {

void
run (const Inputs& in)
{
    // One cell in z: the flood grid is two-dimensional.
    const amrex::Box domain(amrex::IntVect(0, 0, 0),
                            amrex::IntVect(in.n_cell[0] - 1, in.n_cell[1] - 1, 0));
    const amrex::RealBox real_box({in.prob_lo[0], in.prob_lo[1], amrex::Real(0.0)},
                                  {in.prob_hi[0], in.prob_hi[1], amrex::Real(1.0)});
    const amrex::Array<int, 3> is_periodic{{0, 0, 0}};
    const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(in.max_grid_size, in.max_grid_size, 1));
    const amrex::DistributionMapping dm(ba);

    // Water depth, m. No ghost cells: nothing reads a neighbour yet.
    amrex::MultiFab depth(ba, dm, 1, 0);
    depth.setVal(amrex::Real(0.0));

    const int step = 0;
    const amrex::Real time = amrex::Real(0.0);
    const std::string name = amrex::Concatenate(in.plot_file, step, 5);
    amrex::WriteSingleLevelPlotfile(name, depth, {"depth"}, geom, time, step);

    amrex::Print() << "FLUME skeleton: " << in.n_cell[0] << " x " << in.n_cell[1]
                   << " cells in " << ba.size() << " box(es); wrote " << name << "\n";
}

} // namespace flume
