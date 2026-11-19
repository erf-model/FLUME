#include "FLUME_Inputs.H"
#include "FLUME_Run.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>

namespace {

// A key that nothing reads stops the run when it ends, unless the user sets
// amrex.abort_on_unused_inputs = 0. This catches a misspelled namespace
// (flum.n_cell) and AMReX keys FLUME does not use (geometry.prob_lo), which
// the check of the flume namespace cannot see.
void
set_amrex_defaults ()
{
    amrex::ParmParse pp("amrex");
    if (!pp.contains("abort_on_unused_inputs")) {
        pp.add("abort_on_unused_inputs", 1);
    }
}

} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, set_amrex_defaults);
    {
        const flume::Inputs in = flume::read_inputs();
        flume::run(in);
    }
    amrex::Finalize();
    return 0;
}
