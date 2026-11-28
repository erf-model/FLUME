#include "FLUME_Inputs.H"
#include "FLUME_Flood.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <ostream>
#include <string>

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

// The stop itself, made here and not left to AMReX: AMReX prints the names of
// the unused keys only when the ParmParse verbosity is above 0, so with
// amrex.verbose = 0 its abort would not say which key it is. The list goes to
// stdout and into the abort message (stderr), so a log of either names the
// keys. As in AMReX, the I/O rank decides: a key read on that rank alone is
// not unused.
void
abort_on_unused_inputs ()
{
    bool stop = true;
    amrex::ParmParse pp("amrex");
    pp.query("abort_on_unused_inputs", stop);
    if (!stop || !amrex::ParallelDescriptor::IOProcessor() || !amrex::ParmParse::hasUnusedInputs()) {
        return;
    }
    std::string msg = "Unused key(s): nothing read them (amrex.abort_on_unused_inputs = 0"
                      " turns this check off):";
    for (const auto& s : amrex::ParmParse::getUnusedInputs()) {
        msg += "\n  ";
        msg += s;
    }
    amrex::OutStream() << msg << std::endl;
    amrex::Abort(msg);
}

} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, set_amrex_defaults);
    // AMReX reads amrex.parmparse.verbose (and .v) the first time this is
    // called. Left to amrex::Finalize(), that is after it has looked for
    // unused keys: the key would count as unused and stop every run given it.
    amrex::ignore_unused(amrex::ParmParse::Verbose());
    // Run by amrex::Finalize(), where AMReX makes its own check: after the run,
    // since AMReX reads some of its keys only when they are first needed
    // (DistributionMapping.strategy), and after the profilers have written
    // their output (tiny_profiler.output_file is read there). Functions given
    // to ExecOnFinalize run last-in first-out, so this one runs before
    // ParmParse::Finalize, which amrex::Initialize registered.
    amrex::ExecOnFinalize(abort_on_unused_inputs);
    {
        const flume::Inputs in = flume::read_inputs();
        flume::Flood flood(in);
        flood.evolve();
    }
    amrex::Finalize();
    return 0;
}
