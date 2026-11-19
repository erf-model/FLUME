FLUME
=====

FLUME is a two-dimensional shallow-water flood solver for the
`Energy Research and Forecasting (ERF) <https://github.com/erf-model/ERF>`__
atmospheric code, built on `AMReX <https://github.com/AMReX-Codes/amrex>`__.
The expansion of the name, FLood and Urban inundation ModulE, is a draft.

**Status.** The repository holds the build system, the test harness, the
continuous integration and this documentation. The executable reads an inputs
file, builds the grid, the bed and the initial state, and runs a time loop with
a CFL time step that writes plotfiles. The loop does not change the state yet:
the flux kernel is not written.

The plan is a stand-alone code first, which uses AMReX only: the solver,
prescribed rainfall, losses, restart and benchmarks. A build with ERF, with
rainfall from ERF's microphysics, comes after that. ERF is
already a submodule of this repository, because the stand-alone build takes
AMReX from the commit ERF pins.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   Building
   Inputs
   Model
   Testing
