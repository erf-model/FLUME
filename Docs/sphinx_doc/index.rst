FLUME: FLood and Urban inundation ModulE
=========================================

FLUME is a two-dimensional shallow-water flood solver for the
`Energy Research and Forecasting (ERF) <https://github.com/erf-model/ERF>`__
atmospheric code, built on `AMReX <https://github.com/AMReX-Codes/amrex>`__.
It solves the full shallow-water equations over a fixed bed with a
well-balanced finite-volume scheme, wetting and drying, Manning friction and
prescribed rain. It runs with MPI and compiles for CUDA, HIP and SYCL.

FLUME runs stand-alone today and uses only AMReX, taken from the commit ERF
pins. A build inside ERF, with rain from ERF's microphysics, is planned. The
expansion of the name is a draft.

.. list-table:: Capabilities
   :header-rows: 1
   :widths: 22 42 36

   * - Area
     - Available
     - Planned
   * - Scheme
     - HLLC flux with hydrostatic reconstruction; first and second order
     -
   * - Wet and dry ground
     - Moving shorelines, non-negative depth, mass conserved
     -
   * - Boundaries
     - Walls, periodic sides, outflow segments
     - Prescribed discharge and stage
   * - Friction
     - Manning, one coefficient or a map
     -
   * - Rain
     - Constant rate, hyetograph, gridded rain, moving storm cell
     - ERF microphysics
   * - Terrain
     - Bed, buildings and Manning map from ESRI ASCII grids
     - Other raster formats through GDAL
   * - Losses
     -
     - Infiltration, storm-drain inlets
   * - Output
     - AMReX plotfiles
     - Checkpoint and restart, flood maps, water-balance log

.. toctree::
   :caption: USER GUIDE
   :maxdepth: 1

   Building
   Running
   Inputs
   InputFiles
   Testing

.. toctree::
   :caption: THEORY
   :maxdepth: 1

   theory/Equations
   theory/Terrain
   theory/Discretization
   theory/WettingDrying
   theory/SourceTerms
   theory/Boundaries
   theory/TimeStep

.. toctree::
   :caption: VERIFICATION
   :maxdepth: 1

   Verification
   TestCatalog
   References
