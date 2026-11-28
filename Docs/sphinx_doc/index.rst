FLUME: FLood and Urban inundation ModulE
=========================================

FLUME is a two-dimensional shallow-water flood solver for the
`Energy Research and Forecasting (ERF) <https://github.com/erf-model/ERF>`__
atmospheric code, built on `AMReX <https://github.com/AMReX-Codes/amrex>`__.
It solves the full shallow-water equations over a fixed bed with a
well-balanced finite-volume scheme, wetting and drying, Manning friction and
prescribed rain. It runs with MPI and compiles for CUDA, HIP and SYCL.

FLUME runs stand-alone, using only AMReX, taken from the commit ERF pins, or
inside ERF as its flood engine (``FLUME_WITH_ERF``, :ref:`sec:InsideERF`),
where its rain can be the surface precipitation of ERF's microphysics
(``flume.rain.type = microphysics``). The expansion of the name is a draft.

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
     - Constant rate, hyetograph, gridded rain, moving storm cell; inside
       ERF, the surface precipitation of ERF's microphysics
     -
   * - Terrain
     - Bed, buildings and Manning map from ESRI ASCII grids
     - Other raster formats through GDAL
   * - Losses
     - Infiltration: Green-Ampt and the runoff curve number; storm-drain
       inlets (weir and orifice capacity, or capacity from the approach
       flow), one by one or from a list file
     -
   * - Output
     - AMReX plotfiles, flood maps, water-balance log, gauge time series,
       checkpoint and restart
     -

.. toctree::
   :caption: USER GUIDE
   :maxdepth: 1

   Building
   Running
   CoupledExample
   TownExample
   GowanusExample
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
   theory/WaterBalance
   theory/TimeStep

.. toctree::
   :caption: VERIFICATION
   :maxdepth: 1

   Verification
   TestCatalog
   References
