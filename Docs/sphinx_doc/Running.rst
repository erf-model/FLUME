Running
=======

FLUME reads an AMReX inputs file, given as the first argument. Keys can be
added or overridden on the command line:

.. code-block:: bash

   mpiexec -n 4 ./flume inputs flume.max_grid_size=32

A first case
------------

``Tests/DamBreak/inputs`` holds still water 2 m and 1 m deep on the two sides
of a dam, on a flat bed between walls:

.. code-block:: text

   flume.n_cell        = 400 4
   flume.prob_hi       = 100.0 1.0
   flume.max_grid_size = 64
   flume.plot_file     = dam
   flume.stop_time     = 5.0

   flume.init.surface_elevation     = 2.0
   flume.init.dam_position          = 50.0
   flume.init.dam_surface_elevation = 1.0

.. code-block:: bash

   mpiexec -n 2 build-rel/flume Tests/DamBreak/inputs flume.order=2

The run writes ``dam00000`` (initial state) and the plotfile of the last
step. :doc:`Verification` compares this case with its exact solution.

A flood case adds a bed raster, buildings, friction and rain:

.. code-block:: text

   flume.n_cell    = 400 300
   flume.prob_lo   = 500000 4100000
   flume.prob_hi   = 502000 4101500
   flume.stop_time = 3600
   flume.order     = 2
   flume.plot_int  = 500

   flume.bed.file                = bed.asc
   flume.buildings.file          = buildings.asc
   flume.friction.manning_file   = manning.asc

   flume.rain.type = hyetograph
   flume.rain.file = storm.txt

   flume.bc.segments   = outlet
   flume.bc.outlet.side = xhi
   flume.bc.outlet.type = outflow

:doc:`Inputs` lists every key and :doc:`InputFiles` the file formats.

How a run ends
--------------

The run stops at ``flume.stop_time`` or after ``flume.max_step`` steps,
whichever comes first. With neither key it writes the initial state and takes
no step. The last step ends exactly on ``flume.stop_time``.

The run stops with a message naming the input when a key is unknown, out of
range or inconsistent with another, and when a depth becomes negative or not
finite (:ref:`sec:InputChecks`).

.. _sec:Output:

Output
------

The run writes AMReX plotfiles named ``<prefix><step>``, with the step in at
least five digits: the initial state, every ``flume.plot_int`` steps, and the
last state. An existing plotfile of the same name is first renamed to
``<prefix><step>.old.<number>``.

.. list-table:: Plotfile variables
   :header-rows: 1
   :widths: 24 16 60

   * - Variable
     - Unit
     - Meaning
   * - ``depth``
     - m
     - Water depth, :math:`h`.
   * - ``discharge_x``
     - m\ :sup:`2`/s
     - Discharge per unit width in x, :math:`hu`.
   * - ``discharge_y``
     - m\ :sup:`2`/s
     - Discharge per unit width in y, :math:`hv`.
   * - ``bed_elevation``
     - m
     - Bed elevation on the vertical datum, :math:`z_b`, buildings included.
   * - ``surface_elevation``
     - m
     - Water-surface elevation, :math:`z_b + h`; the bed elevation in a dry
       cell.
   * - ``infiltrated_depth``
     - m
     - Depth of water each cell has lost to the soil since the start
       (:ref:`sec:Infiltration`). Written only with
       ``flume.infiltration.type`` other than ``none``.

The plotfile is three-dimensional with one cell in z, from 0 to 1 m; the z
extent carries no meaning. ParaView, VisIt and yt read the format.

The log reports, at start-up, the boundary cells of each segment, the active
source terms, the rain, the infiltration parameters, the range of the Manning map and the
number of building cells.
