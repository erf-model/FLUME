Inputs
======

All of FLUME's keys are in the ``flume`` namespace. Every key has a range
check that stops the run with a message naming the key; no value is clamped.
An unknown key in the ``flume`` namespace stops the run at start-up. Lengths
are in metres, times in seconds, elevations on one vertical datum
(:ref:`sec:Datum`).

Grid
----

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.n_cell``
     - 2 integers
     - required
     - Cells in x and y. Each from 1 to 2\ :sup:`30`.
   * - ``flume.prob_lo``
     - 2 reals
     - ``0 0``
     - Lower corner of the domain. Finite.
   * - ``flume.prob_hi``
     - 2 reals
     - required
     - Upper corner. Finite and greater than ``flume.prob_lo``.
   * - ``flume.periodic``
     - 2 integers
     - ``0 0``
     - Periodic in x and in y. Each 0 or 1.
   * - ``flume.max_grid_size``
     - integer
     - ``64``
     - Largest box edge, in cells. From 1 to 2147483647.

In each direction the corners must be at most one eighth of the largest
floating-point number in magnitude, and the cell size must be a normal number
of at least machine epsilon times the larger corner (:ref:`sec:Precision`).

Time loop and scheme
--------------------

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.stop_time``
     - real
     - none
     - Time at which the run stops. Finite, at least 0. Held in double
       precision in every build.
   * - ``flume.max_step``
     - integer
     - none
     - Largest number of steps. From 0 to 2147483647.
   * - ``flume.cfl``
     - real
     - ``0.5``
     - CFL number (:doc:`theory/TimeStep`). Greater than 0 and at most 1; at
       most 0.5 with ``flume.order = 2``.
   * - ``flume.order``
     - integer
     - ``1``
     - Order of the scheme in space and time, 1 or 2
       (:doc:`theory/Discretization`). Order 2 costs about twice the work
       per step.
   * - ``flume.gravity``
     - real
     - ``9.81``
     - Acceleration of gravity, m/s\ :sup:`2`. Finite, greater than 0.
   * - ``flume.dry_depth``
     - real
     - ``1e-6``
     - Depth at or below which water counts as dry
       (:doc:`theory/WettingDrying`). Finite, at least 0. In single
       precision raise it for beds above about 10 m, where the bed is
       rounded to more than the default.

Output
------

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.plot_file``
     - string
     - ``plt``
     - Plotfile name prefix. Not empty.
   * - ``flume.plot_int``
     - integer
     - ``0``
     - Steps between plotfiles, from 0 to 2147483647. The first and the last
       state are always written.
   * - ``flume.log_int``
     - integer
     - ``1``
     - Steps between lines of the water-balance log ``water_balance.txt``
       and of the gauge files, from 0 to 2147483647; 0 writes no log. With
       a log, the first and the last step are always written
       (:ref:`sec:WaterBalanceLog`). Each line adds a reduction over the
       ranks; with outflow segments, every stage adds one reduction per
       segment; ``0`` adds none.
   * - ``flume.check_file``
     - string
     - ``chk``
     - Checkpoint name prefix. Not empty; needs ``flume.check_int > 0`` and
       must differ from ``flume.plot_file`` (:ref:`sec:Restart`).
   * - ``flume.check_int``
     - integer
     - ``0``
     - Steps between checkpoints, from 0 to 2147483647; 0 writes none. With
       checkpoints, the last step is always written (:ref:`sec:Restart`).
   * - ``flume.restart``
     - string
     - none
     - Checkpoint directory to continue from; one value, not empty. The
       grid, domain, periodic sides, sources, bed, water-balance columns,
       outflow segments, gauge cells and precision must be those of the run
       that wrote it; with a log (``flume.log_int > 0``) it must have been
       written with a log, and with ``flume.flood_maps = 1`` it must hold
       maps (:ref:`sec:Restart`).
   * - ``flume.flood_maps``
     - integer
     - ``0``
     - 1 adds the flood maps to the plotfiles (:ref:`sec:FloodMaps`); 0 or
       1.
   * - ``flume.flood_maps.arrival_depth``
     - real
     - ``flume.dry_depth``
     - Depth, m, above which water has arrived in a cell. Finite and
       :math:`\ge` ``flume.dry_depth`` (a shallower cell is dry); needs
       ``flume.flood_maps = 1``.
   * - ``flume.gauges.names``
     - strings
     - none
     - Names of the gauges, each of letters, digits and ``_``, none twice.
       Needs ``flume.log_int > 0`` (:ref:`sec:Gauges`).
   * - ``flume.gauges.<name>.location``
     - 2 reals
     - none
     - x and y of the gauge, m, inside the domain or on its edge. Required
       for each name in ``flume.gauges.names``; a key for a name not listed
       is unknown.

Boundaries
----------

A side of a direction that is not periodic is a wall, except where a segment
gives it another type (:doc:`theory/Boundaries`).

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.bc.segments``
     - strings
     - none
     - Names of the segments, each once; letters, digits and ``_``.
   * - ``flume.bc.<name>.side``
     - string
     - required
     - ``xlo``, ``xhi``, ``ylo`` or ``yhi``; not a periodic side.
   * - ``flume.bc.<name>.type``
     - string
     - required
     - ``wall`` or ``outflow``. ``discharge`` and ``stage`` are reserved and
       stop the run.
   * - ``flume.bc.<name>.range``
     - 2 reals
     - whole side
     - From and to along the side. Finite, ordered, inside the domain. The
       segment takes the boundary cells whose centre is at the first value
       or above and below the second: at least one, and no cell in two
       segments.

Bed and buildings
-----------------

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.bed.file``
     - string
     - none
     - ESRI ASCII grid of the bed elevation (:ref:`sec:RasterFiles`),
       interpolated bilinearly. Not with any other ``flume.bed`` key.
   * - ``flume.buildings.file``
     - string
     - none
     - ESRI ASCII grid of building heights above the bed, or a mask with
       ``flume.buildings.height``. Values at least 0; a cell above 0 is a
       building (:ref:`sec:Buildings`).
   * - ``flume.buildings.height``
     - real
     - none
     - Height given to every building cell of a mask. Finite, greater than
       0. Needs ``flume.buildings.file``.
   * - ``flume.bed.elevation``
     - real
     - ``0``
     - Bed elevation at the lower corner of the domain. Finite.
   * - ``flume.bed.slope``
     - 2 reals
     - ``0 0``
     - Slope of the plane bed in x and y (rise over run). Finite.
   * - ``flume.bed.continue_slope``
     - integer
     - ``0``
     - ``1`` continues the plane across periodic sides, an endless inclined
       channel; ``0`` repeats the bed. ``1`` needs a periodic direction with
       a slope.
   * - ``flume.bed.bowl_depth``
     - real
     - ``0``
     - Depth of a paraboloid bowl at the centre of the domain; negative for
       a mound. Finite. Needs ``flume.bed.bowl_radius``. For tests.
   * - ``flume.bed.bowl_radius``
     - real
     - none
     - Radius at which the bowl crosses the plane. Finite, greater than 0.
       Needs a ``flume.bed.bowl_depth`` that is not 0.
   * - ``flume.bed.noise``
     - real
     - ``0``
     - Amplitude of a fixed pseudo-random roughness of the bed. Finite, at
       least 0. For tests.

Friction
--------

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.friction.manning``
     - real
     - ``0``
     - Manning coefficient, s/m\ :sup:`1/3`, the same in every cell. Finite,
       at least 0; ``0`` is no friction. Typical: 0.01 to 0.02 smooth
       concrete, 0.03 natural channel, 0.05 to 0.1 vegetated floodplain.
   * - ``flume.friction.manning_file``
     - string
     - none
     - ESRI ASCII grid of the coefficient; each cell takes the raster cell
       that holds its centre. Values at least 0. Not with
       ``flume.friction.manning``.

Rain
----

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.rain.type``
     - string
     - ``none``
     - ``none``, ``constant``, ``hyetograph`` (:ref:`sec:Rain`),
       ``gridded`` or ``moving_cell`` (:ref:`sec:RainSpace`).
       ``microphysics`` is reserved and stops the run.
   * - ``flume.rain.rate``
     - real
     - required with ``constant``
     - Rain rate, mm/h. Finite, at least 0. Only with ``constant``.
   * - ``flume.rain.file``
     - string
     - required with ``hyetograph`` or ``gridded``
     - Hyetograph file (:ref:`sec:HyetographFiles`) or index of rain grids
       (:ref:`sec:RainIndexFiles`). Only with ``hyetograph`` or ``gridded``.
   * - ``flume.rain.interpolation``
     - string
     - ``step``
     - ``step`` holds each rate (or grid) until the next time; ``linear``
       interpolates in time. Only with ``hyetograph`` or ``gridded``.
   * - ``flume.rain.cell.peak_rate``
     - real
     - required with ``moving_cell``
     - Rate at the centre of the cell, mm/h. Finite, at least 0.
   * - ``flume.rain.cell.sigma``
     - real
     - required with ``moving_cell``
     - Standard deviation of the Gaussian, m. Finite, above 0.
   * - ``flume.rain.cell.center``
     - 2 reals
     - required with ``moving_cell``
     - Centre at :math:`t = 0`, m, in the coordinates of ``flume.prob_lo``.
       Both finite.
   * - ``flume.rain.cell.velocity``
     - 2 reals
     - ``0 0``
     - Velocity of the centre, m/s. Both finite. A step over which the
       centre moves more than 2500 ``sigma`` stops the run. The step bound
       stays the peak rate after the cell has left the domain, so a dry run
       needs ``flume.stop_time`` or ``flume.max_step``.

Each ``flume.rain.cell.*`` key is allowed only with ``moving_cell``;
``flume.rain.file`` and ``flume.rain.interpolation`` are not allowed with
it.

Infiltration
------------

The models are in :ref:`sec:Infiltration`. Each soil parameter is one value
for every cell or an ESRI ASCII grid (``_file``; each cell takes the raster
cell that holds its centre, and every cell of the domain must take a value in
the range), not both. A parameter is required with its model and allowed
only with it. Building cells (``flume.buildings.file``) infiltrate nothing.

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.infiltration.type``
     - string
     - ``none``
     - ``none``, ``green_ampt`` or ``curve_number``. ``curve_number`` needs
       rain (``flume.rain.type``).
   * - ``flume.infiltration.green_ampt.conductivity``,
       ``_conductivity_file``
     - real or string
     - required with ``green_ampt``
     - Saturated hydraulic conductivity :math:`K_s`, mm/h. Finite, at least
       0; 0 is impervious.
   * - ``flume.infiltration.green_ampt.suction``,
       ``_suction_file``
     - real or string
     - required with ``green_ampt``
     - Suction head at the wetting front :math:`\psi_s`, mm. Finite, at
       least 0.
   * - ``flume.infiltration.green_ampt.deficit``,
       ``_deficit_file``
     - real or string
     - required with ``green_ampt``
     - Initial moisture deficit :math:`\theta_d`, the porosity less the
       initial moisture content (a fraction). Finite, from 0 to 1.
   * - ``flume.infiltration.cn.number``, ``_number_file``
     - real or string
     - required with ``curve_number``
     - Runoff curve number :math:`CN`. Finite, above 0 and at most 100;
       100 is impervious. It must match ``ia_ratio`` (see below).
   * - ``flume.infiltration.cn.ia_ratio``
     - real
     - ``0.2``
     - Initial abstraction ratio :math:`\lambda = I_a / S`. Finite, at least
       0 and below 1. Only with ``curve_number``.

Initial state
-------------

Without ``flume.init.surface_elevation`` the domain starts dry
(:ref:`sec:InitialState`).

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.init.surface_elevation``
     - real
     - none
     - Elevation of the initial level surface. Finite. Cells whose bed is at
       or above it start dry.
   * - ``flume.init.velocity``
     - 2 reals
     - ``0 0``
     - Initial velocity where there is water, m/s. Finite. Needs
       ``flume.init.surface_elevation``.
   * - ``flume.init.dam_position``
     - real
     - none
     - Coordinate of a dam along ``flume.init.dam_direction``. Finite, with
       a cell centre on each side. Needs
       ``flume.init.dam_surface_elevation``.
   * - ``flume.init.dam_surface_elevation``
     - real
     - none
     - Surface elevation in the cells whose centre is at or above the dam
       position. Finite. Required with a dam.
   * - ``flume.init.dam_velocity``
     - 2 reals
     - ``flume.init.velocity``
     - Velocity on the high side of the dam, m/s. Finite. Needs a dam.
   * - ``flume.init.dam_direction``
     - ``x`` or ``y``
     - ``x``
     - Direction the dam position is measured along. Needs a dam.
   * - ``flume.init.surface_slope``
     - 2 reals
     - ``0 0``
     - Slope of the initial surface about the centre of the domain. Finite.
       Not with a dam. For tests.
   * - ``flume.init.surface_curvature``
     - real
     - ``0``
     - Curvature of the initial surface about the centre, 1/m. Finite. Not
       with a dam. For tests.
   * - ``flume.init.wave_amplitude``
     - real
     - ``0``
     - Amplitude of one sine wave of the depth, with the velocity of a
       simple wave. Finite. Needs a periodic direction along the wave; not
       with a dam, a surface slope or a curvature. For tests.
   * - ``flume.init.wave_direction``
     - ``x`` or ``y``
     - ``x``
     - Direction the wave varies along. Needs a wave.

.. _sec:InputChecks:

Checks
------

* A key given the wrong number of values stops the run.
* A value of the wrong type, or a key with no value, is stopped by AMReX's
  parser; its message may not name the key. The parser accepts expressions:
  ``1e3`` is a valid integer.
* A key in another namespace that nothing reads stops the run when it ends,
  after the output is written (``amrex.abort_on_unused_inputs``, which FLUME
  turns on; ``= 0`` turns it off).
* The bed and the initial state must be finite in every cell.
* Before every step and every plotfile the depth must be finite and not
  negative.
* The time step must be finite, above zero and large enough to change the
  time.
* A dry domain with ``flume.max_step`` alone and no rain ahead stops with a
  message asking for ``flume.stop_time`` (:ref:`sec:DryDomain`).

The integer keys have no joint limit on the number of boxes: a very large
``flume.n_cell`` with a small ``flume.max_grid_size`` can ask for more boxes
than memory holds.
