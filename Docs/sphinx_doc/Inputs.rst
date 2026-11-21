Inputs
======

FLUME reads an AMReX inputs file, given as the first argument. Keys can also be
given or overridden on the command line:

.. code-block:: bash

   mpiexec -n 4 ./flume inputs flume.max_grid_size=32

All of FLUME's keys are in the ``flume`` namespace. Every key has a range check
that stops the run with a message naming the key; no value is clamped. A key in
the ``flume`` namespace that FLUME does not know also stops the run at start-up,
so a misspelled key is not ignored.

A key in any other namespace that nothing reads (``flum.n_cell``, or an AMReX
key FLUME does not use, such as ``geometry.prob_lo``) is listed and stops the
run when it ends, after the output is written. FLUME turns on AMReX's
``amrex.abort_on_unused_inputs`` for this; ``amrex.abort_on_unused_inputs = 0``
turns it off.

A value of the wrong type (text for a number, ``1.5`` for an integer) is
stopped by AMReX's parser. Its message shows the value and may not name the
key. The parser accepts expressions, so ``1e3`` is a valid integer and ``16/2``
a valid number.

.. list-table::
   :header-rows: 1
   :widths: 24 12 14 50

   * - Key
     - Type
     - Default
     - Meaning and allowed values
   * - ``flume.n_cell``
     - 2 integers
     - required
     - Number of cells in x and y. Each must be from 1 to 1073741824
       (2\ :sup:`30`).
   * - ``flume.prob_lo``
     - 2 reals
     - ``0 0``
     - Lower corner of the domain in x and y, in metres. Must be finite.
   * - ``flume.prob_hi``
     - 2 reals
     - required
     - Upper corner of the domain in x and y, in metres. Must be finite and
       greater than ``flume.prob_lo`` in both directions. See the two checks
       below the table.
   * - ``flume.periodic``
     - 2 integers
     - ``0 0``
     - Whether the domain is periodic in x and in y. Each must be 0 or 1. A
       direction with 0 has two sides made of the segments below, walls where
       no segment covers them.
   * - ``flume.bc.segments``
     - strings
     - none
     - Names of the segments of the sides that are not periodic, at least
       one. A name holds only letters, digits and ``_`` and is listed once.
       Without this key every side that is not periodic is a wall
       (:doc:`Model`). The keys of a name that is not listed are unknown.
   * - ``flume.bc.<name>.side``
     - string
     - required
     - The side of segment ``<name>``: ``xlo``, ``xhi``, ``ylo`` or ``yhi``.
       Must not be a side of a periodic direction.
   * - ``flume.bc.<name>.type``
     - string
     - required
     - ``wall`` or ``outflow`` (zero gradient: :doc:`Model`). ``discharge``
       and ``stage`` are kept for later versions and stop the run.
   * - ``flume.bc.<name>.range``
     - 2 reals
     - the whole side
     - From and to, in metres, along the side (y on ``xlo`` and ``xhi``, x
       on ``ylo`` and ``yhi``). Must be finite, the first below the second,
       both within ``flume.prob_lo`` and ``flume.prob_hi``. The segment takes
       the boundary cells whose centre is at the first value or above and
       below the second; it must hold at least one, and no cell may belong
       to two segments. The run prints the cells of each segment.
   * - ``flume.max_grid_size``
     - integer
     - ``64``
     - Largest box edge, in cells, in x and y. Must be from 1 to 2147483647.
   * - ``flume.plot_file``
     - string
     - ``plt``
     - Plotfile name prefix. Must not be empty. The plotfile of step 0 is
       ``<prefix>00000``.
   * - ``flume.plot_int``
     - integer
     - ``0``
     - Steps between plotfiles. Must be from 0 to 2147483647. The first and
       the last state are always written; ``0`` writes only those.
   * - ``flume.stop_time``
     - real
     - none
     - Time at which the run stops, in seconds. Must be finite and at least 0.
       The last step ends at this time. Read in double precision in every
       build.
   * - ``flume.max_step``
     - integer
     - none
     - Largest number of steps. Must be from 0 to 2147483647.
   * - ``flume.cfl``
     - real
     - ``0.5``
     - CFL number of the time step, with the wave rates of the two directions
       added (:doc:`Model`). Must be greater than 0 and at most 1, and at
       most 0.5 with ``flume.order = 2``. For any
       value in the range the water that leaves a cell in a step is limited
       to the water the cell holds (:doc:`Model`).
   * - ``flume.order``
     - integer
     - ``1``
     - Order of the scheme in space and time. Must be 1 or 2. ``1`` is the
       first-order step. ``2`` reconstructs the depth, the surface and the
       velocity in each cell with a minmod limiter and takes two stages per
       step (:doc:`Model`); it reads two cells on each side of a cell and
       costs about twice the work of ``1`` per step.
   * - ``flume.gravity``
     - real
     - ``9.81``
     - Acceleration of gravity, in m/s\ :sup:`2`. Must be finite and greater
       than 0.
   * - ``flume.dry_depth``
     - real
     - ``1e-6``
     - Depth at or below which water counts as dry, in metres. Must be finite
       and at least 0. A cell at or below it has no velocity, holds no
       discharge after a step (the initial state is kept as given) and does not enter the time step; a face that
       sees no more than this depth of a cell sees none (:doc:`Model`). The
       water of such a cell is kept. ``0`` treats every depth above zero as
       wet. In single precision the bed elevation is rounded to about
       1e-6 m at elevations of 10 m, the default: raise the value for
       higher beds.
   * - ``flume.friction.manning``
     - real
     - ``0``
     - Manning coefficient of the bed, in s/m\ :sup:`1/3`, the same in every
       cell. Must be finite and at least 0. ``0`` is no friction. Typical
       values: 0.01 to 0.02 for smooth concrete, about 0.03 for a natural
       channel, 0.05 to 0.1 for a floodplain with vegetation. The update is
       implicit in the discharge and sets no limit on the time step
       (:doc:`Model`). Not with ``flume.friction.manning_file``.
   * - ``flume.friction.manning_file``
     - string
     - none
     - ESRI ASCII grid of the Manning coefficient, in s/m\ :sup:`1/3`
       (`Raster files`_). Must not be empty. Each cell takes the value of the
       raster cell that holds its centre; every value read must be at least
       0. Not with ``flume.friction.manning``. The run prints the smallest
       and the largest coefficient.
   * - ``flume.rain.type``
     - string
     - ``none``
     - ``none``, ``constant`` or ``hyetograph``: no rain, rain at a constant
       rate, or a rate that changes in time from a file. The rain is the same
       in every cell and adds depth without momentum (:doc:`Model`).
       ``gridded``, ``moving_cell`` and ``microphysics`` are kept for later
       versions and stop the run. The run prints the rain at start-up.
   * - ``flume.rain.rate``
     - real
     - required with ``constant``
     - Rain rate, in mm/h. Must be finite and at least 0. Only with
       ``flume.rain.type = constant``.
   * - ``flume.rain.file``
     - string
     - required with ``hyetograph``
     - Hyetograph file (`Hyetograph files`_). Must not be empty. Only with
       ``flume.rain.type = hyetograph``. A relative path is taken from the
       directory the run starts in.
   * - ``flume.rain.interpolation``
     - string
     - ``step``
     - ``step`` holds each rate of the file until the next time; ``linear``
       interpolates linearly in time between the entries. Before the first
       time and after the last the rate is the first or the last one, in
       both modes. Only with ``flume.rain.type = hyetograph``.
   * - ``flume.bed.file``
     - string
     - none
     - ESRI ASCII grid of the bed elevation, in metres on the vertical datum
       (`Raster files`_). Must not be empty. Each cell takes the raster
       interpolated bilinearly between raster cell centres at its centre.
       Not with ``flume.bed.elevation``, ``flume.bed.slope``,
       ``flume.bed.continue_slope``, ``flume.bed.bowl_depth``,
       ``flume.bed.bowl_radius`` or ``flume.bed.noise``: the file gives the
       whole bed.
   * - ``flume.buildings.file``
     - string
     - none
     - ESRI ASCII grid of building heights above the bed, in metres, or a
       mask with ``flume.buildings.height`` (`Raster files`_). Must not be
       empty. Each cell takes the value of the raster cell that holds its
       centre; every value read must be at least 0. A cell with a value
       above 0 has a building: the bed is raised by the value (or by
       ``flume.buildings.height``), and water enters it only over the top
       (:doc:`Model`). The run prints the number of building cells.
   * - ``flume.buildings.height``
     - real
     - none
     - Height, in metres, given to every cell whose ``flume.buildings.file``
       value is above 0: the file is then a mask (0 or 1, say). Must be
       finite and greater than 0. Needs ``flume.buildings.file``.
   * - ``flume.bed.elevation``
     - real
     - ``0``
     - Bed elevation at the lower corner of the domain, in metres on the
       vertical datum. Must be finite.
   * - ``flume.bed.slope``
     - 2 reals
     - ``0 0``
     - Slope of the bed in x and y (rise over run). Must be finite. The bed is
       a plane: see :doc:`Model`.
   * - ``flume.bed.continue_slope``
     - integer
     - ``0``
     - ``0``: across a periodic side the bed repeats, with a step of the
       slope times the length of the domain. ``1``: the plane of the bed
       goes on across the periodic sides, an endless inclined channel
       (:doc:`Model`). Must be 0 or 1; ``1`` needs ``flume.periodic = 1`` in
       a direction where ``flume.bed.slope`` is not 0.
   * - ``flume.bed.bowl_depth``
     - real
     - ``0``
     - Depth of a paraboloid bowl added to the bed, in metres, at the centre
       of the domain. Must be finite. A negative value makes a mound; ``0``
       is no bowl. Needs ``flume.bed.bowl_radius`` unless it is 0
       (:doc:`Model`). For tests.
   * - ``flume.bed.bowl_radius``
     - real
     - none
     - Distance from the centre of the domain at which the bowl crosses the
       plane bed, in metres. Must be finite and greater than 0, with a
       square that is a normal number. Required
       with a ``flume.bed.bowl_depth`` that is not 0, and needs one.
   * - ``flume.bed.noise``
     - real
     - ``0``
     - Amplitude of a rough bed, in metres. Must be finite and at least 0.
       Each cell's bed moves by a fixed pseudo-random amount between minus
       and plus this value (:doc:`Model`). For tests.
   * - ``flume.init.surface_elevation``
     - real
     - none
     - Elevation of the initial, level water surface, in metres on the same
       datum as the bed. Must be finite. Cells whose bed is at or above it
       start dry. Without this key the whole domain starts dry. With a dam it
       is the surface on the low-coordinate side.
   * - ``flume.init.surface_slope``
     - 2 reals
     - ``0 0``
     - Slope of the initial water surface in x and y (rise over run). Must be
       finite. ``flume.init.surface_elevation`` is then the elevation at the
       centre of the domain (:doc:`Model`). Needs
       ``flume.init.surface_elevation``; not with a dam. For tests.
   * - ``flume.init.surface_curvature``
     - real
     - ``0``
     - Curvature of the initial water surface, in 1/m: the surface rises by
       this value times the square of the distance from the centre of the
       domain. Must be finite. Needs ``flume.init.surface_elevation``; not
       with a dam. For tests.
   * - ``flume.init.wave_amplitude``
     - real
     - ``0``
     - Amplitude of one sine wave of the initial depth along
       ``flume.init.wave_direction``, in metres, with the velocity of a
       simple wave (:doc:`Model`). Must be finite; ``0`` is no wave. Needs
       ``flume.init.surface_elevation`` and ``flume.periodic = 1`` along the
       wave; not with a dam, ``flume.init.surface_slope`` or
       ``flume.init.surface_curvature``. The state is a simple wave only
       over a flat bed. For tests.
   * - ``flume.init.wave_direction``
     - ``x`` or ``y``
     - ``x``
     - Direction along which the wave varies. Needs a
       ``flume.init.wave_amplitude`` that is not 0.
   * - ``flume.init.velocity``
     - 2 reals
     - ``0 0``
     - Initial velocity in x and y where there is water, in m/s. Must be
       finite. Needs ``flume.init.surface_elevation``. With a dam it is the
       velocity on the low-coordinate side.
   * - ``flume.init.dam_position``
     - real
     - none
     - Coordinate of a dam, in metres, along the direction
       ``flume.init.dam_direction``. Must be finite and leave a cell centre
       on each side: greater than the first cell centre and at most the last
       one in that direction. Along a periodic direction the two sides also
       meet at the side of the domain, which makes a second dam there.
       Needs ``flume.init.surface_elevation`` and
       ``flume.init.dam_surface_elevation``. Without this key there is no dam.
   * - ``flume.init.dam_surface_elevation``
     - real
     - none
     - Elevation of the initial water surface in the cells whose centre is at
       or above the dam position, in metres on the same datum. Must be
       finite. Required with ``flume.init.dam_position``, and needs it.
   * - ``flume.init.dam_velocity``
     - 2 reals
     - ``flume.init.velocity``
     - Initial velocity in x and y in the cells whose centre is at or above
       the dam position, in m/s. Must be finite. Needs
       ``flume.init.dam_position``.
   * - ``flume.init.dam_direction``
     - ``x`` or ``y``
     - ``x``
     - Direction along which the dam position is measured: the dam lies
       across it. Needs ``flume.init.dam_position``.

The keys that take two values stop the run when given any other number of
values, and so do the keys that take one.

Two more checks use ``flume.prob_lo``, ``flume.prob_hi`` and ``flume.n_cell``
together, in each direction:

* The corners must be at most one eighth of the largest floating-point number
  in magnitude (2.2e307 in double precision, 4.2e37 in single).
* The cell size, ``(prob_hi - prob_lo) / n_cell``, must be a normal
  floating-point number and at least machine epsilon times the larger corner
  in magnitude. Below that, neighbouring cell edges are the same number. In
  double precision the limit is about 1e-9 m at coordinates of 4.4e6 m. In
  single precision it is about 0.5 m there.

The run stops at ``flume.stop_time`` or after ``flume.max_step`` steps,
whichever comes first. With neither key the run takes no step and writes the
initial state only.

Four checks are made on the fields, after the inputs are read, because finite
inputs can give a number that is not finite:

* The bed elevation must be finite in every cell (``flume.bed.slope`` times
  the size of the domain, ``flume.bed.noise`` added to the bed, or a bowl of
  a small radius can pass the largest floating-point number; in single
  precision so can a value of ``flume.bed.file`` or ``flume.buildings.file``).
* The initial depth and discharges must be finite in every cell.
* Before every step and before every plotfile the depth must be finite and
  not negative, and the wave speed finite, in every cell. The step limits
  the water that leaves a cell to the water it holds, so a negative depth
  after the first step is not expected for any ``flume.cfl``; the check
  stays as a guard.
* At every step the time step must be finite, above zero and large enough to
  change the time.

A domain with no depth above ``flume.dry_depth`` has no wave to set the time
step, and the run prints, once for each such spell, that the fluxes are not
computed; the discharges of its cells are set to zero, as a step with fluxes
does in a dry cell. With rain ahead the rain sets the step (:doc:`Model`). Without rain the run takes one step to
``flume.stop_time``; with ``flume.max_step`` alone it stops with a message
asking for ``flume.stop_time``, after the initial plotfile is written. This
also holds later in a run, if the water spreads until no cell is above the
tolerance.

A key with no value (``flume.cfl =``) is stopped by AMReX's parser, with its
own message.

The integer keys have no joint limit on the number of boxes. A very large
``flume.n_cell`` with a small ``flume.max_grid_size`` can ask for more boxes
than memory holds.

Raster files
------------

``flume.bed.file``, ``flume.buildings.file`` and ``flume.friction.manning_file``
name ESRI ASCII grids (``.asc``), the plain-text raster format that GIS tools
such as QGIS and ArcGIS write. A relative name is taken from the directory the
run starts in. A name with spaces goes in double quotes in the inputs file
(``flume.bed.file = "survey 2024/bed.asc"``). The name must be a regular file,
not a directory. The file is a header followed by the values:

.. code-block:: text

   ncols         4
   nrows         3
   xllcorner     500000.0
   yllcorner     4100000.0
   cellsize      10.0
   NODATA_value  -9999
   12.1 12.4 12.9 13.3
   11.8 12.0 12.5 12.8
   11.6 11.7 12.1 12.4

* The header keywords may be in any case and in any order, each once.
  ``ncols`` and ``nrows`` are whole numbers from 1 to 2\ :sup:`30`. The lower
  left corner is ``xllcorner`` and ``yllcorner``, or the centre of the lower
  left cell, ``xllcenter`` and ``yllcenter``. The cell size is ``cellsize``,
  or ``dx`` and ``dy`` for cells that are not square; it must be greater
  than 0. ``NODATA_value`` is optional. Any other keyword stops the run.
* The values follow, ``ncols`` per row, the northern row first. There must be
  exactly ``ncols`` times ``nrows`` of them, each a finite number: ``nan``
  and ``inf`` stop the run with the position of the value. Mark missing
  data with ``NODATA_value`` instead.
* The raster must be on the coordinates of the run, in metres, and must hold
  the whole domain (``flume.prob_lo`` to ``flume.prob_hi``, within one
  millionth of a raster cell). FLUME does not reproject or change units.
  Only the raster cells that hold the domain, and one more on each side, are
  kept; the rest of the file may be much larger than the domain. The I/O
  rank reads the whole file and holds it as text while it parses it, so the
  file must fit in the memory of one rank; the part kept, at most
  2\ :sup:`31` - 1 cells, is then sent to every rank.
* The bed is interpolated bilinearly between raster cell centres. Within half
  a raster cell of the raster's edges the values of the edge cells are held.
  A plane given at the raster cell centres is reproduced exactly.
* Buildings and the Manning coefficient take the raster cell that holds the
  cell centre. A centre on a raster cell edge takes the cell above it, in x
  and in y. "On the edge" allows for round-off: a position, in raster cells,
  within 16 epsilons of the build's precision of a whole number counts as
  that number.
* On a periodic domain the raster is not wrapped: the bed is what the file
  gives on each side, so a raster that is not periodic itself leaves a step
  in the bed across the periodic seam.
* A cell that reads a ``NODATA_value`` cell (with a weight above 0, for the
  bed) stops the run with the number of such cells.

Each message about a file names the key and the file. A value beyond the
range of the build's precision (above about 3.4e38 in single precision) stops
the run with its position. Rasters through GDAL (GeoTIFF and other formats)
are kept for a later version.

Hyetograph files
----------------

``flume.rain.file`` is a text file of lines of two numbers: a time in s and
a rain rate in mm/h. For example::

   # time (s)  rate (mm/h)
   0      36
   600    72   # the rate doubles
   1200   0

* ``#`` starts a comment that runs to the end of the line. Blank lines are
  skipped. Numbers are separated by spaces or tabs.
* Each line that is not blank holds exactly two numbers, both finite; the
  rate is at least 0; the times increase strictly from one line to the next.
  The file holds at least one line. Times may be negative and need not start
  at 0.
* The rate between the times is held (``flume.rain.interpolation = step``)
  or linear (``linear``); before the first time and after the last it is
  the first or the last rate (:doc:`Model`). A storm that starts after the
  run starts needs a first line ``0 0`` (or a line at its start with rate
  0); a storm that ends needs a last line with rate 0. When the first time
  is after 0 with a rate above 0, the run prints a note that this rate falls
  from the start of the run.

One rank reads the file and sends the table to the others. A file that
cannot be opened, a path that is not a regular file, and each failed check
stop the run on every rank with one message that names the key, the file
and the line.

Single precision and large coordinates
--------------------------------------

The corners are stored in the precision of the build. In single precision a
corner at UTM-size coordinates is rounded to the nearest representable number,
which is 0.5 m apart at 4.4e6 m: ``4400123.37`` becomes ``4400123.5``. The cell
size then differs from the one requested. Use double precision, or corners near
the origin, when this matters. Raster files are read in double precision; the
position of each cell in the raster is formed from the difference between the
domain corner, as the build stores it, and the raster's corner in double
precision before it is rounded to the precision of the build, so map
coordinates do not cost digits in the sampling. In single precision the
domain corner itself is rounded as above, and the cells sit where the rounded
corner puts them.

Output
------

The run writes AMReX plotfiles named ``<prefix><step>``, with the step in at
least five digits: the initial state, every ``flume.plot_int`` steps, and the
last state. If a plotfile of that name exists, AMReX first renames it to
``<prefix><step>.old.<number>``. Plotfiles of other steps from an earlier run
in the same directory are left as they are.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

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
     - Bed elevation on the vertical datum, :math:`z_b`.
   * - ``surface_elevation``
     - m
     - Water-surface elevation on the same datum, :math:`z_b + h`. In a dry
       cell it is the bed elevation.

The plotfile is three-dimensional with one cell in z. Its z extent, 0 to 1 m,
carries no meaning.
