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

The run writes ``dam00000`` (initial state), the plotfile of the last step
and the water-balance log ``water_balance.txt`` (:ref:`sec:WaterBalanceLog`)
in the working directory, and ends with a line such as
``FLUME: water balance (m3): stored 150 (at step 0 150), residual ...``.
:doc:`Verification` compares this case with its exact solution.

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
   * - ``max_depth``, ``max_speed``, ``max_depth_speed``, ``arrival_time``
     - m, m/s, m\ :sup:`2`/s, s
     - The flood maps, written only with ``flume.flood_maps = 1``
       (:ref:`sec:FloodMaps`).

The plotfile is three-dimensional with one cell in z, from 0 to 1 m; the z
extent carries no meaning. ParaView, VisIt and yt read the format.

The log reports, at start-up, the boundary cells of each segment, the active
source terms, the rain, the infiltration parameters, the range of the Manning map and the
number of building cells.

.. _sec:FloodMaps:

Flood maps
~~~~~~~~~~

With ``flume.flood_maps = 1`` every plotfile also holds, for each cell, what
the run has reached up to that plotfile:

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Variable
     - Unit
     - Meaning
   * - ``max_depth``
     - m
     - Largest depth.
   * - ``max_speed``
     - m/s
     - Largest speed, :math:`\sqrt{u^2 + v^2}`; 0 while the cell is dry.
   * - ``max_depth_speed``
     - m\ :sup:`2`/s
     - Largest depth times speed, the magnitude of the discharge
       :math:`\sqrt{(hu)^2 + (hv)^2}`; 0 while the cell is dry.
   * - ``arrival_time``
     - s
     - First time the depth was above ``flume.flood_maps.arrival_depth``;
       -1 where it has not been.

A cell is dry when its depth is at or below ``flume.dry_depth``. The maps
follow the state at the start and at the end of every step, not the stages
inside a step, so an arrival time is the end of a step, up to one step late.
A cell wet at the start has arrival time 0 when its depth is above the
arrival depth, which is at least ``flume.dry_depth``. The maximum of depth
times speed is taken by itself: it need not come at the time of the largest
depth or speed. The maps are held in the precision of the state: in a
single-precision build an arrival time near :math:`10^6` s is rounded to
about 0.06 s.

.. _sec:WaterBalanceLog:

Water-balance log
~~~~~~~~~~~~~~~~~

Unless ``flume.log_int = 0``, the I/O rank writes ``water_balance.txt`` in
the working directory: one line at the start, one every ``flume.log_int``
steps, and one at the last step. An existing file of that name is
overwritten (plotfiles are renamed instead), except on a restart
(:ref:`sec:Restart`). At the end the run prints the
last line as ``FLUME: water balance (m3): ...``. Each line holds, in m\ :sup:`3` and summed from the start of the run,
the water on the surface, the water each source has added or removed, the
water that has left through each outflow segment, and the residual of the
balance (:doc:`theory/WaterBalance`). The header names the columns:

.. code-block:: text

   # FLUME water balance: volumes in m3, cumulative from step 0
   # epsilon 2.2204460492503131e-16 (machine epsilon of the double-precision state)
   # residual = stored - stored at step 0 - (added) + (removed) + (out)
   # added: rain
   # removed: infiltrated
   # out: out_out
   # step time stored rain infiltrated out_out residual
   0 0 0 0 0 0 0

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Column
     - Meaning
   * - ``step``, ``time``
     - Step number and time, s.
   * - ``stored``
     - Water on the surface, :math:`V`.
   * - ``rain``
     - Rain added to the domain (with ``flume.rain.type`` other than
       ``none``).
   * - ``infiltrated``
     - Water lost to the soil, as a positive number (with
       ``flume.infiltration.type`` other than ``none``).
   * - ``out_<segment>``
     - Water that has left through outflow segment ``<segment>``; below 0
       where more has come in than gone out.
   * - ``residual``
     - Stored water minus that at step 0, minus the water added, plus the
       water removed and out; zero in exact arithmetic.

The ``# added:``, ``# removed:`` and ``# out:`` lines list the columns that
enter the residual with each sign. ``# epsilon`` is the machine epsilon of the build.
A domain without outflow segments, rain or infiltration has only ``stored``
and ``residual``. ``Tests/check_balance.py`` checks a log against the bound
of :doc:`theory/WaterBalance`; the regression tests run it on every log.

.. _sec:Gauges:

Gauges
~~~~~~

Each gauge of ``flume.gauges.names`` writes ``gauge_<name>.txt`` in the
working directory (overwriting a file of that name, except on a restart) at
the lines of the
water-balance log:

.. code-block:: text

   # FLUME gauge face at x = 24 m, y = 7 m, cell (12, 3)
   # step time depth surface_elevation velocity_x velocity_y
   0 0 0 5 0 0

Depth and surface elevation are in m, the velocities in m/s. A gauge reads
the cell that holds its point; there is no interpolation. A point on a face
between two cells, within 4 machine epsilons of the largest coordinate
(and at most :math:`10^{-3}` of a cell), reads the cell above it in that
direction, and a point on the high edge of the domain
the last cell. The start-up log names the cell of each gauge. The velocity
is 0 in a dry cell. The values are printed with the digits that give back
the stored number.

.. _sec:Restart:

Checkpoint and restart
~~~~~~~~~~~~~~~~~~~~~~

With ``flume.check_int > 0`` the run writes a checkpoint directory
``<flume.check_file><step>`` (the step in at least five digits) every
``flume.check_int`` steps and at the last step. An existing directory of
that name is first renamed, as plotfiles are. A checkpoint holds the step
and the time, the state and the bed, the infiltrated depth (and, for the
curve number, the rain fallen on each cell), the flood maps, and the
water-balance totals with their compensation terms. The fields are AMReX
``VisMF`` files of the valid cells; ``Header`` is a text file of keys and
values, with times and volumes written to 17 digits.

``flume.restart=<dir>`` continues the run from that checkpoint:

.. code-block:: bash

   mpiexec -n 2 ./flume inputs flume.check_int=400 flume.max_step=550
   mpiexec -n 1 ./flume inputs flume.check_int=400 flume.restart=chk00400

The restart reads the same inputs. ``flume.stop_time`` and
``flume.max_step`` count from the start of the run, not from the
checkpoint. A restart with a later ``flume.stop_time`` than the first part
continues the run, but the last step of the first part was shortened to end
at its ``flume.stop_time``, so the result differs from that of one run to
the later time.

The restart writes no plotfile, log line or checkpoint at the step it
starts from. ``water_balance.txt`` and the gauge files are cut back in
place: the lines of later steps are removed, and so is a line cut short by
a run that was stopped while writing it (a line with fewer fields than the
header names). The line of the checkpoint's step is kept when that step is
a multiple of ``flume.log_int``, as in the run without a break. A missing
file is started again with its header, from the checkpoint's step;
``Tests/check_balance.py``, which needs the line of step 0, does not accept
such a file. Plotfiles and checkpoints of later steps that the first part
wrote are not removed; the restart overwrites those of the steps it reaches
(the old ones are renamed), and a plotfile the first part wrote at its last
step, off the plot interval, stays. ``flume.check_file`` must differ from
``flume.plot_file``.

The run stops with a message naming ``flume.restart`` when the checkpoint
does not match the inputs: another grid, domain or set of periodic sides,
another list of sources, a bed that differs in any cell, other water-balance
columns (rain, infiltration and outflow segments), outflow segments on other
faces, gauges in other cells, a build of the other precision, a checkpoint
written with ``flume.log_int = 0`` when the restart keeps a log, or one
without maps when ``flume.flood_maps = 1``. It also stops on a checkpoint of
another version, one with a key given twice, and a state that is not finite
in some cell. Other inputs, such as the CFL number, the Manning coefficient
or ``flume.bed.continue_slope`` (which sets the bed outside the domain), are
not checked: a change there changes the result from the checkpoint on.

The restart may use another rank count and other boxes. The state, the
infiltrated depth, the maps, the plotfiles and the gauge files are then
those of the run without a break, to the last bit. The water-balance totals
of the outflow and the rain are sums over the ranks: on another rank count
they are continued from their sum over the ranks, and the log agrees with
the run without a break to round-off (below :math:`10^{-14}` of the largest
value of a column in the tests), not to the bit. On the same rank count and
boxes the log is the same byte for byte.
