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
   * - ``drained_depth``
     - m
     - Depth of water each cell has lost to a storm-drain inlet since the
       start (:ref:`sec:Drains`); 0 in a cell without an inlet. Written only
       with ``flume.drains.names``.
   * - ``max_depth``, ``max_speed``, ``max_depth_speed``, ``arrival_time``
     - m, m/s, m\ :sup:`2`/s, s
     - The flood maps, written only with ``flume.flood_maps = 1``
       (:ref:`sec:FloodMaps`).

The plotfile is three-dimensional with one cell in z, from 0 to 1 m; the z
extent carries no meaning. ParaView, VisIt and yt read the format.

The log reports, at start-up, the boundary cells of each segment, the active
source terms, the rain, the infiltration parameters, each drain inlet with its cell and
the depth above which it flows as an orifice, the range of the Manning map and the
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
A cell wet at the start has the start time as its arrival time (0, or
ERF's time when FLUME starts inside ERF from a checkpoint without flood
state) when its depth is above the arrival depth, which is at least
``flume.dry_depth``. The maximum of depth times speed is taken by itself: it need not come at the time of the largest
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
water that has left through each outflow segment, the water that has
entered through each discharge segment, and the residual of the
balance (:doc:`theory/WaterBalance`). The header names the columns:

.. code-block:: text

   # FLUME water balance: volumes in m3, cumulative from step 0
   # epsilon 2.2204460492503131e-16 (machine epsilon of the double-precision state)
   # residual = stored - stored at step 0 - (added) + (removed) + (out) - (in)
   # added: rain
   # removed: infiltrated
   # out: out_out
   # in: in_inlet
   # step time stored rain infiltrated out_out in_inlet residual
   0 0 0 0 0 0 0 0

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
   * - ``drained``
     - Water taken by the storm-drain inlets, as a positive number (with
       ``flume.drains.names``).
   * - ``out_<segment>``
     - Water that has left through outflow segment ``<segment>``; below 0
       where more has come in than gone out.
   * - ``in_<segment>``
     - Water that has entered through discharge segment ``<segment>``: the
       integral of its hydrograph from the start, to round-off.
   * - ``residual``
     - Stored water minus that at step 0, minus the water added and in,
       plus the water removed and out; zero in exact arithmetic.

The ``# added:``, ``# removed:``, ``# out:`` and ``# in:`` lines list the columns that
enter the residual with each sign; ``# in:`` (and the ``- (in)`` of the
residual line) appear only with a discharge segment. ``# epsilon`` is the machine epsilon of the build.
A domain without outflow or discharge segments, rain, infiltration or drain inlets has only ``stored``
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
curve number, the rain fallen on each cell), the drained depth with the inlets it belongs to, the flood maps, and the
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
columns (rain, infiltration, drain inlets, outflow and discharge segments),
outflow or discharge segments on other faces, gauges in other cells, drain inlets with another
name, cell, perimeter, area or coefficient, a build of the other precision, a checkpoint
written with ``flume.log_int = 0`` when the restart keeps a log, or one
without maps when ``flume.flood_maps = 1``. It also stops on a checkpoint of
another version, one with a key given twice, and a state that is not finite
in some cell. Other inputs, such as the CFL number, the Manning coefficient,
``flume.bed.continue_slope`` (which sets the bed outside the domain), or the
hydrograph file and interpolation of a discharge segment, are not checked: a
change there changes the result from the checkpoint on.

The restart may use another rank count and other boxes. The state, the
infiltrated depth, the drained depth, the maps, the plotfiles and the gauge files are then
those of the run without a break, to the last bit. The water-balance totals
of the outflow, the inflow and the rain are sums over the ranks: on another rank count
they are continued from their sum over the ranks, and the log agrees with
the run without a break to round-off (below :math:`10^{-14}` of the largest
value of a column in the tests), not to the bit. On the same rank count and
boxes the log is the same byte for byte.

.. _sec:InsideERF:

Inside ERF
----------

With ``-DFLUME_WITH_ERF=ON`` (:doc:`Building`) FLUME is compiled into
``erf_exec`` as a flood engine. ERF runs it with
``erf.flood.enable = true`` and ``erf.flood.engine = flume``:

.. code-block:: bash

   mpiexec -n 2 build-erf/ERF/Exec/erf_exec Tests/ERF/inputs_erf \
       FILE=Tests/ERF/flume_plane

ERF builds the flood grid. It covers the region of the anchor level
(``erf.flood.anchor_level``): the whole domain on level 0, or a finer level's
refined region, such as a catchment. That region is refined in x and y by
``erf.flood.grid_ratio``. The flood boxes are the anchor level's boxes,
refined, and each sits on the rank of the atmospheric columns above it.

FLUME reads its own ``flume.*`` keys as in a stand-alone run, except the keys
that ERF sets. Each of these stops the run with a message that names ERF's
key instead:

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - FLUME key
     - ERF sets it with
   * - ``flume.n_cell``, ``flume.prob_lo``, ``flume.prob_hi``,
       ``flume.periodic``, ``flume.max_grid_size``
     - ``geometry.prob_lo``, ``geometry.prob_hi`` (or
       ``geometry.prob_extent``), ``geometry.is_periodic``,
       ``amr.n_cell``, ``amr.max_grid_size``, ``erf.flood.grid_ratio``,
       ``erf.flood.anchor_level``
   * - ``flume.plot_file``, ``flume.plot_int``
     - ``erf.flood.plot_file``, ``erf.flood.plot_int``,
       ``erf.flood.plot_per``
   * - ``flume.check_file``, ``flume.check_int``
     - ``erf.check_file``, ``erf.check_int``, ``erf.check_per``
   * - ``flume.restart``
     - ``erf.restart``
   * - ``flume.stop_time``, ``flume.max_step``
     - ``stop_time``, ``max_step``

FLUME starts at ERF's start time, and its times (rain tables, hydrographs)
are on ERF's clock. In each ERF step of the anchor level it takes its own
steps, each set by its CFL number, the rain and ``flume.max_dt``, and lands
on the end of the ERF step. The rain is FLUME's own (``flume.rain.*``), or,
with ``flume.rain.type = microphysics``, ERF's:

- In each ERF step, ERF hands over the surface precipitation of each
  atmospheric column over the step, kg/m\ :sup:`2`: the liquid part (rain)
  and the frozen part (snow, graupel, hail). FLUME takes both, as liquid
  water: 1 kg/m\ :sup:`2` is a depth of 1 mm on every flood cell under the
  column. There is no snowpack: the frozen part becomes water on the flood
  grid at once, so the timing of snowmelt is not modelled.
- ERF's precipitation of a step can be negative in a column (a correction
  from a finer level, or rounding). FLUME does not lower the depth for it:
  the amount waits in the cell's carry and is taken from the next rain. The
  ``rain`` column of the water balance leaves it out until then.
- FLUME holds the rate constant over the ERF step. Each of its steps inside
  the ERF step gets the share of its length, so the steps of one ERF step
  together get the ERF step's precipitation. The rain-limited step
  (:ref:`sec:RainStep`) uses that rate.
- ERF stops the run at the start if its microphysics accumulates no surface
  precipitation (``erf.moisture_model`` such as ``Kessler``, ``SAM``,
  ``Morrison``, ``WSM6`` or ``WDM6`` does). Every other ``flume.rain.*``
  key (``rate``, ``file``, ``interpolation``, ``cell.*``) stops the run
  with ``microphysics``.
- ERF's checkpoint keeps its reference of the accumulated precipitation, so
  a restart continues with the same rain. From a checkpoint written without
  the flood, the reference is the accumulation ERF restores with its
  microphysics, so the rain also continues from there.

The bed is FLUME's own (``flume.bed.*``),
not ERF's terrain, so ``flood_bed_elevation`` can differ from it.

ERF writes the flood plotfiles: FLUME's plotfile (:ref:`sec:Output`), named
``<erf.flood.plot_file><ERF step>`` (``plt_flood_`` by default), with
FLUME's time and step in the header. Without ``erf.flood.plot_int`` and
``erf.flood.plot_per`` they follow ERF's ``erf.plot_int_1`` and
``erf.plot_per_1``. There is one at step 0 and one at each ERF plot interval; there is
none at the last step unless it falls on the interval. The ERF plotfile can
also carry ``flood_depth``, ``flood_discharge_x``, ``flood_discharge_y``,
``flood_bed_elevation`` and ``flood_surface_elevation``, averaged over each
atmospheric column. FLUME writes ``water_balance.txt`` and the gauge files
in the working directory, as in a stand-alone run, and the water-balance
summary at ERF's last step: the step ``max_step`` sets, or the one that
ends within 10\ :sup:`-6` of a step of ``stop_time``, as ERF stops.

ERF's checkpoints carry FLUME's, in ``<checkpoint>/Level_<lev>/FLUME``.
``erf.restart`` continues FLUME from there, and cuts the logs back as a
stand-alone restart does (:ref:`sec:Restart`), on the same or another
number of ranks. A checkpoint written without the flood
(``erf.flood.enable = false``) holds no FLUME state: FLUME then starts from
its initial state at ERF's time. Its logs and maps start at that time (the
arrival time of a cell wet at the start is ERF's time), and it counts its
steps on from ERF's, which it reads from the checkpoint's ``Header``, so
its summary is written at ``max_step`` as in a run from step 0.

FLUME's steps inside an ERF step are its own: each is the CFL, rain and
``flume.max_dt`` step, except the one that would pass the end of the ERF
step, which is shortened to end on it. An ERF step shorter than FLUME's own
step is one FLUME step. So no FLUME step is longer than those limits allow
(the step that lands can pass them by under 10\ :sup:`-6` of itself, as in a
stand-alone run on ``flume.stop_time``). In ``erf_flume_plane_step_lengths``
(ERF steps from 1.5 to 15 s) no run inside ERF is further from the
time-converged flood than the stand-alone run with FLUME's own steps, within
2 %, and each differs from that run by at most half of its time error. ERF
does not shorten its step for its plotfiles or checkpoints, so the flood
does not depend on where they fall: a run with a checkpoint at an ERF step
after several FLUME substeps, and the run continued from it, give the flood
of the run without one, to the last bit
(``erf_flume_plane_substeps_restart``). With a variable ERF step
(``erf.cfl``) the flood follows ERF's steps and is not the same to the last
bit as with another sequence of steps.

With ``erf.fixed_dt`` equal to ``flume.max_dt``, and a flood whose CFL and
rain steps are not shorter, each ERF step is one FLUME step. The flood is
then that of the stand-alone run with the same ``flume.max_dt``, to the last
bit: the plotfiles, the water balance and the gauge files. The test
``erf_flume_plane`` checks this (``Tests/ERF``). The anchor level's step
(``erf.fixed_dt`` over the substeps of the levels down to it) should be a
number that sums exactly in binary, such as 0.25 s, and the stand-alone
``flume.stop_time`` a whole number of steps. With a step such as 0.1 s the
two runs differ by round-off.

.. _sec:RainReplay:

Replay of ERF's rain
~~~~~~~~~~~~~~~~~~~~

With ``flume.rain.replay_file = <dir>/<stem>.txt`` and
``flume.rain.type = microphysics``, FLUME writes the rain it takes from ERF
as gridded rain (:ref:`sec:RainIndexFiles`). A stand-alone run of the same
flood with ``flume.rain.type = gridded`` and ``flume.rain.file`` set to that
index reads it back without the atmosphere. It gives the coupled flood to
round-off, and it can run on through the hours of recession after the rain
without advancing ERF.

- Each ERF step with rain from :math:`t_0` to :math:`t_1` adds the line
  ``t0 <stem>_<n>.asc`` to the index, where ``<n>`` counts the grids from
  ``000001``. The grid, in ``<dir>``, holds the rate of each flood cell in
  mm/h, :math:`D / (t_1 - t_0) \times 3.6 \times 10^6` for a depth
  :math:`D` in m, with 17 significant digits.
- The index always ends with ``t <stem>_zero.asc``, a grid of no rain at
  the end of the last step with rain (at the start, the start time). Read
  with ``flume.rain.interpolation = step`` (the default), the replay gets
  ERF's rain step by step and no rain after it. An ERF step of no rain adds
  no line.
- The index is written in place: a step with rain writes its line over the
  closing line before it and a new closing line after it. A run that stops
  between two steps leaves an index that ends with its closing line. At
  the start of a run the directory is made and an existing index is
  replaced (grids of an earlier, longer run stay in the directory, but the
  index does not name them).
- ERF's checkpoint holds the state of the index: its size, its closing line
  and the number of grids of rain in it. A restart cuts the index back to
  that state and goes on with the same names. The restart stops if the
  index is missing or does not hold what the run had written when it wrote
  the checkpoint (for example, the index of a run that went on from an older
  checkpoint and was stopped), or if the checkpoint was written without
  ``flume.rain.replay_file``.
- A negative depth of a step (a correction from a finer level) is not
  written in its step: it is taken from the rain the cell gets in the steps
  that follow, and the part still held goes into the checkpoint. The
  coupled run takes such a deficit from the depth at the next step with
  rain; a rate cannot be negative, so the replay takes it from that rain
  only. The replay is deeper by up to the deficit until the rain has made it
  up; the totals agree.
- The stand-alone run sets the flood grid that ERF sets in the coupled run
  (``flume.n_cell``, ``flume.prob_lo``, ``flume.prob_hi``,
  ``flume.periodic``; ``Tests/ERFReplay/inputs_replay`` is an example),
  keeps ``flume.rain.interpolation = step``, and does not set
  ``flume.rain.replay_file``. The times of the index are ERF's: a flood
  that starts with ERF at a time :math:`T > 0` has an index that starts at
  :math:`T`, and a replay from 0 s is dry until then. To end the replay
  where the coupled run ended, set ``flume.stop_time`` to the time of the
  last line of the index, written with 17 digits.
- The replay runs on any boxes and ranks. With the same steps as the
  coupled run (``flume.max_dt`` equal to ERF's step, as in
  ``Tests/ERFReplay``), it gives the coupled depths to round-off: each
  step's rain comes back within 2 units of round-off, and the flow can carry
  that on, most in thin layers of water. With FLUME steps inside each ERF
  step, or stand-alone steps that do not land on the times of the index, the
  two runs agree to within the difference in their steps.
- The grids are ESRI ASCII text: about 25 bytes per flood cell per step
  with rain. Each step with rain gathers the whole flood grid on one rank
  to write it, and the stand-alone run reads every grid at its start; the
  grid numbers have six digits, so names sort in order up to 999999 grids.
