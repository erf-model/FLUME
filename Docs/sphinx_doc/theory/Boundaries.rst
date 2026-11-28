Boundary conditions
===================

Each direction is periodic or not (``flume.periodic``). The state and the
bed carry one layer of ghost cells at first order and two at second order.

A side that is not periodic is made of segments (``flume.bc.segments``),
each a run of boundary cells with one type. A boundary cell belongs to a
segment when its centre is at the first coordinate of
``flume.bc.<name>.range`` or above, and below the second. Cells no segment
covers are walls.

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Type
     - Ghost cells
   * - ``wall``
     - Mirror the cells inside: same depth, bed and tangential discharge;
       normal discharge with its sign changed. The mass flux through a wall
       is exactly zero.
   * - ``outflow``
     - Take the boundary cell next to them, the bed included (zero
       gradient), where the discharge of that cell is zero or points out of
       the domain; mirror the cells inside, as a wall, where it points in.
       Water leaves and does not enter (:ref:`sec:Outflow`).
   * - ``zero_gradient``
     - Take the boundary cell next to them, the bed included, whatever the
       direction of the flow. Water leaves or enters, and nothing bounds
       what enters (:ref:`sec:ZeroGradient`).
   * - ``discharge``
     - Take the boundary cell next to them (zero gradient), as for
       ``zero_gradient``; the fluxes of the boundary faces are then replaced
       by those of the inflow hydrograph (:ref:`sec:Inflow`).
   * - periodic
     - Take the cell they cover.

``stage`` is reserved for a prescribed water level.

.. _sec:Outflow:

Outflow
-------

An ``outflow`` segment lets water leave and lets none in. Let
:math:`q_n` (m\ :sup:`2`/s) be the discharge of a boundary cell along the
normal that points into the domain, :math:`h` its depth, :math:`q_t` its
discharge along the side, and let the subscript :math:`m` mark the
:math:`m`-th cell inside the side. Ghost layer :math:`m` of that boundary
cell takes

.. math::

   (h, q_n, q_t)_{\text{ghost}} =
   \begin{cases}
   (h_1, \; q_{n,1}, \; q_{t,1}) & q_{n,1} \le 0 \quad \text{(zero gradient)} \\
   (h_m, \; -q_{n,m}, \; q_{t,m}) & q_{n,1} > 0 \quad \text{(wall)}
   \end{cases}

and the bed of the boundary cell in both cases.

**Why.** With the copy the two states of the boundary face are equal, so
the face passes the discharge of the boundary cell, which points out of the
domain or is zero. With the mirror the mass flux of the face is zero, as
through a wall. So the mass flux of an outflow face never points into the
domain, and the ``out_<segment>`` column of the water balance never falls
(to the rounding of one product per step, as for a wall).

**What it costs.** A cell whose flow turns towards the domain sees a wall
until it turns back: a wave that reaches the side while the flow there
points inward is reflected. The choice is made cell by cell and stage by
stage from the state alone, so it needs no memory and a restart repeats it.
Water from outside, such as a tide or a surge, cannot enter through an
``outflow`` segment; a prescribed level (``stage``) is reserved and not
implemented.

Where the flow leaving is supercritical the side sets nothing the flow
inside depends on. Where it is subcritical, zero gradient is an
approximation, with these limitations:

* A shock that leaves sends a weak wave back. In the dam-break test it
  changes the depth by 0.2 % and the velocity by 0.3 % at first order, and
  by 0.7 % and 2 % at second order. It does not shrink with the cell size.
* Still water next to an outflow side is not stable where the bed of the
  boundary cell is lower than that of the cell inside it: a disturbance
  grows from round-off over thousands of steps (up to 9 % per step over a
  rough bed in a 1-D model of the scheme). The disturbance can only drain
  the water: the pond empties through the side (197.8 of 200.8
  m\ :sup:`3` in an hour at first order, 197.5 at second, in the case of
  ``Tests/Outflow/inputs_step``, whose water starts at 10\ :sup:`-4` m/s;
  with no initial velocity the same pond loses less than
  10\ :sup:`-9` m\ :sup:`3` in 8 h). Water that should stay ponded against a side needs a wall
  there.

* On a sloping bed, subcritical flow toward an outflow side backs up. The
  ghost cells copy the bed, so the boundary face passes the discharge of the
  cell, while the faces inside carry the bed-slope term of the hydrostatic
  reconstruction as well. In a channel of slope 0.02 with
  :math:`n` = 0.1 s/m\ :sup:`1/3`, fed with 1 m\ :sup:`2`/s (normal depth
  0.812 m, Froude number 0.44) in cells of 0.5 m, the depth at the outflow
  end at 1500 s is 1.28 m at first order and 1.83 m at second, and the
  backwater reaches 55 m and 86 m up the channel. In a steep (supercritical)
  channel only the last few cells are bent.

.. _sec:ZeroGradient:

Zero gradient
-------------

A ``zero_gradient`` segment copies the boundary cell into its ghost cells
whatever the direction of the flow. It was the ``outflow`` type before the
one above, and gives the same result wherever no boundary cell has a
discharge that points into the domain. Where one has, water enters, and
nothing bounds how much:

* **A bed that falls into the domain.** The ghost cells hold the depth and
  the discharge of the boundary cell, as if the slope went on outside. Rain
  of 50 mm/h on a plane of 200 m by 200 m (cells of 2 m, :math:`n` = 0.015
  s/m\ :sup:`1/3`, 0.05 m of bed noise) that falls by 1 cm per metre in
  :math:`x` and 0.5 cm per metre in :math:`y`, open on four sides, draws
  0.13 m\ :sup:`3`/s through its upper side in :math:`x` once the flow is
  steady, 23 % of the 0.56 m\ :sup:`3`/s of rain.
* **A bed that steps down to the side under ponded water.** Let the bed of
  the boundary cell lie :math:`d` below that of the cell inside, under a
  level surface, with the water moving into the domain at the speed
  :math:`u`. The boundary face passes :math:`h u`, the discharge of the
  whole depth. The face inside passes :math:`(h - d) u`: the hydrostatic
  reconstruction gives it the depth above the higher bed. The boundary cell
  gains :math:`d u` per unit width and time, and

  .. math::

     \Delta x \, \frac{\partial (h u)}{\partial t} = d \, u^2 > 0 ,

  so its momentum is not opposed. Its depth rises, the ghost cells rise
  with it, and the inflow grows: no level outside holds it back. Still
  water 0.5 m deep over a step of 0.1 m, moving inward at
  10\ :sup:`-4` m/s, grows from 200.8 m\ :sup:`3` to 8974 m\ :sup:`3` in
  300 s (``Tests/Outflow/inputs_step``). Rain ponding against a wall
  between two ``zero_gradient`` sides takes in water through one side and
  passes it out of the other: after an hour the domain stores
  1.1 :math:`\times` 10\ :sup:`7` m\ :sup:`3` from 1000 m\ :sup:`3` of
  rain (``Tests/Outflow/inputs_pond``, run to 3600 s). The residual of the
  water balance stays at round-off, since the water that enters is counted
  as negative outflow; a falling ``out_<segment>`` column is the sign.

Use ``zero_gradient`` where the flow is known to leave over the whole
segment for the whole run, or to compare with ``outflow``.

.. _sec:Inflow:

Inflow hydrograph
-----------------

A ``discharge`` segment takes the discharge :math:`Q(t)` (m\ :sup:`3`/s)
of a file (``flume.bc.<name>.file``) and spreads it evenly over the width of
the segment. The faces of the segment carry fluxes set from it; the rest of
the scheme sees them like any other face fluxes.

.. math::

   q = \frac{Q}{W} , \qquad
   F_h = q , \qquad
   F_{n} = q \, u_b + \tfrac{1}{2} g h_b^2 , \qquad
   F_{t} = 0 ,

with :math:`W` the width of the segment (m), :math:`q` the discharge per
unit width (m\ :sup:`2`/s), :math:`F_h` the mass flux into the domain,
:math:`F_n` the flux of the momentum along the normal into the domain,
:math:`F_t` that of the momentum along the side, :math:`g` the gravity,
:math:`h_b` the depth at the face (m) and :math:`u_b = q / h_b` the speed at
the face (m/s).

**The depth at the face.** The cell inside the face has the depth :math:`h`
and the velocity :math:`u_n` along the normal into the domain; a cell at or
below ``flume.dry_depth`` counts as dry, :math:`h = u_n = 0`. Where the
inflow is subcritical one characteristic, of speed :math:`u_n - c`, leaves
the domain and brings the invariant :math:`R = u_n - 2c`,
:math:`c = \sqrt{g h}`, to the face. The face state keeps it:

.. math::

   \frac{q}{h_b} - 2 \sqrt{g h_b} = R .

The left side falls from :math:`+\infty` to :math:`-\infty` as
:math:`h_b` grows, so the root is unique. It is subcritical exactly when it
lies above the critical depth

.. math::

   h_c = \left( \frac{q}{\sqrt{g}} \right)^{2/3} ,

that is when :math:`R < -\sqrt{g h_c}`. Otherwise both characteristics
enter the domain and the face state comes from outside alone; FLUME then
takes :math:`h_b = h_c`, the depth at which a given discharge carries the
least energy. That covers inflow into a dry cell. The root is found by
Newton's method in :math:`s = \sqrt{h_b}` from :math:`\sqrt{h_c}`. With
:math:`q = 0` the face keeps :math:`R` with :math:`u_b = 0`, so
:math:`h_b = \max(-R/2, 0)^2 / g`, which is :math:`h` for still water: a
lake at rest next to a segment with no discharge stays at rest.

**Why the mass is exact.** Over a step from :math:`t_0` to :math:`t_1` both
stages of the update take the same flux,

.. math::

   q = \frac{I(t_1) - I(t_0)}{\Delta t \, W} , \qquad
   I(t) = \int_0^t Q(\tau) \, d\tau ,

with :math:`I` the integral of the hydrograph as interpolated (``step`` or
``linear``, ``flume.bc.<name>.interpolation``) and :math:`\Delta t` the
step. The water that enters over the run is therefore the integral of the
hydrograph, to round-off, whatever the steps. The water balance counts it in
the column ``in_<name>`` (:ref:`sec:WaterBalanceLog`).

**The time step.** The step is bounded by the wave rate of the face states
as well as by that of the cells: :math:`(u_b + c_b) / \Delta n +
c_b / \Delta s`, with :math:`c_b = \sqrt{g h_b}`, :math:`\Delta n` the
cell size along the normal and :math:`\Delta s` along the side, for the
largest discharge of the hydrograph from the current time to the end of the
longest step the cells and the stop time allow. Over a dry domain without a
stop time that is the rest of the hydrograph: a later peak then sets the
step from the start, which is safe but smaller than needed. Water entering
a dry domain therefore sets the step from the first one; a domain that is dry and
gets no rain and no more inflow still needs ``flume.stop_time``.

**Limitations.** The discharge is uniform across the segment; a segment
that should carry more water in its middle is split into several segments.
The water enters along the normal. Flow inside that moves toward the
segment (:math:`u_n < 0`) is met with the depth of the invariant, as if the
two states were joined by a rarefaction; the exact connection is a shock,
so the face depth is too high: by 0.7 % at a Froude number of 0.3 and 10 %
at 1.6 for 1 m of water, about as much as at a wall. The critical depth for supercritical
inflow is a choice: a hydrograph does not fix the depth at which water
enters a steep channel, and the flow a short distance inside is set by the
slope and the friction. Because of the limitation of the outflow side
above, a subcritical channel on a slope is tested here only with a
supercritical (steep) slope.

Periodic sides and a sloping bed
--------------------------------

By default the bed repeats across a periodic side, with a step of the slope
times the length of the domain. With ``flume.bed.continue_slope = 1`` the
ghost cells take the bed of the cell they cover plus the slope times the
distance between the two: an endless inclined channel. The depth and the
discharges still repeat. Only the plane is continued, not a bowl or a
raster.
