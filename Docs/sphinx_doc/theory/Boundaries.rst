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
       gradient). Nothing is prescribed.
   * - periodic
     - Take the cell they cover.

``discharge`` and ``stage`` are reserved for a prescribed inflow and a
prescribed water level.

Outflow
-------

Where the flow leaving is supercritical the side sets nothing the flow
inside depends on. Where it is subcritical, zero gradient is an
approximation, with two limitations:

* A shock that leaves sends a weak wave back. In the dam-break test it
  changes the depth by 0.2 % and the velocity by 0.3 % at first order, and
  by 0.7 % and 2 % at second order. It does not shrink with the cell size.
* Still water next to an outflow side is not stable where the bed of the
  boundary cell is lower than that of the cell inside it: a disturbance
  grows from round-off over thousands of steps (up to 9 % per step over a
  rough bed in a 1-D model of the scheme). Runs whose water comes to rest
  against an outflow side over a bed that steps down to the side should not
  stay long in that state.

Periodic sides and a sloping bed
--------------------------------

By default the bed repeats across a periodic side, with a step of the slope
times the length of the domain. With ``flume.bed.continue_slope = 1`` the
ghost cells take the bed of the cell they cover plus the slope times the
distance between the two: an endless inclined channel. The depth and the
discharges still repeat. Only the plane is continued, not a bowl or a
raster.
