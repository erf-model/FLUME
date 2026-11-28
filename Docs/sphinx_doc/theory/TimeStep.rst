Time step
=========

With :math:`c = \sqrt{g h}`, the largest rate at which waves cross a cell is

.. math::

   r = \max_{h > h_{\mathrm{dry}}}
       \left(\frac{|u| + c}{\Delta x} + \frac{|v| + c}{\Delta y}\right) ,

taken over all boxes and all MPI ranks, and the time step is

.. math::

   \Delta t = \min\left(\frac{C}{r},\; \Delta t_{\max}\right) ,

with the CFL number :math:`C` (``flume.cfl``) and the longest step
:math:`\Delta t_{\max}` (``flume.max_dt``, no limit if not set). The rates of the two
directions are added because the scheme is unsplit.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - ``flume.order``
     - Largest :math:`C`
     - Reason
   * - 1
     - 1
     - Linear stability bound.
   * - 2
     - 0.5
     - A side of a cell holds up to 1.5 times the cell's depth. Above 0.5
       the outflow limit, first order where it acts, becomes part of the
       scheme.

The step is shortened to end on ``flume.stop_time`` (inside ERF, on the end
of each ERF step), and lengthened onto it when it falls short by round-off
of the summed steps, or by less than 10\ :sup:`-6` of the step, as ERF ends
its own run. A step of ``flume.max_dt`` a hair below ERF's step thus lands
on the end of the ERF step, with no extra step of 10\ :sup:`-9` s after it. The time is held in double precision in
every build. The run stops if the step is zero, not finite, or too small to
change the time.

.. _sec:RainStep:

Rain
----

Rain adds depth after the explicit stages, so the step also covers the state
at its end. With :math:`\hat r` the largest rain rate over the step, the
step is the largest :math:`\Delta t`, no larger than the CFL step, with

.. math::

   \Delta t \left(R + S \sqrt{\hat r\, \Delta t}\right) \le C ,
   \qquad
   S = \sqrt{g}\left(\frac{1}{\Delta x} + \frac{1}{\Delta y}\right) ,
   \qquad
   R = \max\left(r,\; S \sqrt{h_{\mathrm{dry}}}\right) .

The bound follows from :math:`\sqrt{g\,(h + \hat r \Delta t)} \le
\sqrt{g h} + \sqrt{g\, \hat r\, \Delta t}`. It is solved for
:math:`\sqrt{\Delta t}` by Newton's method.

An inflow inside the domain (:ref:`sec:InflowSource`) raises the depth of
its cells by at most :math:`f_{\max} \hat Q / (\Delta x\, \Delta y)` per
second, with :math:`\hat Q` its largest discharge over the step and
:math:`f_{\max}` the largest part of it that one cell takes. These rates
are added to :math:`\hat r`, so the bound holds where rain and inflows meet
in one cell, and a dry domain fed only by an inflow gets a finite step.
The rates are added even for inflows in cells apart, so the step is smaller
than it needs to be there: with :math:`N` equal inflows on a dry domain,
by about :math:`N^{1/3}`.

With rain from ERF's microphysics (``flume.rain.type = microphysics``,
:ref:`sec:InsideERF`), the precipitation :math:`D` of an ERF step from
:math:`T_0` to :math:`T_1` is a depth on each cell, held at the constant
rate :math:`D / (T_1 - T_0)`. A FLUME step from :math:`t_0` to :math:`t_1`
inside it adds

.. math::

   D \, \bigl(f(t_1) - f(t_0)\bigr) ,
   \qquad
   f(t) = \min\Bigl(1,\; \max\Bigl(0,\; \frac{t - T_0}{T_1 - T_0}\Bigr)\Bigr) ,

so that the steps of the ERF step add :math:`D` together, and :math:`\hat r`
is the largest :math:`\max(D, 0) / (T_1 - T_0)` over the cells.

The rain bound gives a dry domain under rain a finite step. On water of
0.1 m and more it costs under 2 % of the CFL step. On a domain that is all
thin film it shortens the step to about half (0.6 at 50 mm/h and 0.45 at
200 mm/h, on films of 10\ :sup:`-4` m over 5 m cells).

.. _sec:DryDomain:

Dry domain
----------

If no cell is wet, :math:`r = 0`. The fluxes are not computed (the run
prints so, once for each such spell) and the sources still act.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Case
     - Step
   * - Rain ahead
     - The rain bound sets the step.
   * - No rain, ``flume.max_dt`` given
     - Steps of ``flume.max_dt``, the last one shortened onto the stop
       time if there is one.
   * - No rain, ``flume.stop_time`` given, no ``flume.max_dt``
     - One step to the stop time.
   * - No rain, ``flume.max_step`` alone
     - The run stops with a message asking for ``flume.stop_time`` or
       ``flume.max_dt``.
