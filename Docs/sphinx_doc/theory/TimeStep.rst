Time step
=========

With :math:`c = \sqrt{g h}`, the largest rate at which waves cross a cell is

.. math::

   r = \max_{h > h_{\mathrm{dry}}}
       \left(\frac{|u| + c}{\Delta x} + \frac{|v| + c}{\Delta y}\right) ,

taken over all boxes and all MPI ranks, and the time step is

.. math::

   \Delta t = \frac{C}{r} ,

with the CFL number :math:`C` (``flume.cfl``). The rates of the two
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

The step is shortened to end on ``flume.stop_time``, and lengthened when it
falls short of it by round-off only. The time is held in double precision in
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
   * - No rain, ``flume.stop_time`` given
     - One step to the stop time.
   * - No rain, ``flume.max_step`` alone
     - The run stops with a message asking for ``flume.stop_time``.
