Model
=====

This page describes what the code holds and computes today. The flux kernel,
the boundary conditions and the source terms are not written yet, so a time
step advances the time and leaves the state as it is.

State
-----

FLUME works on one level of cells in the horizontal plane. Each cell holds
three conserved quantities, all at the cell centre:

.. list-table::
   :header-rows: 1
   :widths: 16 16 68

   * - Symbol
     - Unit
     - Meaning
   * - :math:`h`
     - m
     - Water depth. Zero in a dry cell.
   * - :math:`hu`
     - m\ :sup:`2`/s
     - Discharge per unit width in x.
   * - :math:`hv`
     - m\ :sup:`2`/s
     - Discharge per unit width in y.

The bed elevation :math:`z_b` is a fourth field on the same cells. The state
and the bed are AMReX ``MultiFab`` objects on one ``BoxArray`` with its own
``DistributionMapping``. The boxes are one cell thick.

Vertical datum
--------------

The bed elevation :math:`z_b` and the water-surface elevation :math:`\eta`
are heights above one vertical datum, in metres. The inputs give the bed and
the initial surface as elevations on that datum. The depth is never an input.
It follows from the two elevations when the state is built (next section), and
afterwards the surface follows from the bed and the depth:

.. math::

   \eta = z_b + h .

In a dry cell :math:`\eta = z_b`. The plotfile holds :math:`h`, :math:`hu`,
:math:`hv`, :math:`z_b` and :math:`\eta`.

In single precision an elevation :math:`z` is held to about
:math:`6 \times 10^{-8}\,|z|`, which is 0.24 mm at 4000 m: a depth taken
as the difference of two elevations is no finer than that.

Bed and initial state
---------------------

The bed is a plane through the elevation :math:`z_0`
(``flume.bed.elevation``) at the lower corner :math:`(x_0, y_0)` of the domain
(``flume.prob_lo``), with slopes :math:`s_x` and :math:`s_y`
(``flume.bed.slope``):

.. math::

   z_b(x, y) = z_0 + s_x\,(x - x_0) + s_y\,(y - y_0) ,

evaluated at the cell centres. The distance from the corner is computed from
the cell index, :math:`x - x_0 = (i + 1/2)\,\Delta x`, so it does not lose
digits at large coordinates.

The initial water surface is level, at :math:`\eta_0`
(``flume.init.surface_elevation``), and the water moves with one velocity
:math:`(u_0, v_0)` (``flume.init.velocity``):

.. math::

   h = \max(\eta_0 - z_b,\; 0), \qquad hu = h\,u_0, \qquad hv = h\,v_0 .

A cell whose bed is at or above :math:`\eta_0` starts dry.

Time step
---------

A gravity wave moves at :math:`c = \sqrt{g h}` relative to the water, where
:math:`g` is the acceleration of gravity (``flume.gravity``). With the
velocity :math:`u = hu / h`, :math:`v = hv / h`, the largest rate at which
waves cross a cell is

.. math::

   r = \max_{\text{cells with } h > 0}
       \left(\frac{|u| + c}{\Delta x} + \frac{|v| + c}{\Delta y}\right) ,

where :math:`\Delta x` and :math:`\Delta y` are the cell sizes. The maximum
is taken over all boxes and all MPI ranks. The time step is

.. math::

   \Delta t = \frac{C}{r} ,

with the CFL number :math:`C` (``flume.cfl``, :math:`0 < C \le 1`).

The rates of the two directions are added because the scheme will update a
cell with the fluxes of both directions in one step (unsplit). For a
first-order upwind scheme of that kind, :math:`C \le 1` is the linear
stability bound. A second-order scheme that keeps the depth positive needs
about half of that, which is the default, 0.5. The limit for the scheme that
is chosen is settled when it is added.

A step that reaches ``flume.stop_time`` is shortened to end on it. A step that
falls short of it by round-off only (4 machine epsilons of the stop time, in
double precision) is lengthened to end on it, so that no step of that size
follows. The time and ``flume.stop_time`` are held in double precision in
every build: a single-precision sum of steps stops growing after about
:math:`2^{24}` steps.

A cell with :math:`h = 0` has no wave and does not enter the maximum; the
velocity is not computed there. If no cell holds water, :math:`r = 0` and no
step follows from the formula. The run then takes one step to
``flume.stop_time``; with ``flume.max_step`` alone it stops with a message
that asks for ``flume.stop_time``.

The run also stops with a message if the step is zero, not finite, or too
small to change the time.
