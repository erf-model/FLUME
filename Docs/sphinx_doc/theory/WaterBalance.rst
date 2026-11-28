Water balance
=============

The water on the surface changes only through the sources (rain in,
infiltration and drain inlets out) and through the faces of the outflow,
zero-gradient and discharge segments. The
balance log (:ref:`sec:WaterBalanceLog`) counts each of these where it
happens and checks that they add up.

Residual
--------

With :math:`V(t)` the volume on the surface at time :math:`t`,
:math:`V_0 = V(0)`, :math:`S_n(t)` the volume source :math:`n` has added
since the start (below 0 for a sink) and :math:`Q_m(t)` the volume that has
left through outflow or zero-gradient segment :math:`m` (below 0 where
water has come in through a zero-gradient segment)
and :math:`P_l(t)` the volume that has entered through discharge segment
:math:`l`, the residual is

.. math::

   R(t) = V(t) - V_0 - \sum_n S_n(t) + \sum_m Q_m(t) - \sum_l P_l(t) ,

which is zero in exact arithmetic.

* :math:`V = A \sum_{i,j} h_{ij}`, with :math:`A = \Delta x \, \Delta y` the
  cell area and the sum over the valid cells.
* :math:`S_n` is counted by the source itself, not from :math:`V`: the rain
  from the depth the provider gives each cell, less the part a compensated
  sum still holds back from the state (at most the rounding of a depth); the
  infiltration from the depth each cell has lost to the soil, with the part
  its compensated sum carries (:ref:`sec:Infiltration`); the drain inlets
  likewise (:ref:`sec:Drains`).
* :math:`Q_m = \sum_k \sum_{f \in m} \sigma_f \, F_f \, \ell_f \, w_k
  \Delta t_k`, with :math:`F_f` the mass flux (m\ :sup:`2`/s) the update
  uses on boundary face :math:`f` of the segment, after the outflow limit,
  :math:`\ell_f` the length of the face, :math:`\sigma_f = +1` on a high
  side and :math:`-1` on a low side, :math:`\Delta t_k` the step and
  :math:`w_k` the weight of each stage: 1 at first order and
  :math:`\tfrac{1}{2}` for each of the two stages at second order, whose
  mean is the step.
* :math:`P_l` is formed as :math:`Q_m` on the faces of the discharge
  segment, with the opposite sign: :math:`\sigma_f = -1` on a high side
  and :math:`+1` on a low side.

**Why count each term where it happens.** :math:`R` then measures how well
the update, the outflow limit and the sources keep the water they account
for. A term taken from the change of :math:`V` would make :math:`R` zero by
construction. Walls are not counted: the mass flux through a wall face is
zero up to the rounding of one product, and what it carries shows in
:math:`R`.

All sums are formed in double precision, also in a single-precision build.
The rain, outflow and inflow totals add a volume each step to a compensated (Kahan)
sum, so their rounding does not grow with the steps; the infiltration total
is summed again from the cells at each line.

Bound
-----

Each step rounds the depth of every cell it changes. Between two lines
:math:`a` and :math:`b` of the log, :math:`n` steps apart, the residual may
change by at most

.. math::

   |R_b - R_a| \le F \Bigl( \varepsilon \, n \bigl( \max(|V_a|, |V_b|) + M \bigr)
   + \varepsilon_d \, C \Bigr) ,

.. math::

   M = \sum_{c} |c_b - c_a| , \qquad
   C = |V_0| + \sum_{c} \max(|c_a|, |c_b|) ,

with :math:`c` running over the source and segment columns,
:math:`\varepsilon` the machine epsilon of the state, :math:`\varepsilon_d`
that of double precision and :math:`F = 64`.

**What the terms are.** :math:`n (V + M)`: each step may leave a few
roundings of the water it holds and moves. :math:`C`: :math:`R` is formed
from cumulative totals held in double precision, and each is rounded to a
relative :math:`\varepsilon_d`; this term does not grow with :math:`n`.

**Why per pair of lines.** A bound on the change between lines, not on
:math:`R` itself, keeps the test as sharp late in a run as early: a term left
out of the balance fails as soon as it moves more water between two lines
than the bound. The largest ratio of a change to
:math:`\varepsilon n (V + M) + \varepsilon_d C` over the regression tests
is 14.5 in double and 15.3 in single precision (:doc:`../Verification`),
against :math:`F = 64`.

**What it does not detect.** A term that moves less water between two
lines than the bound: per step, less than about
:math:`F \varepsilon (V + m)`, with :math:`m` the water moved in a step.
Such a term passes every pair of lines, and over :math:`T` steps it can add
up to :math:`T F \varepsilon (V + m)`. On the every-term test (5315 steps,
:math:`V = 0.49` m\ :sup:`3`) that is :math:`4 \times 10^{-11}`
m\ :sup:`3` in double precision and 0.02 m\ :sup:`3` in single
precision, against final residuals of :math:`6 \times 10^{-16}` and
:math:`10^{-7}` m\ :sup:`3`. The bound grows with the steps between lines,
so a larger ``flume.log_int`` detects less; the every-term test is also run
with a line every 25 steps.

Open sides that draw water in
-----------------------------

A ``zero_gradient`` side copies the boundary cell into its ghost cells
whatever the direction of the flow (:ref:`sec:ZeroGradient`). Where the
discharge of a boundary cell points into the domain, water enters; the log
then shows a falling :math:`Q_m`. On the sloped, walled domain of
``Tests/Infiltration/inputs_balance`` with a ``zero_gradient`` segment on
the low-x side, 10.2 m\ :sup:`3` come in through it over the run. The
residual stays at round-off: the water is counted. With an ``outflow``
segment there, 2.75 m\ :sup:`3` leave and none enters: the column of an
``outflow`` segment never falls (:ref:`sec:Outflow`).
