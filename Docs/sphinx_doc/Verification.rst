Verification
============

Each case below is a CTest regression test that compares a run with a closed
form. The numbers are what the runs give in double precision; the tests bound
them slightly above. :doc:`TestCatalog` lists every test.

Lake at rest
------------

A level surface over a rough bed (a plane moved by up to 0.5 m in every
cell), between walls, for 50 steps (``Tests/LakeAtRest``). Depth and
discharges stay at their initial values within 16 machine epsilons of the
surface elevation, at both orders, with dry islands, with a dry depth of
0.05 m, with outflow segments, with Manning friction, and over raster
buildings. Without the bed-slope term of the hydrostatic reconstruction the
discharges reach 1 m\ :sup:`2`/s.

Dam break
---------

Still water 2 m deep on one side of a dam, on a flat bed, 400 cells of
0.25 m, compared at 5 s with the exact solutions of Delestre et al. (2013,
sections 4.1.1 and 4.1.2) (``Tests/DamBreak``).

**Wet bed** (1 m on the low side). The solution is a rarefaction, a constant
state :math:`(h_m, u_m)` and a shock of speed :math:`S`:

.. math::

   u_m = 2\,(c_1 - c_m) = (h_m - h_0)\sqrt{\frac{g\,(h_m + h_0)}{2\,h_m h_0}} ,
   \qquad S = \frac{h_m u_m}{h_m - h_0} ,

with :math:`h_1`, :math:`h_0` the deep and shallow depths. Here
:math:`h_m` = 1.4538 m, :math:`u_m` = 1.3058 m/s, :math:`S` = 4.1831 m/s.

**Dry bed.** Ritter's solution; the front has run 44.3 m at 5 s.

.. list-table:: Dam break: L1 errors
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Case
     - Depth, order 1
     - Discharge, order 1
     - Depth, order 2
     - Discharge, order 2
   * - Wet bed
     - 0.80 %
     - 1.59 %
     - 0.24 %
     - 0.50 %
   * - Dry bed
     - 0.46 %
     - 0.39 %
     - 0.13 %
     - 0.12 %
   * - Dry bed, outflow side, 9 s
     - 0.45 %
     - 0.30 %
     - 0.10 %
     - 0.07 %
   * - Wet bed, outflow side, 10 s
     - 0.67 %
     - 1.21 %
     - 0.20 %
     - 0.32 %

On the wet bed the constant state is within 0.19 % (order 1) and 0.03 %
(order 2), and the shock within 0.2 and 0.1 cell of its place. Mass is
conserved to 16 machine epsilons. On the dry bed no depth is negative and no
water lies beyond the exact front; the thin front of the first-order run
lags the exact one, the depth passing 1 mm about 7 m behind it.

The first step of the dry-bed run has a closed form, which pins the dry-bed
wave speeds: with ``flume.cfl`` = 0.5 the first dry cell ends 1/3 m deep and
the last wet cell 5/3 m deep, within 64 machine epsilons.

With the outflow side the dry-bed front leaves supercritical and Ritter's
solution still holds. On the wet bed the flow behind the shock is
subcritical, and the wave the side sends back (:doc:`theory/Boundaries`)
leaves the constant state 0.17 % low in depth and 0.28 % high in velocity at
first order, 0.68 % and 1.97 % at second order.

Thacker's solutions in a bowl
-----------------------------

Water in a paraboloid bowl (:math:`h_0` = 0.1 m, :math:`a` = 1 m) with a
shoreline moving over the dry bed, compared with the solutions of Delestre
et al. (2013, section 4.2.2) (``Tests/Bowl``). With :math:`X`, :math:`Y`
measured from the centre and :math:`z = -h_0 (1 - r^2 / a^2)`:

**Planar** (a plane surface that turns about the axis; :math:`\eta` = 0.5 m,
period 4.4857 s):

.. math::

   h = \max\left(\frac{\eta h_0}{a^2}\left(2 X \cos\omega t
       + 2 Y \sin\omega t - \eta\right) - z,\; 0\right) , \qquad
   (u, v) = \eta\,\omega\,(-\sin\omega t,\; \cos\omega t) , \qquad
   \omega = \frac{\sqrt{2 g h_0}}{a} .

**Radial** (a paraboloid surface that rises and falls; :math:`r_0` = 0.8 m,
period 2.2429 s), with :math:`A = (a^2 - r_0^2) / (a^2 + r_0^2)` and
:math:`D = 1 - A \cos\omega t`:

.. math::

   h = \max\left(h_0\left[\frac{\sqrt{1 - A^2}}{D} - 1
       - \frac{r^2}{a^2}\left(\frac{1 - A^2}{D^2} - 1\right)\right] - z,\; 0\right) ,
   \qquad
   (u, v) = \frac{\omega A \sin\omega t}{2 D}\,(X, Y) , \qquad
   \omega = \frac{\sqrt{8 g h_0}}{a} .

.. list-table:: Thacker's solutions: L1 errors over the wet cells, relative
                to the mean initial depth (and :math:`\sqrt{g h_0}`)
   :header-rows: 1
   :widths: 30 14 14 14 14 14

   * - Run
     - Cells
     - Depth, order 1
     - Discharge, order 1
     - Depth, order 2
     - Discharge, order 2
   * - Planar, half a period
     - 100 x 100
     - 13.5 %
     - 19.7 %
     - 1.71 %
     - 2.12 %
   * - Planar, half a period
     - 200 x 200
     - 7.3 %
     - 10.8 %
     - 0.61 %
     - 0.71 %
   * - Planar, three periods
     - 100 x 100
     - 56.5 %
     - 76.8 %
     - 3.67 %
     - 5.56 %
   * - Radial, half a period
     - 100 x 100
     - 4.2 %
     - 2.8 %
     - 1.02 %
     - 0.71 %
   * - Radial, half a period
     - 200 x 200
     - 2.2 %
     - 1.8 %
     - 0.29 %
     - 0.20 %
   * - Radial, three periods
     - 100 x 100
     - 20.7 %
     - 4.5 %
     - 3.80 %
     - 1.54 %

The first-order scheme damps the motion; its error halves with the cell
size. At second order the error falls by 2.8 (planar) and 3.5 (radial) when
the cell size is halved, rates of 1.5 and 1.8. In every run no depth is
negative and mass is conserved to 64 machine epsilons.

Simple wave: rate of convergence
--------------------------------

One sine wave of 0.1 m on water 1 m deep over a flat bed, 10 m long and
periodic, in water moving at (0.5, -0.3) m/s, compared at 1 s with the exact
solution (``Tests/Wave``). The depth at :math:`x` at time :math:`t` is the
initial depth at :math:`x_0`, where

.. math::

   x_0 + \bigl(u_0 + 3 c(x_0) - 2 c_H\bigr)\, t = x , \qquad
   c_H = \sqrt{g H} .

The solution is smooth until 3.2 s.

.. list-table:: Simple wave: mean error, relative to the amplitude
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Run
     - Cells along
     - Steps
     - Depth
     - Discharge
   * - Order 1, x
     - 100
     - 154
     - 0.0353
     - 0.0408
   * - Order 2, x
     - 100
     - 154
     - 0.00350
     - 0.00413
   * - Order 2, x
     - 200
     - 307
     - 0.00094
     - 0.00111
   * - Order 2, x
     - 400
     - 614
     - 0.00025
     - 0.00030
   * - Order 2, y
     - 100
     - 142
     - 0.00271
     - 0.00252
   * - Order 2, y
     - 200
     - 283
     - 0.00073
     - 0.00068
   * - Order 2, x, cells 250 m wide
     - 100
     - 82
     - 0.00364
     - 0.00428

The second-order error falls by 3.7 each time the cell size is halved, a
rate of 1.9; the minmod limiter keeps it below 2. The first-order run is ten
times further from the solution on the same cells.

Uniform flow with friction
--------------------------

Water 1 m deep starts at rest in an endless channel that falls by 1 mm per
metre, with :math:`n` = 0.03 s/m\ :sup:`1/3` (``Tests/Friction``). The depth
stays 1 m and the discharge follows

.. math::

   q(t) = q_n \tanh\left(\frac{t}{T_r}\right) , \qquad
   q_n = \frac{h^{5/3} S^{1/2}}{n} , \qquad T_r = \frac{q_n}{g h S} ,

with :math:`q_n` = 1.05409 m\ :sup:`2`/s and :math:`T_r` = 107.451 s.

.. list-table:: Friction channel: largest relative error of the discharge
   :header-rows: 1
   :widths: 60 20 20

   * - Run
     - Order 1
     - Order 2
   * - Transient, 100 s
     - 5.64e-4
     - 1.66e-4
   * - Transient, 100 s, half the step
     -
     - 8.31e-5
   * - Steady, 1500 s
     - 2.50e-4
     - 1.2e-13
   * - Thin film: 1 mm on a slope of 0.1, :math:`n` = 0.05, 10 s
     -
     - 1.1e-9
   * - Manning map, one coefficient per row (0.02 to 0.06)
     - 5.63e-4
     - 1.91e-4

The steady error at first order is the bias of the first-order bed-slope
term, :math:`1 - \sqrt{1 - S \Delta x / (2 h)}` = 2.50031e-4
(:ref:`sec:Friction`). The transient error halves with the step: friction is
first order in time. In the thin film the friction time is 0.0645 s and the
step about 2 s; an explicit friction term would reverse the discharge.

Rain on a plane: kinematic wave
-------------------------------

Rain of :math:`i` = 50 mm/h on a dry, impermeable plane 100 m long, falling
at :math:`S_0` = 0.01 to an outflow side, with :math:`n` = 0.01
s/m\ :sup:`1/3` (``Tests/Rain/inputs_plane``). With
:math:`\alpha = S_0^{1/2} / n` the kinematic-wave solution is

.. math::

   h(x, t) = \begin{cases} i\, t, & t < t_x \\
             (i\, x / \alpha)^{3/5}, & t \ge t_x \end{cases}
   \qquad
   q = \alpha\, h^{5/3} , \qquad
   t_x = \left(\frac{x}{\alpha\, i^{2/3}}\right)^{3/5} .

The plane reaches equilibrium at 349 s, with 4.8 mm at the outlet and a
Froude number of 1.3.

.. list-table:: Kinematic wave: largest error to 900 s
   :header-rows: 1
   :widths: 40 30 30

   * - Order, cell, CFL
     - Depth, L1
     - Outlet discharge
   * - 1, 0.2 m, 0.5
     - 10.2 %
     - 11.1 %
   * - **1, 0.1 m, 0.5 (test)**
     - **5.5 %**
     - **5.7 %**
   * - 2, 0.2 m, 0.5
     - 2.8 %
     - 3.6 %
   * - 2, 0.1 m, 0.5
     - 1.2 %
     - 2.6 %
   * - **2, 0.1 m, 0.125 (test)**
     - **0.22 %**
     - **1.3 %**

At order 1 the error falls with the cell size and is largest at the upstream
end, where the film is thinner than the bed step of a cell
(:ref:`sec:HydrostaticReconstruction`). At order 2 and CFL 0.5 the error is
mostly that of friction acting once per step, which halves with the step. At
order 1 the outlet discharge at equilibrium equals :math:`i L` to round-off.

A coarser set-up (1 m cells, :math:`n` = 0.03) is outside what the scheme
resolves: the bed falls by more than the film depth per cell (depth error
27 % at equilibrium) and the outflow is subcritical.

Rain volume
-----------

Rain on a closed, flat, dry domain (``Tests/Rain/inputs_flat``): every cell
holds the rain fallen by the plotfile's time, equal to the last bit in
double precision, for a constant rate and for step and linear hyetographs.
On 1 m of still water, rain of 0.36 mm/h (1.6e-8 m per step) is kept in
single precision by the compensated sum. In a closed basin with moving
water the volume grows by the rain fallen (8.9e-13 of an added 6.4e-3).

Rain that varies in space (``Tests/RainField/inputs_flat``, 24 by 16 cells
of 2 m in 6 boxes), against a reference computed from the files:

.. list-table::
   :header-rows: 1
   :widths: 64 36

   * - Case
     - Largest difference
   * - Every cell, one step, grids of 3.7 by 2.9 m, linear in time
     - 9.9e-16 of the largest depth
   * - Every cell, one step, grid of 1.3 m
     - 4.3e-16
   * - Every cell, 1 m grid on solver edges, NODATA ring outside
     - 17 eps
   * - Water to 1200 s, three grids, step and linear
     - 7.1e-16 of the water (3 eps)
   * - Moving cell at rest, every cell
     - 5.8e-16
   * - Moving cell, 29 m in one step (20 sub-intervals)
     - 1.1e-13 (1.3e-2 with one)
   * - Moving cell leaving the domain, water to 1200 s
     - 1.4e-16

.. _sec:VerifyInfiltration:

Infiltration
------------

Closed, flat domains where every cell is the same, against the closed forms
of :ref:`sec:Infiltration` (``Tests/Infiltration``). Soil parameters from
Rossman and Huber (2016, Table 4-7).

.. list-table::
   :header-rows: 1
   :widths: 64 36

   * - Case
     - Largest difference
   * - Green-Ampt, 0.5 m ponded on silt loam, 6 h, 4506 steps
     - 2.9e-11 m (a time shift of 3.7e-6 s)
   * - Green-Ampt, 36 mm/h on dry sandy loam, 2 h
     - 7.6e-17 m
   * - Curve number 75, 36 mm/h, 2 h (72 mm, loss 50.30 mm)
     - 1.4e-17 m
   * - Curve number 65.31 with :math:`\lambda = 0.05`
     - 1.4e-17 m
   * - Curve-number map with flow and buildings, per cell
     - 3.5e-18 m

The ponded case differs only by a shift in time: the first step takes the
ponded water as rain and saturates the surface after
:math:`K_s M \Delta t^2/h^2 \approx 7.5 \times 10^{-6}` s, where exact
ponding takes half that, so :math:`F` lags by
:math:`K_s M \Delta t^2/(2h^2) \approx 3.7 \times 10^{-6}` s; from there
the integrated equation is exact. In every case with a closed form, the
depth plus the infiltrated depth of each cell is the water it got to 4
machine epsilons, also over the 279217 steps of a curve-number run on 2 cm
cells. On a sloped, walled
domain with a hyetograph, a conductivity map and buildings (second order, 3
ranks), the water on the surface plus the infiltrated water equals the rain
to :math:`5.7 \times 10^{-14}` of 28.8 m\ :sup:`3` per m\ :sup:`2` of cell.

.. _sec:VerifyDrains:

Storm-drain inlets
------------------

The inlets of :ref:`sec:Drains` (``Tests/Drains``) are checked against the
exact solution of :math:`dh/dt = -\min(a h^{3/2}, b h^{1/2})`, which the
check evaluates through its inverse, the time to reach a depth. With level,
still water over an inlet in each of 4 cells of 100 m\ :sup:`2`
(:math:`P = 2` m, :math:`A = 0.5` m\ :sup:`2`), every cell holds that
solution at the time of each plotfile:

.. list-table::
   :header-rows: 1
   :widths: 64 36

   * - Case
     - Largest difference
   * - From 0.8 m: orifice to :math:`d^* = 0.4469` m at 30.4 s, then weir, to
       8.34 mm at 600 s
     - 0.8 epsilons of 0.8 m
   * - From 0.3 m, below :math:`d^*`: weir only, to 7.20 mm at 600 s
     - 1.0 epsilon of 0.3 m
   * - :math:`C_w = 1.5`, :math:`C_o = 0.6` (:math:`d^* = 0.4429` m)
     - 3.4 epsilons of 0.8 m
   * - Drained volume of the first case, 316.6641978 m\ :sup:`3`
     - 1e-15 relative

One inlet (:math:`P = 0.5` m, :math:`A = 0.05` m\ :sup:`2`) in a closed
basin of 256 m\ :sup:`2`, still water 0.3 m deep, drains 0.979 to 0.987
times the volume of the closed form for one cell of the basin's area over
20 min (0.976 to 0.984 at second order). The lumped form is an upper bound
while the inlet cell is the lowest point of the water, since the rate grows
with the depth; the check takes 0.95 as the lower bound, from the run. The run is the
same to the last bit on 3 ranks in 16 boxes, and after a restart from a
checkpoint on 2 ranks. A 1 mm film over an inlet of :math:`P = 100` m on a
1 m cell, which an explicit step at the weir rate would let take 13 mm,
keeps every depth at or above 0; the depth it drains in the first step is
that of the exact solution to 0.2 epsilons of 1 mm. In every plotfile of these runs the depth
plus the drained depth summed over the cells is that of the start to 8
epsilons; in the basin, where the water moves and every step rounds each
wet depth, to :math:`8 + \sqrt{n}` epsilons after :math:`n` steps (the
runs drift by 1 epsilon in double precision and 17 in single precision
over 6095 steps).

The approach-flow law is checked against the exact solution of
:math:`dh/dt = -c h^{p}`, through its inverse
:math:`t = (h_0^{q} - h^{q}) / (q c)`. Sixteen periodic cells of 25
m\ :sup:`2`, each with an inlet of 0.4 m\ :sup:`2` from a list file, hold
water 0.5 m deep that moves at (1.2, 0.5) m/s. The water stays level and
keeps its velocity, so every cell holds that solution at the time of each
plotfile, and its discharge is its depth times the velocity:

.. list-table::
   :header-rows: 1
   :widths: 64 36

   * - Case
     - Largest difference
   * - :math:`k = 0.302`, :math:`\alpha = -0.816`: to 69.43 mm at 60 s
     - 1.3 epsilons of 0.5 m
   * - :math:`k = 0.226`, :math:`\alpha = -0.5`: to 171.3 mm at 60 s
     - 2.9 epsilons of 0.5 m
   * - :math:`\alpha = 0`: to 419.5 mm at 60 s
     - 0.6 epsilons of 0.5 m
   * - To 120 s in 145 steps: every cell empty after :math:`t_e = 87.05` s,
       none before the closed form reaches the dry depth, no depth below 0
     - 1.6 epsilons of 0.5 m
   * - The same in 6492 steps (a CFL number of 0.011)
     - 11.6 epsilons of 0.5 m
   * - Still water, :math:`\alpha = -1`: to 64.05 mm at 60 s
     - 0.7 epsilons of 0.5 m
   * - Two lists with other :math:`k` and :math:`\alpha` in one run, first
       step: each cell the depth of its own inlet
     - 0.2 epsilons of 0.5 m
   * - Drained volume of the first case, 172.2272660 m\ :sup:`3`
     - 1e-12 relative (the bound of the check)

The checks take 64 epsilons: the powers of this closed form round more than
the square roots of the weir law, and the error grows with the number of
steps. In the run of 6492 steps a step lands in the last 0.04 s before
:math:`t_e`, where the closed form is at or below the dry depth; the step
takes that film and the cells end empty. With still water and the default exponent
no inlet takes water, and every depth stays that of the start to the last
bit. Two inlets of 0.2 m\ :sup:`2` in every cell give the plotfiles of the
first case to the last bit. One inlet of 0.2 m\ :sup:`2` in a closed basin
of 256 m\ :sup:`2`, in water 0.3 m deep that starts at (0.5, 0.2) m/s and
sloshes between the walls, drains 3.950 m\ :sup:`3` of 76.8 in 60 s (4.310
at second order), only from its cell; the run is the same to the last bit
on 3 ranks in 16 boxes with the inlet from a list file, and after a restart
on 2 ranks (the restart reads the list from a file with a byte-order mark
and line ends of CR and LF). This case has no closed form: 3.950 and 4.310
are the values of this grid at a CFL number of 0.5, and the first is kept
as a regression value. The weir and orifice law is as it
was: the basin above with ``type = weir_orifice`` given, with its inlet
from a list file, and with the inlet as two halves in one cell of a list
file, gives the same plotfiles to the last bit.

The approach-flow tests check the step against its own differential
equation. One published pair of numbers is checked: for their grate of
0.3375 m\ :sup:`2` in water 0.4 m deep, Xia et al. (2022, section 3.4)
give 0.1166 m\ :sup:`3`/s at 0.1 m/s and 0.1570 at 0.5 m/s; the rate of
the check script gives 0.1166 and 0.1567 (:math:`g = 9.81` m/s\ :sup:`2`),
and one step of FLUME with that grate at 0.5 m/s follows the closed form of
that rate to 64 epsilons. Otherwise the tests do not compare the law with
measurements: Annis et al.
(2026) calibrated it against street sensors, and Xia et al. (2022) against
laboratory runs.

Water balance and flood maps
----------------------------

Every regression run that leaves a water-balance log has its residual
checked (:doc:`theory/WaterBalance`), and its residual is recomputed from
the columns. The largest change of the residual between two lines is 14.5
times :math:`\varepsilon n (V + M) + \varepsilon_d C` in double precision
(``rain_dry_start_discharge``) and 15.3 times in single precision
(``infiltration_green_ampt_rain``), against the bound of 64. The columns are pinned against closed forms
(``Tests/WaterBalance``):

.. list-table::
   :header-rows: 1
   :widths: 46 30 24

   * - Case
     - Closed form
     - Difference
   * - Rain on a closed basin, 36 mm/h for 600 s on 1536 m\ :sup:`2`
     - 9.216 m\ :sup:`3` rained and stored
     - 7.1e-15 m\ :sup:`3`
   * - Ponded Green-Ampt infiltration, 6 h on 102400 m\ :sup:`2`
     - 9395.3755 m\ :sup:`3` infiltrated
     - 1.1e-6 m\ :sup:`3`
   * - Uniform flow at 10 m/s, 1 m deep, out of a 20 m side for 2 s
     - 400 m\ :sup:`3` out
     - 4.5e-13 m\ :sup:`3`
   * - The same out of segments of 5 m and 15 m (second order)
     - 100 and 300 m\ :sup:`3`
     - 4.5e-13 m\ :sup:`3`
   * - Uniform flow in a periodic domain
     - 12288 m\ :sup:`3` stored
     - 0

The infiltration difference is the 3.7e-6 s lag of the first step
(Infiltration, above). With rain, infiltration and outflow
together on the kinematic-wave plane (second order, 5315 steps), the
residual ends at -5.6e-16 m\ :sup:`3`. The log without any one of its three
terms fails the check on at least 90 % of the pairs of lines between which
that term changes (in fact on all of them in double precision, and on all
but 4 of 5268 for the outflow in single precision), also with a line every
25 steps.

On the wet-bed dam break, with an arrival depth of 1.25 m between the two
depths, the arrival time of each cell behind the shock is within 0.024 s
(order 1) and 0.014 s (order 2) of :math:`(x - 50\,\mathrm{m})/S`, against
a cell crossing time of 0.060 s plus one step. The maximum depth is within 0.40 % and
0.027 % of :math:`h_m`, and the maximum speed and depth times speed within
1.5 % and 0.11 % of :math:`u_m` and :math:`h_m u_m`. On the rain basin every
cell's maximum depth is the rain fallen, its speeds are 0, and water arrives
one step after the time the closed form gives. On a sloped bed with water in
motion, three gauges in three boxes on 2 ranks hold exactly the values of
their cells in the plotfile of the same step.

Partial outflow
---------------

The dry-bed dam break on 200 by 40 cells, with the side ahead an outflow
over half its length and a wall over the rest. In the open rows the L1 error
against Ritter's solution is 0.75 % (order 1) and 0.20 % (order 2) of the
depth. In the wall rows the depth at the wall is within 4.3 % of the bore
depth from the jump conditions.

Inflow hydrograph
-----------------

A channel 200 m long and 4 m wide, two cells across and 0.5 m cells along
it (1 m in the uniform and normal-depth runs), flat and without friction
unless stated, fed through a discharge segment on one side
(:ref:`sec:Inflow`, ``Tests/Inflow``), with :math:`g` = 9.81 m/s\ :sup:`2`.

* **Uniform.** Water 1 m deep at 1 m/s, fed with 1 m\ :sup:`2`/s, is a
  steady state. After 300 s (3419 steps) on 4 ranks every cell holds it
  within 250 machine epsilons in the depth and 1030 in the discharge at
  both orders.
* **Rarefaction.** The same flow when the discharge drops to
  0.5 m\ :sup:`2`/s. The invariant :math:`u - 2\sqrt{g h}` is kept across
  the rarefaction, so next to the side
  :math:`0.5 / h_b - 2\sqrt{g h_b} = 1 - 2\sqrt{g}` and
  :math:`h_b` = 0.869019 m. The same on each of the four sides.
* **Dry fan.** 1 m\ :sup:`2`/s onto a dry bed. The face is at the critical
  depth :math:`h_c` = 0.467136 m and a rarefaction runs from it to the dry
  front, :math:`c = c_c - x / (3t)`, :math:`c_c = \sqrt{g h_c}`, compared
  at 30 s over the first 60 m.
* **Normal depth.** 1 m\ :sup:`2`/s down a slope of 0.02 with
  :math:`n` = 0.03 s/m\ :sup:`1/3`, 1 m cells, from dry to 600 s:
  :math:`h_n = (q n / \sqrt{S})^{3/5}` = 0.394424 m (Froude number 1.29),
  compared from 40 m to 190 m.

.. list-table:: Inflow channel: largest relative error
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Run
     - Depth, order 1
     - Discharge, order 1
     - Depth, order 2
     - Discharge, order 2
   * - Rarefaction, 20 s, first 50 m
     - 2.2e-4
     - 6.8e-4
     - 5.6e-5
     - 2.8e-4
   * - Dry fan, 30 s, first 60 m
     - 1.07e-2
     - 2.5e-3
     - 7.5e-3
     - 1.1e-3
   * - Normal depth, 600 s
     - 7.7e-3
     - 4.6e-13
     - 2.1e-3
     - 3.5e-3

The fan error is largest next to the side, where the fan starts at the
critical speed. The first-order normal depth is raised by the bias of the
first-order bed-slope term, :math:`g S (h - S \Delta x / 2)`: with 1 m
cells it gives :math:`h` = 0.397451 m, 0.767 % above :math:`h_n`; at second
order the depth error falls to 1.1e-3 and 6e-4 with cells of 0.5 m and
0.25 m.

On a closed basin of 40 by 40 m with 0.5 m of water and two segments, the
water that has entered by 300 s equals the integral of each hydrograph:

.. list-table:: Inflow basin: water entered and stored at 300 s, m\ :sup:`3`
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Run
     - Segment a
     - Segment b
     - Stored
     - Difference
   * - Triangle linear, steps held, 4 ranks
     - 200
     - 260
     - 1260
     - 2.3e-13
   * - Triangle held, steps linear
     - 280
     - 205
     - 1285
     - 0
   * - First entry at 50 s
     - 600
     - 0
     - 1400
     - 0

The difference is that of the stored water from 800 m\ :sup:`3` plus the
two integrals. The steps do not land on the times of the files.

Rasters
-------

A bed raster holding :math:`z = 2 + 0.03 x - 0.02 y + 0.001 x y` at the
centres of 5 m cells is reproduced on 1 m cells to 8.9e-16 m
(``Tests/Raster``). A dam break through streets between raster buildings 3
to 8 m high keeps every building cell dry and conserves mass to 64 machine
epsilons; with the buildings 0.3 m high the water goes over them and mass is
conserved.

Decomposition
-------------

Runs on 1, 2 and 4 ranks, in boxes down to one cell wide, give the same
cells. An outflow side given as one segment or as two that meet inside a box
gives plotfiles identical to the last bit.
