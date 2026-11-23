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

Partial outflow
---------------

The dry-bed dam break on 200 by 40 cells, with the side ahead an outflow
over half its length and a wall over the rest. In the open rows the L1 error
against Ritter's solution is 0.75 % (order 1) and 0.20 % (order 2) of the
depth. In the wall rows the depth at the wall is within 4.3 % of the bore
depth from the jump conditions.

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
