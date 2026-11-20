Model
=====

This page describes what the code holds and computes today: the
shallow-water equations over a fixed bed, advanced by a first-order
finite-volume scheme between walls or periodic sides. Wetting and drying
fronts, friction, rainfall and open boundaries are not written yet.

Equations
---------

FLUME solves the shallow-water equations for the depth :math:`h` and the
discharges :math:`hu` and :math:`hv` over a bed of elevation :math:`z_b`:

.. math::

   \frac{\partial h}{\partial t} + \frac{\partial (hu)}{\partial x}
      + \frac{\partial (hv)}{\partial y} &= 0 , \\
   \frac{\partial (hu)}{\partial t}
      + \frac{\partial}{\partial x}\left(hu^2 + \tfrac{1}{2} g h^2\right)
      + \frac{\partial (huv)}{\partial y} &= -g h \frac{\partial z_b}{\partial x} , \\
   \frac{\partial (hv)}{\partial t} + \frac{\partial (huv)}{\partial x}
      + \frac{\partial}{\partial y}\left(hv^2 + \tfrac{1}{2} g h^2\right)
      &= -g h \frac{\partial z_b}{\partial y} ,

where :math:`u` and :math:`v` are the depth-averaged velocities and :math:`g`
is the acceleration of gravity (``flume.gravity``). The only source term so
far is the bed slope.

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

Two more features of the initial state build test cases.

**Rough bed.** ``flume.bed.noise`` gives an amplitude :math:`A`, and the bed
of cell :math:`(i, j)` moves by :math:`A\,r_{ij}`, where :math:`r_{ij}` is a
fixed pseudo-random number in :math:`[-1, 1)` computed from the two cell
indices. The bed is the same for any boxes and any number of ranks.

**Dam.** ``flume.init.dam_position`` gives a coordinate :math:`x_d` along x
or y (``flume.init.dam_direction``). Cells whose centre is below :math:`x_d`
start under the surface :math:`\eta_0`; the others start under
``flume.init.dam_surface_elevation``. The velocity is the same on both sides
unless ``flume.init.dam_velocity`` gives another one for the second side.

Finite-volume step
------------------

A step updates every cell :math:`(i, j)` from the fluxes on its four faces
and the bed-slope term, in one stage (forward Euler) and without splitting
the directions. With :math:`U = (h, hu, hv)`:

.. math::

   U_{ij}^{n+1} = U_{ij}^{n}
     - \frac{\Delta t}{\Delta x}\left(F_{i+1/2,j} - F_{i-1/2,j}\right)
     - \frac{\Delta t}{\Delta y}\left(G_{i,j+1/2} - G_{i,j-1/2}\right)
     + \Delta t\, S_{ij} .

:math:`F` and :math:`G` are the fluxes across the faces in x and y, and
:math:`S` is the bed-slope term, both defined below. The scheme is first order
in space and time. The fluxes are computed once per face and stored, so the
water that leaves a cell through a face is the water that enters its
neighbour.

Hydrostatic reconstruction
--------------------------

The bed is a step at each face. The face between a cell with depth
:math:`h` and bed :math:`z_b` and a neighbour with bed :math:`z_b'` sees the
cell through the depth (Audusse et al., 2004)

.. math::

   h^{*} = \max\bigl(0,\; h - \max(0,\; z_b' - z_b)\bigr) .

This is the part of the cell's water that lies above the higher of the two
beds. It is never negative. Only the difference of the two bed elevations
enters, so large elevations cost no digits. The flux on the face is computed
from :math:`h^{*}` on each side and the velocities of the two cells.

The bed-slope term of a cell uses the same depths. With :math:`h^{*}_{E}`,
:math:`h^{*}_{W}`, :math:`h^{*}_{N}` and :math:`h^{*}_{S}` the depths the cell
presents to its high-x, low-x, high-y and low-y faces,

.. math::

   S^{hu}_{ij} = \frac{g}{2\,\Delta x}\left(h^{*2}_{E} - h^{*2}_{W}\right) ,
   \qquad
   S^{hv}_{ij} = \frac{g}{2\,\Delta y}\left(h^{*2}_{N} - h^{*2}_{S}\right) ,

and :math:`S^{h} = 0`.

**Why.** For water at rest under a level surface, the two cells at a face
present the same :math:`h^{*}`, the mass flux is zero and the momentum flux is
:math:`\tfrac{1}{2} g h^{*2}`. The term :math:`S` then cancels the difference
of the momentum fluxes of a cell, and the water stays at rest, to round-off,
whatever the bed (a well-balanced scheme). A cell next to a dry cell whose bed is above
the surface presents no water to that face, and nothing flows.

**Limit.** Where the bed rises at a face by at least the depth of the lower
cell, :math:`h^{*}` is zero and the face sees less of the bed step than there
is. On a steep bed under shallow moving water the scheme then underestimates
the slope: the velocity comes out too low and the depth too high (Delestre et
al., 2012). The bed step is the slope times the cell size, so the effect goes
away as the cells get smaller. Water at rest is not affected.

HLLC flux
---------

At a face, let L and R be the low and high sides, :math:`u_n` the velocity
normal to the face, :math:`u_t` the velocity along it, and
:math:`c = \sqrt{g h^{*}}`. The fastest waves to the two sides are estimated
as

.. math::

   s_L = \min(u_{n,L} - c_L,\; u^{\star} - c^{\star}) , \qquad
   s_R = \max(u_{n,R} + c_R,\; u^{\star} + c^{\star}) ,

.. math::

   u^{\star} = \tfrac{1}{2}(u_{n,L} + u_{n,R}) + c_L - c_R , \qquad
   c^{\star} = \max\bigl(0,\; \tfrac{1}{2}(c_L + c_R)
               + \tfrac{1}{4}(u_{n,L} - u_{n,R})\bigr) .

These follow from taking both waves as rarefactions. Across the wave on the
low side :math:`u_n + 2c` keeps its value, and across the wave on the high
side :math:`u_n - 2c` does. The state between the two waves then has
:math:`u^{\star} + 2 c^{\star} = u_{n,L} + 2 c_L` and
:math:`u^{\star} - 2 c^{\star} = u_{n,R} - 2 c_R`, which gives
:math:`u^{\star}` and :math:`c^{\star}` above. Each estimate is the faster
of the edge of the wave in the cell's own state and its edge in the middle
state.

When one side has no water the middle state is dry, :math:`c^{\star} = 0`,
and the speeds are those of a front on a dry bed:
:math:`s_L = u_{n,L} - c_L`, :math:`s_R = u_{n,L} + 2 c_L` with the high side
dry, and :math:`s_L = u_{n,R} - 2 c_R`, :math:`s_R = u_{n,R} + c_R` with the
low side dry. A face that is dry on both sides has no flux.

With :math:`f = (h^{*} u_n,\; h^{*} u_n^2 + \tfrac{1}{2} g h^{*2})` the mass
and normal-momentum fluxes of each side, the flux on the face is

.. math::

   f^{\mathrm{HLL}} =
   \begin{cases}
     f_L , & s_L \ge 0 , \\
     \dfrac{s_R f_L - s_L f_R + s_L s_R (q_R - q_L)}{s_R - s_L} ,
           & s_L < 0 < s_R , \\
     f_R , & s_R \le 0 ,
   \end{cases}

where :math:`q = (h^{*},\; h^{*} u_n)`. The momentum along the face is carried
by the mass flux from the side it comes from:

.. math::

   f_{t} = f^{\mathrm{HLL}}_{h}\; \times
   \begin{cases} u_{t,L} , & f^{\mathrm{HLL}}_{h} \ge 0 , \\
                  u_{t,R} , & f^{\mathrm{HLL}}_{h} < 0 . \end{cases}

This last choice is what the C of HLLC adds: it keeps a shear layer that lies
along a face as sharp as upwinding allows. For the shallow-water equations
the middle wave of the HLLC flux moves with the sign of the HLL mass flux, so
the flux above is the HLLC flux written without the speed of that wave.

**Limit.** The two speeds are exact when both waves are rarefactions. For a
shock that runs into much shallower water they are below the exact shock
speed: for still water next to a very thin layer the estimate is 1.5 times
the celerity of the deep side where the exact front moves at about 1.7
times. In a check over random pairs of wet states, 13 % had an estimated
speed inside the exact wave fan by more than 1 %. The tests of
:doc:`Testing` pass with these speeds. Speeds built from the depth ratio
across the shock remove the gap but need a depth below which a
cell counts as dry, which this version does not have.

The velocity of a cell is :math:`hu / h` where :math:`h > 0` and zero in a dry
cell. There is no depth below which a cell counts as dry yet, so a very thin
layer of water can carry a large velocity. This is treated with the wetting
and drying fronts.

Sides of the domain
-------------------

Each direction is periodic or closed by walls (``flume.periodic``). The state
and the bed carry one layer of ghost cells. Between boxes and across a
periodic side a ghost cell holds the value of the cell it covers. Outside a
wall it mirrors the first cell inside: same depth and bed, same discharge
along the wall, and the discharge normal to the wall with its sign changed.
The mass flux through a wall is then exactly zero. A sloping bed is taken as
it is at a periodic side: the bed has a step there.

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

The rates of the two directions are added because the scheme updates a cell
with the fluxes of both directions in one step (unsplit). :math:`C \le 1` is
the linear stability bound of the first-order scheme. The depth stays
non-negative when no cell loses more water in a step than it holds. This is
not guaranteed for every :math:`C \le 1`: the wave speeds :math:`s_L` and
:math:`s_R` at a face can be larger than the rate of either cell. A cell of
still water of depth :math:`h` with dry cells on its four sides loses
:math:`2 c h / 3` through each face and ends the step with
:math:`h\,(1 - 4C/3)` on square cells, which is negative for
:math:`C > 0.75`. For the default :math:`C = 0.5` a search over random
states found no negative depth; this is not a proof. The depth is not
clipped: if it becomes negative in any cell the run stops with a message
that asks for a lower ``flume.cfl``. The state is checked before every step and
before every plotfile.

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

References
----------

* E. Audusse, F. Bouchut, M.-O. Bristeau, R. Klein and B. Perthame (2004). A
  fast and stable well-balanced scheme with hydrostatic reconstruction for
  shallow water flows. *SIAM Journal on Scientific Computing* 25(6),
  2050-2065. https://doi.org/10.1137/S1064827503431090
* O. Delestre, S. Cordier, F. Darboux and F. James (2012). A limitation of
  the hydrostatic reconstruction technique for Shallow Water equations.
  *Comptes Rendus Mathematique* 350(13-14), 677-681.
  https://doi.org/10.1016/j.crma.2012.08.004
* A. Ritter (1892). Die Fortpflanzung der Wasserwellen. *Zeitschrift des
  Vereines Deutscher Ingenieure* 36(33), 947-954.
