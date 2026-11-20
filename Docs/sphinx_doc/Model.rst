Model
=====

This page describes what the code holds and computes today: the
shallow-water equations over a fixed bed, advanced by a first-order
finite-volume scheme with wetting and drying, between walls or periodic
sides. Friction, rainfall and open boundaries are not written yet.

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

Four more features of the bed and the initial state build test cases.

**Rough bed.** ``flume.bed.noise`` gives an amplitude :math:`A`, and the bed
of cell :math:`(i, j)` moves by :math:`A\,r_{ij}`, where :math:`r_{ij}` is a
fixed pseudo-random number in :math:`[-1, 1)` computed from the two cell
indices. The bed is the same for any boxes and any number of ranks.

**Bowl.** ``flume.bed.bowl_depth`` (:math:`h_0`) and
``flume.bed.bowl_radius`` (:math:`a`) add a paraboloid to the plane bed. With
:math:`X`, :math:`Y` the distance of a cell centre from the centre of the
domain and :math:`r^2 = X^2 + Y^2`, the bed moves by

.. math::

   -h_0 \left(1 - \frac{r^2}{a^2}\right) .

It is lowered by :math:`h_0` at the centre and crosses the plane at
:math:`r = a`. A negative :math:`h_0` makes a mound.

**Tilted and curved surface.** ``flume.init.surface_slope``
(:math:`\sigma_x`, :math:`\sigma_y`) and ``flume.init.surface_curvature``
(:math:`\kappa`) replace the level surface by

.. math::

   \eta(x, y) = \eta_0 + \sigma_x X + \sigma_y Y + \kappa\, r^2 ,

so :math:`\eta_0` is then the elevation at the centre of the domain. With the
bowl these two give the initial states of Thacker's solutions
(:doc:`Testing`).

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
neighbour. Before the update the fluxes of a cell that would give away more
water than it holds are scaled down (`Wetting and drying`_).

Hydrostatic reconstruction
--------------------------

The bed is a step at each face. The face between a cell with depth
:math:`h` and bed :math:`z_b` and a neighbour with bed :math:`z_b'` sees the
cell through the depth (Audusse et al., 2004)

.. math::

   h^{*} = \max\bigl(0,\; h - \max(0,\; z_b' - z_b)\bigr) .

This is the part of the cell's water that lies above the higher of the two
beds. It is never negative. A value at or below the dry depth
:math:`h_{\mathrm{dry}}` (``flume.dry_depth``) is replaced by zero
(`Wetting and drying`_). Only the difference of the two bed elevations
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

The flux has the two-wave form of Harten, Lax and van Leer (HLL). The form,
the wave speeds and the dry-bed speeds used here are written out for the
shallow-water equations by Satyaprasad et al. (2024, section 3.4). The
tangential momentum is carried by the mass flux from the upwind side, as in
the solver of Delestre et al. (2014, section 4.4), which follows the HLLC
extension of Toro and co-workers.

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
Delestre et al. (2014) pick the side by the sign of
:math:`u_{n,L} + u_{n,R}`; here it is picked by the sign of the mass flux,
so the tangential momentum always moves with the water that crosses the face.

**Limit.** The two speeds are exact when both waves are rarefactions. For a
shock that runs into much shallower water they are below the exact shock
speed: for still water next to a very thin layer the estimate is 1.5 times
the celerity of the deep side where the exact front moves at about 1.7
times. In a check over random pairs of wet states, 13 % had an estimated
speed inside the exact wave fan by more than 1 %. The tests of
:doc:`Testing` pass with these speeds. Speeds built from the depth ratio
across the shock remove the gap; they are not in this version.

Wetting and drying
------------------

Water can reach dry ground and leave it. Four rules keep the depth
non-negative and the velocity bounded where the water is thin. They use one
input, the dry depth :math:`h_{\mathrm{dry}}` (``flume.dry_depth``, default
10\ :sup:`-6` m).

**Dry cells.** A cell with :math:`h \le h_{\mathrm{dry}}` is dry. Its velocity
is zero in the fluxes; elsewhere the velocity is :math:`hu / h`. After each
step the discharges of a dry cell are set to zero. Its water is kept, so mass
is conserved.

**Dry faces.** The depth a cell presents to a face is

.. math::

   h^{*} = \begin{cases}
       h - \max(0,\; z_b' - z_b) & \text{if this is above } h_{\mathrm{dry}} , \\
       0 & \text{otherwise.}
   \end{cases}

The flux and the bed-slope term use this depth. The cut is made on
:math:`h^{*}` and not on :math:`h`: under a level surface the two cells at a
face present the same :math:`h^{*}`, so both are kept or both are cut, and
water at rest stays at rest for any :math:`h_{\mathrm{dry}}`. The two values
are computed from different numbers and agree to round-off only, so a face
whose :math:`h^{*}` is within round-off of :math:`h_{\mathrm{dry}}` could be
cut on one side alone. That would start a flow of
:math:`\tfrac{2}{3} h^{*} \sqrt{g h^{*}}` in water at rest. Where one side is
above :math:`h_{\mathrm{dry}}`, the other is not, and the two differ by no
more than 4 machine epsilons of the larger cell depth, both are cut. The
decision is made once for a face and used by its flux and by the bed-slope
terms of its two cells.

**Outflow limit.** With :math:`F^{h}` and :math:`G^{h}` the mass fluxes, the
water that leaves cell :math:`(i, j)` in a step through its four faces is

.. math::

   o_{ij} = \frac{\Delta t}{\Delta x}\left[\max(F^{h}_{i+1/2,j}, 0)
            + \max(-F^{h}_{i-1/2,j}, 0)\right]
          + \frac{\Delta t}{\Delta y}\left[\max(G^{h}_{i,j+1/2}, 0)
            + \max(-G^{h}_{i,j-1/2}, 0)\right] ,
   \qquad
   \phi_{ij} = \min\left(1,\; \frac{h_{ij}}{o_{ij}}\right) .

The three components of the flux on a face are multiplied by the
:math:`\phi` of the cell the mass flux leaves. A face with no mass flux is
not changed. The new depth is then at least
:math:`h_{ij} - \phi_{ij}\, o_{ij} \ge 0`, for any ``flume.cfl``. Each face
still has one flux, so mass is conserved. The sum :math:`o_{ij}` counts the
water that leaves and not the water that arrives in the same step, so a cell
can be limited that the plain step would have left with water. The bed-slope
term is not scaled; water at rest has no mass flux and is not limited. A depth that the rounding of this
sum leaves below zero by at most 16 machine epsilons of the old depth is set
to zero.

**Why.** Without the limit a wet cell between dry ones can lose more than it
holds: still water of depth :math:`h` with dry cells on its four sides loses
:math:`2 c h / 3` through each face and would end the step with
:math:`h\,(1 - 4C/3)` on square cells, negative for a CFL number
:math:`C > 0.75`. A column of such cells, dry on two sides only, loses
:math:`4C/3` times :math:`\Delta y / (\Delta x + \Delta y)`: up to twice as
much on cells that are long across the flow. Where :math:`\phi < 1` the scheme is no longer the plain
HLLC scheme in that cell for that step.

**Time step.** Dry cells do not enter the largest wave rate (`Time step`_),
so a film that carries a large velocity does not shorten the step.

**Limits.** The scheme is first order, and at a moving shoreline it damps
the motion: see the bowl tests in :doc:`Testing`. A shore that retreats leaves a film on the bed, which stays
until water flows over it again: a cell at or below the dry depth passes no
water on. With
``flume.dry_depth = 0`` a front that advances over a dry bed also carries a
trace of water ahead of it, one cell per step. The
limit of the hydrostatic reconstruction on a steep bed (above) applies to
thin water on a slope in particular.

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

   r = \max_{\text{cells with } h > h_{\mathrm{dry}}}
       \left(\frac{|u| + c}{\Delta x} + \frac{|v| + c}{\Delta y}\right) ,

where :math:`\Delta x` and :math:`\Delta y` are the cell sizes. The maximum
is taken over all boxes and all MPI ranks. The time step is

.. math::

   \Delta t = \frac{C}{r} ,

with the CFL number :math:`C` (``flume.cfl``, :math:`0 < C \le 1`).

The rates of the two directions are added because the scheme updates a cell
with the fluxes of both directions in one step (unsplit). :math:`C \le 1` is
the linear stability bound of the first-order scheme. The depth stays
non-negative when no cell loses more water in a step than it holds. The
wave speeds :math:`s_L` and :math:`s_R` at a face can be larger than the rate
of either cell, so :math:`C \le 1` alone does not give this; the outflow
limit does (`Wetting and drying`_). The state is still checked before every
step and before every plotfile, and the run stops with a message on a depth
that is negative or not finite.

A step that reaches ``flume.stop_time`` is shortened to end on it. A step that
falls short of it by round-off only (4 machine epsilons of the stop time, in
double precision) is lengthened to end on it, so that no step of that size
follows. The time and ``flume.stop_time`` are held in double precision in
every build: a single-precision sum of steps stops growing after about
:math:`2^{24}` steps.

A dry cell (:math:`h \le h_{\mathrm{dry}}`) has no wave and does not enter
the maximum; the velocity is not computed there. If no cell is wet,
:math:`r = 0` and no step follows from the formula. The run then takes one step to
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
* O. Delestre, C. Lucas, P.-A. Ksinant, F. Darboux, C. Laguerre, T.-N.-T. Vo,
  F. James and S. Cordier (2013). SWASHES: a compilation of shallow water
  analytic solutions for hydraulic and environmental studies. *International
  Journal for Numerical Methods in Fluids* 72(3), 269-300.
  https://doi.org/10.1002/fld.3741 (open copy: https://arxiv.org/abs/1110.0288)
* O. Delestre, F. Darboux, F. James, C. Lucas, C. Laguerre and S. Cordier
  (2014). FullSWOF: A free software package for the simulation of shallow
  water flows. arXiv:1401.4125. https://arxiv.org/abs/1401.4125
* D. Satyaprasad, S. N. Kuiry and S. Sundar (2024). A meshless geometric
  conservation weighted least square method for solving the shallow water
  equations. arXiv:2410.05968. https://arxiv.org/abs/2410.05968
