Model
=====

This page describes what the code holds and computes today: the
shallow-water equations over a fixed bed, advanced by a first-order or a second-order
finite-volume scheme with wetting and drying, between walls, outflow sides
or periodic sides, with Manning friction on the bed. Rainfall and inflow
boundaries are not written yet.

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
is the acceleration of gravity (``flume.gravity``). With Manning friction
(`Manning friction`_) the two momentum equations also have the friction
term on their right-hand side; without it the only source term is the bed
slope.

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

Five more features of the bed and the initial state build test cases.

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

**Wave.** ``flume.init.wave_amplitude`` (:math:`A`) adds one sine wave to
the level surface along x or y (``flume.init.wave_direction``), in a domain
that is periodic along it. With :math:`s` the distance of a cell centre from
the low side of the domain over its length, :math:`H = \eta_0 - z_b` and
:math:`c = \sqrt{g h}`,

.. math::

   h = \max\bigl(0,\; H + A \sin 2\pi s\bigr) , \qquad
   u = u_0 + 2\left(c - \sqrt{g \max(0, H)}\right)

along the wave, where :math:`u_0` is the component of
``flume.init.velocity`` along it; the velocity across it is the other
component. Over a flat bed this is a simple wave, whose exact solution is
known until it breaks (:doc:`Testing`).

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
:math:`S` is the bed-slope term, both defined below. This scheme is first order
in space and time, and is the default (``flume.order = 1``). ``flume.order = 2``
selects the scheme of `Second order`_, which is built from the same parts. The fluxes are computed once per face and stored, so the
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
still has one flux, so mass is conserved, apart from the flux through an
outflow side (`Sides of the domain`_). The sum :math:`o_{ij}` counts the
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

**Limits.** The first-order scheme damps the motion at a moving shoreline;
the second-order one damps it less: see the bowl tests in :doc:`Testing`. A shore that retreats leaves a film on the bed, which stays
until water flows over it again: a cell at or below the dry depth passes no
water on. With
``flume.dry_depth = 0`` a front that advances over a dry bed also carries a
trace of water ahead of it, one cell per step. The
limit of the hydrostatic reconstruction on a steep bed (above) applies to
thin water on a slope in particular.

Second order
------------

``flume.order = 2`` replaces the constant state of a cell by a linear one in
each direction, and the single stage by two. The first-order scheme stays
available and is not changed by this option.

**Limited differences.** Along x, with :math:`\eta = h + z_b` the surface
elevation, cell :math:`i` gets a difference of the depth, of the surface and
of the two velocities from its two neighbours (the same along y):

.. math::

   \delta h_i = \mathrm{mm}(h_i - h_{i-1},\; h_{i+1} - h_i) , \qquad
   \delta \eta_i = \mathrm{mm}(\eta_i - \eta_{i-1},\; \eta_{i+1} - \eta_i) ,

and :math:`\delta u_i`, :math:`\delta v_i` alike, with the minmod function

.. math::

   \mathrm{mm}(a, b) = \begin{cases}
       a & \text{if } ab > 0 \text{ and } |a| \le |b| , \\
       b & \text{if } ab > 0 \text{ and } |a| > |b| , \\
       0 & \text{otherwise.}
   \end{cases}

The difference of the surface is computed as
:math:`(h_i - h_{i-1}) + (z_{b,i} - z_{b,i-1})`, so large elevations cost no
digits. The velocity of a dry neighbour counts as zero. A dry cell
(:math:`h_i \le h_{\mathrm{dry}}`) is not reconstructed: all its differences
are zero.

**Sides of a cell.** The cell presents to its low-x face (W) and to its
high-x face (E)

.. math::

   h_{W} = h_i - \tfrac{1}{2}\delta h_i , \qquad
   h_{E} = h_i + \tfrac{1}{2}\delta h_i , \qquad
   z_{W} = z_{b,i} - \tfrac{1}{2}(\delta\eta_i - \delta h_i) , \qquad
   z_{E} = z_{b,i} + \tfrac{1}{2}(\delta\eta_i - \delta h_i) ,

.. math::

   u_{W} = u_i - \frac{h_{E}}{2 h_i}\,\delta u_i , \qquad
   u_{E} = u_i + \frac{h_{W}}{2 h_i}\,\delta u_i ,

and :math:`v_W`, :math:`v_E` alike. The bed on a side is what the surface and
the depth leave: :math:`h_W + z_W = \eta_i - \tfrac{1}{2}\delta\eta_i`. Three
properties follow.

* :math:`|\delta h_i| \le h_i` when the three depths are not negative (the
  minmod value is no larger than :math:`h_i - h_{i-1} \le h_i` when positive,
  and than :math:`h_i - h_{i+1} \le h_i` in magnitude when negative), so
  :math:`h_W, h_E \ge h_i / 2 \ge 0`.
* :math:`\tfrac{1}{2}(h_W u_W + h_E u_E) = h_i u_i`: the mean of the two
  sides is the discharge of the cell. The weights :math:`h_E / 2h_i` and
  :math:`h_W / 2h_i` lie between 1/4 and 3/4.
* Under a level surface, :math:`\delta\eta_i = 0`: the surface is level on
  the two sides whatever the bed and the depth do.

**Faces.** The face between cells :math:`i` and :math:`i+1` applies the
`Hydrostatic reconstruction`_ and the `HLLC flux`_ to the high side of cell
:math:`i` and the low side of cell :math:`i+1`: depths :math:`h_{E,i}` and
:math:`h_{W,i+1}`, beds :math:`z_{E,i}` and :math:`z_{W,i+1}`, and the
velocities of these two sides. This gives :math:`h^{*}_{E}` for cell
:math:`i`, with the same cut at the dry depth.

**Bed-slope term.** Inside the cell the bed goes from :math:`z_W` to
:math:`z_E`, and the hydrostatic reconstruction leaves a step at each face.
The term has one part for each (Audusse et al., 2004):

.. math::

   S^{hu}_{i} = \frac{g}{2\,\Delta x}\left[
       \left(h^{*2}_{E} - h_{E}^{2}\right) - \left(h^{*2}_{W} - h_{W}^{2}\right)
       - (h_W + h_E)(z_E - z_W) \right] .

The first two parts are the steps at the two faces, as at first order. The
last is :math:`-g h\, \partial z_b / \partial x` over the cell, with the mean
of the two depths. With :math:`z_E - z_W = \delta\eta_i - \delta h_i`,
:math:`h_E - h_W = \delta h_i` and :math:`h_W + h_E = 2 h_i`, the terms
without a star are

.. math::

   -\left(h_E^2 - h_W^2\right) - (h_W + h_E)(z_E - z_W)
     = -2 h_i\,\delta h_i - 2 h_i\,(\delta\eta_i - \delta h_i)
     = -2 h_i\,\delta\eta_i ,

and the term is

.. math::

   S^{hu}_{i} = \frac{g}{2\,\Delta x}\left(h^{*2}_{E} - h^{*2}_{W}\right)
     - \frac{g}{\Delta x}\, h_i\,\delta\eta_i .

The first part is the first-order term with the depths of the sides. The
second is the weight of the water on the slope of the surface inside the
cell.

**Why.** For water at rest under a level surface :math:`\delta\eta_i = 0`
and the surface on every side of every cell is the same. The two sides of a
face then present the same :math:`h^{*}`, as at first order, and the momentum
flux :math:`\tfrac{1}{2} g h^{*2}` on the two faces of a cell is cancelled by
the first part of :math:`S`; the second part is zero to round-off. The water
stays at rest to round-off.

**Two stages.** With :math:`\mathcal{E}(U, \Delta t)` one step of the
`Finite-volume step`_ with the fluxes and the bed-slope term above, including
the outflow limit and the dry-cell rule,

.. math::

   U^{(1)} = \mathcal{E}(U^{n}, \Delta t) , \qquad
   U^{(2)} = \mathcal{E}(U^{(1)}, \Delta t) , \qquad
   U^{n+1} = \tfrac{1}{2}\left(U^{n} + U^{(2)}\right) .

This is Heun's method. The step :math:`\Delta t` is computed once, from
:math:`U^{n}`. Each stage keeps the depth non-negative (the outflow limit
uses the cell depth and applies to any flux), and a mean of two non-negative
depths is non-negative. Each stage conserves mass, apart from the flux
through an outflow side, so the mean does. The
discharges of a cell whose mean depth is at or below
:math:`h_{\mathrm{dry}}` are set to zero.

**Limits.** The minmod limiter sets the differences to zero at a maximum or
a minimum of a field, where the scheme is then first order. The error of a
smooth solution with such points falls more slowly than the square of the
cell size in the largest norm (:doc:`Testing` gives the measured rates).
Next to a dry cell the depth difference is limited by the neighbour's depth
and the cell is close to first order. The outflow limit and the cut at the
dry depth are those of the first-order scheme.

Source terms
------------

A source term changes the state in a cell without a flux between cells.
Each source is a class with a name and one step of :math:`\Delta t`
(``Source/FLUME_Source.H``). The run builds the list of sources its inputs
switch on and prints their names at start-up (``FLUME: sources: ...``, or
``none``). Each source acts once per step, after the fluxes and the bed slope
of all the stages of the step, on the state at the end of the step, in the
order of the list. Today the list is empty or holds Manning friction (next
section). Rain, infiltration and drains are to be added as further sources.

The bed-slope term is not one of these sources: it is part of each explicit
stage, where it balances the fluxes of still water exactly (`Hydrostatic
reconstruction`_).

Manning friction
----------------

With the Manning coefficient :math:`n` (``flume.friction.manning``, in
s/m\ :sup:`1/3`) and the discharge :math:`\mathbf{q} = (hu, hv)`, the two
momentum equations gain the term

.. math::

   \mathbf{S}_f = -\frac{g\, n^2\, |\mathbf{q}|\, \mathbf{q}}{h^{7/3}} .

In uniform flow down a bed of slope :math:`S` it balances gravity,
:math:`g h S = g n^2 q^2 / h^{7/3}`, which gives Manning's formula for the
normal discharge, :math:`q_n = h^{5/3} S^{1/2} / n`.

**Implicit update.** After the explicit terms (fluxes and bed slope) have
given the discharge :math:`\mathbf{q}^*` and the depth :math:`h` at the end
of the step, the friction step is backward Euler in the discharge with
that depth held fixed:

.. math::

   \mathbf{q} + \Delta t\, \frac{g\, n^2\, |\mathbf{q}|}{h^{7/3}}\, \mathbf{q}
      = \mathbf{q}^* .

:math:`\mathbf{q}` has the direction of :math:`\mathbf{q}^*`, and its size
:math:`r |\mathbf{q}^*|` solves :math:`a r^2 + r - 1 = 0` with
:math:`a = \Delta t\, g\, n^2 |\mathbf{q}^*| / h^{7/3}`:

.. math::

   r = \frac{2}{1 + \sqrt{1 + 4 a}} .

**Why this form.** :math:`0 < r \le 1`: friction never reverses the flow and
never adds to it, whatever the step, so it needs no limit on the time step.
When friction is stiff (:math:`a \gg 1`, thin water on a steep slope) the
discharge goes to the balance of the explicit terms and friction in one
step. For uniform flow on a plane, where the explicit terms add
:math:`\Delta t\, g h S` per step, the steady state of the update is the
balance :math:`g h S = g n^2 q^2 / h^{7/3}` for any :math:`\Delta t`. The
formula is written without :math:`h^{14/3}`, which can overflow; a depth
whose power 7/3 underflows to zero gets a zero discharge. A cell at or
below :math:`h_{\mathrm{dry}}` is not touched (its discharge is already zero).

**At second order.** With :math:`L` the explicit terms and :math:`F` the
friction step, a step is

.. math::

   U^* = U^n + \Delta t\, L(U^n), \quad
   M = \tfrac{1}{2}\left(U^n + U^* + \Delta t\, L(U^*)\right), \quad
   U^{n+1} = F(M) :

the two explicit stages and their mean as without friction, and then
friction once over the whole step. If friction acted on each stage and the
mean were taken of the damped stages, the mean with :math:`U^n` would remove
only half of a stiff film's excess discharge per step.

The second stage starts from a discharge that friction has not damped. In a
stiff film it is much larger than the discharge the step ends with (about
30 times in the thin-film test), and its fluxes run at a Courant number of
about 4 for a step chosen from the start of the step. The outflow limit
keeps the depth at zero or above. Damping the first stage as well, and
giving the amount back before the mean, was tried; it changed none of the
test results, including a film of uneven depth (:doc:`Testing`), so it is
not done.

**Limits.** Friction is first order in time at both orders: the error of
the transient in the tests halves with the step (:doc:`Testing`). At first
order the bed-slope term of the hydrostatic reconstruction drives uniform
flow on a plane with :math:`g S (h - S \Delta x / 2)` instead of
:math:`g h S` (Delestre et al., 2012): the steady discharge is low by the
fraction :math:`1 - \sqrt{1 - S \Delta x / (2 h)} \approx S \Delta x / (4 h)`,
for :math:`S \Delta x \le h`. A film thinner than the drop of the bed over
one cell (:math:`S \Delta x > h`) has a dry lower face; it is driven with
:math:`g h^2 / (2 \Delta x)`, a fraction :math:`h / (2 S \Delta x)` of
:math:`g h S`, and its steady discharge is the normal discharge times the
square root of that fraction (7 % for 1 mm on a slope of 0.1 on cells of
1 m). The
second-order reconstruction is exact on a plane. The coefficient is the same in every cell; the solver
holds it as a field per cell, which a roughness map read from a raster will
fill.

Sides of the domain
-------------------

Each direction is periodic or not (``flume.periodic``). The state and the
bed carry one layer of ghost cells at first order and two at second order.
Between boxes and across a periodic side a ghost cell holds the value of the
cell it covers.

A side that is not periodic is made of segments (``flume.bc.segments``):
each a run of boundary cells along the side, with one type. A boundary cell
belongs to a segment when its centre is at the first coordinate of
``flume.bc.<name>.range`` or above, and below the second; the run prints
the cells of each segment. The cells no segment covers are walls, so
without ``flume.bc.segments`` every side that is not periodic is a wall.
Two types exist:

**Wall.** The first ghost layer mirrors the first cell inside and the second
layer the second cell: same depth and bed, same discharge along the wall,
and the discharge normal to the wall with its sign changed. The mass flux
through a wall is then exactly zero.

**Outflow.** Every ghost layer takes the boundary cell next to it, the bed
included (zero gradient). Water leaves, or enters, as the flow at the side
carries it. Where the flow leaving is
supercritical the boundary sets nothing the flow inside depends on, and a
dry-bed front leaves the domain as the exact solution says (:doc:`Testing`).
Where it is subcritical, zero gradient is an approximation: a shock that
leaves sends a weak wave back. In the dam-break test the wave changes the
depth by 0.2 % and the velocity by 0.3 % at first order, 0.7 % and 2 % at
second order. The wave does not get smaller as the cells get smaller, and
it grows with the Froude number of the flow leaving. No discharge or level
is prescribed at an outflow side.

Still water next to an outflow side stays still to round-off over a single
step, since the surface of the ghost cells is that of the boundary cell. It
is not stable where the bed of a boundary cell is lower than that of the
cell inside it. With the hydrostatic reconstruction the face to the inside
then sees the boundary cell :math:`h^*` deep, :math:`h^* < h`, while the
face to the ghost cell sees it :math:`h` deep and carries all of its
discharge. A disturbance of still water grows from round-off. In the lake
at rest over the rough bed (``lake_at_rest_segments``, :doc:`Testing`) the
largest discharge grows from 10\ :sup:`-13` to 10\ :sup:`-3` m\ :sup:`2`/s
in 8000 steps at both orders; with islands it reaches 30 m\ :sup:`2`/s
within 4000 steps. A 1-D model of the same scheme gives growth factors per step up
to 1.09 over rough beds and 1.002 to 1.007 over uniform slopes of 0.001 to
0.1. Manning friction acts on the square of the discharge and does not damp
it. Runs whose water comes to rest next to an outflow side over a bed that
steps down to the side should not run for thousands of steps in that state.
Scaling the discharge of the ghost cell by :math:`h^*/h` removes the growth,
but makes the side act as a weir the height of the last step of the bed:
in the same 1-D model, flow near critical on a slope of 0.01 leaves the
side 48 % too deep. That scaling is not used.

Every box whose ghost layers reach outside a side is filled, also a box that
ends one cell short of it. The flux of the outflow limit (`Wetting and
drying`_) leaves a cell through an outflow side as through any face; water
entering from a ghost cell is not limited. The names ``discharge`` and
``stage`` are kept for a prescribed inflow discharge and a prescribed water
level from a time series, which later versions add; today they stop the run
with a message.

By default a sloping bed is taken as it is at a periodic side: the bed
repeats, and has a step of the slope times the length of the domain there.
With ``flume.bed.continue_slope = 1`` a ghost cell beyond a periodic side
takes the bed of the cell it covers plus the slope times the distance
between the two, so the plane of the bed goes on: an endless inclined
channel. The depth and the discharges still repeat. This is how the
friction tests drive uniform flow (:doc:`Testing`); a lake at rest over such
a bed is not at rest, since its surface does not repeat. Only the plane is
continued: a bowl (``flume.bed.bowl_depth``) still repeats, with a step at
the periodic side where it does not reach zero there.

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

with the CFL number :math:`C` (``flume.cfl``, :math:`0 < C \le 1`, and
:math:`C \le 0.5` at second order).

The rates of the two directions are added because the scheme updates a cell
with the fluxes of both directions in one step (unsplit). :math:`C \le 1` is
the linear stability bound of the first-order scheme. The second-order
scheme uses the same step for its two stages, with :math:`C \le 0.5`: a
side of a cell holds up to 1.5 times the depth of the cell, so above half
the first-order bound the water that leaves through one face can exceed
what the cell holds, and the outflow limit, which is first order where it
acts, becomes part of the scheme. On cells much longer across the wave than
along it, where one direction sets the whole step, the simple wave of the
tests has an error of 0.36 % of its amplitude at :math:`C = 0.5` and of
1.53 % at :math:`C = 1` on 100 cells, and its rate of convergence falls from
1.9 to 1.25 (:doc:`Testing`). The depth stays
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
