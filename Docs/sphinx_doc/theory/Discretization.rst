Discretization
==============

Finite-volume step
------------------

With :math:`U = (h, hu, hv)`, a first-order step (``flume.order = 1``, the
default) updates cell :math:`(i, j)` in one stage, without splitting the
directions:

.. math::

   U_{ij}^{n+1} = U_{ij}^{n}
     - \frac{\Delta t}{\Delta x}\left(F_{i+1/2,j} - F_{i-1/2,j}\right)
     - \frac{\Delta t}{\Delta y}\left(G_{i,j+1/2} - G_{i,j-1/2}\right)
     + \Delta t\, S_{ij} .

:math:`F` and :math:`G` are the face fluxes and :math:`S` the bed-slope
term. Each face has one flux, used by both of its cells. Rain and friction
act after the step (:doc:`SourceTerms`).

.. _sec:HydrostaticReconstruction:

Hydrostatic reconstruction
--------------------------

The bed is a step at each face. A cell with depth :math:`h` and bed
:math:`z_b` presents to a neighbour with bed :math:`z_b'` the depth
(Audusse et al., 2004)

.. math::

   h^{*} = \max\bigl(0,\; h - \max(0,\; z_b' - z_b)\bigr) ,

the part of its water above the higher of the two beds. The flux on the face
is computed from :math:`h^{*}` on each side and the velocities of the two
cells. With :math:`h^{*}_{E}`, :math:`h^{*}_{W}`, :math:`h^{*}_{N}`,
:math:`h^{*}_{S}` the depths a cell presents to its four faces, the bed-slope
term is

.. math::

   S^{hu}_{ij} = \frac{g}{2\,\Delta x}\left(h^{*2}_{E} - h^{*2}_{W}\right) ,
   \qquad
   S^{hv}_{ij} = \frac{g}{2\,\Delta y}\left(h^{*2}_{N} - h^{*2}_{S}\right) .

**Well balanced.** Under a level surface the two cells at a face present the
same :math:`h^{*}`: the mass flux is zero and :math:`S` cancels the
difference of the momentum fluxes. Water at rest stays at rest to round-off
over any bed.

**Limitation.** Where the bed rises at a face by at least the depth of the
lower cell, :math:`h^{*} = 0` and the face sees less of the bed step than
there is. Shallow water moving down a steep bed is then driven too weakly:
the velocity is low and the depth high (Delestre et al., 2012). The effect
falls with the cell size. Water at rest is not affected.

.. _sec:HLLC:

HLLC flux
---------

At a face, let L and R be the low and high sides, :math:`u_n` and
:math:`u_t` the velocities normal to and along the face, and
:math:`c = \sqrt{g h^{*}}`. The wave speeds are estimated from the middle
state of two rarefactions (Satyaprasad et al., 2024, section 3.4):

.. math::

   s_L = \min(u_{n,L} - c_L,\; u^{\star} - c^{\star}) , \qquad
   s_R = \max(u_{n,R} + c_R,\; u^{\star} + c^{\star}) ,

.. math::

   u^{\star} = \tfrac{1}{2}(u_{n,L} + u_{n,R}) + c_L - c_R , \qquad
   c^{\star} = \max\bigl(0,\; \tfrac{1}{2}(c_L + c_R)
               + \tfrac{1}{4}(u_{n,L} - u_{n,R})\bigr) .

With one side dry the speeds are those of a front on a dry bed:
:math:`s_L = u_{n,L} - c_L`, :math:`s_R = u_{n,L} + 2 c_L` with the high
side dry, and :math:`s_L = u_{n,R} - 2 c_R`, :math:`s_R = u_{n,R} + c_R`
with the low side dry. A face dry on both sides has no flux.

With :math:`q = (h^{*},\; h^{*} u_n)` and
:math:`f = (h^{*} u_n,\; h^{*} u_n^2 + \tfrac{1}{2} g h^{*2})`:

.. math::

   f^{\mathrm{HLL}} =
   \begin{cases}
     f_L , & s_L \ge 0 , \\
     \dfrac{s_R f_L - s_L f_R + s_L s_R (q_R - q_L)}{s_R - s_L} ,
           & s_L < 0 < s_R , \\
     f_R , & s_R \le 0 .
   \end{cases}

The momentum along the face is carried by the mass flux from its upwind side
(Delestre et al., 2014, section 4.4):

.. math::

   f_{t} = f^{\mathrm{HLL}}_{h}\; \times
   \begin{cases} u_{t,L} , & f^{\mathrm{HLL}}_{h} \ge 0 , \\
                  u_{t,R} , & f^{\mathrm{HLL}}_{h} < 0 . \end{cases}

For the shallow-water equations the middle wave of the HLLC flux moves with
the sign of the HLL mass flux, so this is the HLLC flux written without the
speed of that wave.

**Limitation.** The speed estimates are exact for two rarefactions. For a
shock into much shallower water they are below the exact shock speed: 1.5
times the celerity of the deep side for still water next to a very thin
layer, where the exact front moves at about 1.7 times.

.. _sec:SecondOrder:

Second order
------------

``flume.order = 2`` replaces the constant state of a cell by a linear one in
each direction and the single stage by two.

**Limited differences.** Along x (and y alike), with
:math:`\eta = h + z_b`:

.. math::

   \delta h_i = \mathrm{mm}(h_i - h_{i-1},\; h_{i+1} - h_i) , \qquad
   \delta \eta_i = \mathrm{mm}(\eta_i - \eta_{i-1},\; \eta_{i+1} - \eta_i) ,

and :math:`\delta u_i`, :math:`\delta v_i` alike, where
:math:`\mathrm{mm}(a, b)` is the minmod function: the smaller of :math:`a`
and :math:`b` in magnitude when they have the same sign, zero otherwise. A
dry cell is not reconstructed, and the velocity of a dry neighbour counts as
zero.

**Sides of a cell.** The cell presents to its low-x face (W) and high-x face
(E)

.. math::

   h_{W,E} = h_i \mp \tfrac{1}{2}\delta h_i , \qquad
   z_{W,E} = z_{b,i} \mp \tfrac{1}{2}(\delta\eta_i - \delta h_i) , \qquad
   u_{W} = u_i - \frac{h_{E}}{2 h_i}\,\delta u_i , \qquad
   u_{E} = u_i + \frac{h_{W}}{2 h_i}\,\delta u_i ,

and :math:`v_W`, :math:`v_E` alike. The side depths are at least
:math:`h_i / 2`, the mean of the two side discharges is :math:`h_i u_i`, and
under a level surface the surface is level on both sides.

**Faces.** The hydrostatic reconstruction and the HLLC flux are applied to
the E side of cell :math:`i` and the W side of cell :math:`i+1`.

**Bed-slope term.** (Audusse et al., 2004)

.. math::

   S^{hu}_{i} = \frac{g}{2\,\Delta x}\left(h^{*2}_{E} - h^{*2}_{W}\right)
     - \frac{g}{\Delta x}\, h_i\,\delta\eta_i .

The first part is the first-order term with the depths of the sides; the
second is the weight of the water on the slope of the surface inside the
cell. Under a level surface :math:`\delta\eta_i = 0` and the scheme is well
balanced as at first order.

**Two stages.** With :math:`\mathcal{E}(U, \Delta t)` one first-order step
with the fluxes and bed-slope term above (Heun's method):

.. math::

   U^{(1)} = \mathcal{E}(U^{n}, \Delta t) , \qquad
   U^{(2)} = \mathcal{E}(U^{(1)}, \Delta t) , \qquad
   U^{n+1} = \tfrac{1}{2}\left(U^{n} + U^{(2)}\right) .

:math:`\Delta t` is computed once, from :math:`U^{n}`. Each stage keeps the
depth non-negative and conserves mass, so the mean does.

**Limitations.** The minmod limiter is first order at a maximum or minimum
of a field, so the error of a smooth solution falls at a rate below 2 in the
largest norm (1.9 measured, :doc:`../Verification`). Next to a dry cell the
scheme is close to first order.
