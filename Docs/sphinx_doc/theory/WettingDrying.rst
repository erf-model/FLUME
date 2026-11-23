Wetting and drying
==================

Water can reach dry ground and leave it. Four rules keep the depth
non-negative and the velocity bounded in thin water. They use one input, the
dry depth :math:`h_{\mathrm{dry}}` (``flume.dry_depth``, default
10\ :sup:`-6` m).

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Rule
     - Statement
   * - Dry cell
     - A cell with :math:`h \le h_{\mathrm{dry}}` has zero velocity in the
       fluxes, and its discharges are set to zero after each step. Its
       water is kept.
   * - Dry face
     - A reconstructed depth :math:`h^{*} \le h_{\mathrm{dry}}`
       (:ref:`sec:HydrostaticReconstruction`) is replaced by zero. The cut
       is decided once per face and used by its flux and by the bed-slope
       terms of its two cells.
   * - Outflow limit
     - The fluxes leaving a cell are scaled so the cell loses no more water
       than it holds (below).
   * - Time step
     - Dry cells do not enter the CFL condition (:doc:`TimeStep`).

The dry-face cut acts on :math:`h^{*}`, which is the same on both sides of a
face under a level surface, so water at rest stays at rest for any
:math:`h_{\mathrm{dry}}`. Where the two values straddle
:math:`h_{\mathrm{dry}}` by round-off only (4 machine epsilons of the larger
cell depth), both are cut.

Outflow limit
-------------

With :math:`F^{h}`, :math:`G^{h}` the mass fluxes, the water that leaves
cell :math:`(i, j)` in a step is

.. math::

   o_{ij} = \frac{\Delta t}{\Delta x}\left[\max(F^{h}_{i+1/2,j}, 0)
            + \max(-F^{h}_{i-1/2,j}, 0)\right]
          + \frac{\Delta t}{\Delta y}\left[\max(G^{h}_{i,j+1/2}, 0)
            + \max(-G^{h}_{i,j-1/2}, 0)\right] ,
   \qquad
   \phi_{ij} = \min\left(1,\; \frac{h_{ij}}{o_{ij}}\right) .

The three components of the flux on a face are multiplied by the
:math:`\phi` of the cell the mass flux leaves. The new depth is then at
least :math:`h_{ij} - \phi_{ij}\, o_{ij} \ge 0` for any ``flume.cfl``. Each
face still has one flux, so mass is conserved. The bed-slope term is not
scaled, and water at rest is not limited.

Without the limit, still water of depth :math:`h` with dry cells on its four
sides would end a step with :math:`h\,(1 - 4C/3)` on square cells: negative
for a CFL number :math:`C > 0.75`.

Limitations
-----------

* Where :math:`\phi < 1` the scheme is no longer the plain HLLC scheme in
  that cell for that step.
* The first-order scheme damps the motion of a moving shoreline; the
  second-order scheme damps it much less (:doc:`../Verification`).
* A retreating shore leaves a film at or below :math:`h_{\mathrm{dry}}` on
  the bed, which stays until water flows over it again.
* With ``flume.dry_depth = 0`` a front on a dry bed carries a trace of water
  ahead of it, one cell per step.
