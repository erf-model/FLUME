Bed, buildings and initial state
================================

Plane bed
---------

Without ``flume.bed.file`` the bed is a plane through the elevation
:math:`z_0` (``flume.bed.elevation``) at the lower corner :math:`(x_0, y_0)`
of the domain, with slopes :math:`s_x`, :math:`s_y` (``flume.bed.slope``):

.. math::

   z_b(x, y) = z_0 + s_x\,(x - x_0) + s_y\,(y - y_0) ,

evaluated at the cell centres. Two options for test cases change it:

* **Bowl** (``flume.bed.bowl_depth`` :math:`h_0`, ``flume.bed.bowl_radius``
  :math:`a`): with :math:`r` the distance from the centre of the domain, the
  bed moves by :math:`-h_0\,(1 - r^2/a^2)`. A negative :math:`h_0` makes a
  mound.
* **Rough bed** (``flume.bed.noise`` :math:`A`): the bed of cell
  :math:`(i, j)` moves by :math:`A\,r_{ij}`, with :math:`r_{ij}` a fixed
  pseudo-random number in :math:`[-1, 1)` computed from the cell indices. It
  is the same for any boxes and any number of ranks.

.. _sec:BedRaster:

Bed from a raster
-----------------

``flume.bed.file`` gives the whole bed as an ESRI ASCII grid
(:ref:`sec:RasterFiles`). Let the raster have :math:`N_x \times N_y` cells of
size :math:`\delta_x \times \delta_y`, lower left corner :math:`(X_0, Y_0)`
and values :math:`Z_{ab}` (column :math:`a` from the west, row :math:`b` from
the south). At a cell centre :math:`(x, y)` the position in raster cells,
clamped to the raster cell centres, is

.. math::

   f_x = \min\Bigl(\max\Bigl(\frac{x - X_0}{\delta_x} - \frac{1}{2},\; 0\Bigr),\;
         N_x - 1\Bigr) ,

and :math:`f_y` alike. With :math:`a = \lfloor f_x \rfloor`,
:math:`w_x = f_x - a`, and :math:`b`, :math:`w_y` likewise, the bed is the
bilinear interpolant

.. math::

   z_b = (1 - w_x)(1 - w_y)\, Z_{a,b} + w_x (1 - w_y)\, Z_{a+1,b}
       + (1 - w_x)\, w_y\, Z_{a,b+1} + w_x w_y\, Z_{a+1,b+1} .

The interpolant is continuous, so a coarse raster on a fine grid leaves no
steps in the bed, and it reproduces a plane exactly. Within half a raster
cell of the raster's edges the edge values are held.

.. _sec:Buildings:

Buildings
---------

``flume.buildings.file`` gives each cell the value :math:`H` of the raster
cell that holds its centre, with no interpolation. Where :math:`H > 0` the
bed is raised by :math:`H`, or by ``flume.buildings.height`` when the file is
a mask. A building is a block of raised bed: the hydrostatic reconstruction
(:ref:`sec:HydrostaticReconstruction`) lets water onto it only when the
surface beside it is above its top, and below that its face acts as a wall.

Limitations:

* No porosity. A building narrower than a cell is lost or widened to a cell.
* The only force on the water at a building face is the hydrostatic pressure
  of the cell's own depth. A wall boundary adds the pressure of the reflected
  flow; the two agree at rest and differ in moving water.

.. _sec:InitialState:

Initial state
-------------

The initial surface is level at :math:`\eta_0`
(``flume.init.surface_elevation``) and the water moves with one velocity
:math:`(u_0, v_0)` (``flume.init.velocity``):

.. math::

   h = \max(\eta_0 - z_b,\; 0), \qquad hu = h\,u_0, \qquad hv = h\,v_0 .

.. list-table:: Options
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Initial state
   * - Dam
     - Cells whose centre is at or above ``flume.init.dam_position`` along
       ``flume.init.dam_direction`` take
       ``flume.init.dam_surface_elevation`` and, if given,
       ``flume.init.dam_velocity``.
   * - Tilted and curved surface
     - :math:`\eta = \eta_0 + \sigma_x X + \sigma_y Y + \kappa r^2`, with
       :math:`X`, :math:`Y`, :math:`r` measured from the centre of the domain
       (``flume.init.surface_slope``, ``flume.init.surface_curvature``).
       With the bowl these give Thacker's solutions (:doc:`../Verification`).
   * - Wave
     - One sine wave of amplitude :math:`A` (``flume.init.wave_amplitude``)
       along a periodic direction. With :math:`s` the fraction of the domain
       length, :math:`H = \eta_0 - z_b` and :math:`c = \sqrt{g h}`:
       :math:`h = \max(0, H + A \sin 2\pi s)`,
       :math:`u = u_0 + 2\,(c - \sqrt{g H})`. Over a flat bed this is a
       simple wave.
