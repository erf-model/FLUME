Governing equations
===================

FLUME solves the shallow-water equations for the depth :math:`h` and the
discharges :math:`hu` and :math:`hv` over a fixed bed of elevation
:math:`z_b`:

.. math::

   \frac{\partial h}{\partial t} + \frac{\partial (hu)}{\partial x}
      + \frac{\partial (hv)}{\partial y} &= R , \\
   \frac{\partial (hu)}{\partial t}
      + \frac{\partial}{\partial x}\left(hu^2 + \tfrac{1}{2} g h^2\right)
      + \frac{\partial (huv)}{\partial y}
      &= -g h \frac{\partial z_b}{\partial x} + S_{f,x} , \\
   \frac{\partial (hv)}{\partial t} + \frac{\partial (huv)}{\partial x}
      + \frac{\partial}{\partial y}\left(hv^2 + \tfrac{1}{2} g h^2\right)
      &= -g h \frac{\partial z_b}{\partial y} + S_{f,y} .

.. list-table:: Symbols
   :header-rows: 1
   :widths: 18 18 64

   * - Symbol
     - Unit
     - Meaning
   * - :math:`h`
     - m
     - Water depth; zero in a dry cell.
   * - :math:`hu`, :math:`hv`
     - m\ :sup:`2`/s
     - Discharge per unit width in x and y.
   * - :math:`u`, :math:`v`
     - m/s
     - Depth-averaged velocity.
   * - :math:`z_b`
     - m
     - Bed elevation.
   * - :math:`g`
     - m/s\ :sup:`2`
     - Acceleration of gravity (``flume.gravity``).
   * - :math:`R`
     - m/s
     - Rain rate (:ref:`sec:Rain`).
   * - :math:`\mathbf{S}_f`
     - m\ :sup:`2`/s\ :sup:`2`
     - Manning friction (:ref:`sec:Friction`).

State
-----

The three conserved quantities and the bed are held at the cell centres of
one level of cells in the horizontal plane, as AMReX ``MultiFab`` objects on
one ``BoxArray``. The boxes are one cell thick.

.. _sec:Datum:

Vertical datum
--------------

The bed elevation and the water-surface elevation :math:`\eta` are heights
above one vertical datum. The inputs give the bed and the initial surface as
elevations; the depth is never an input. Afterwards

.. math::

   \eta = z_b + h ,

and :math:`\eta = z_b` in a dry cell. The scheme uses only differences of
elevations, so a large datum offset costs no digits in double precision
(:ref:`sec:Precision`).
