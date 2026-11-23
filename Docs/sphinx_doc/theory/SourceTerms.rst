Source terms
============

A source term changes the state in a cell without a flux between cells. Each
source is a class with a name and one step of :math:`\Delta t`
(``Source/FLUME_Source.H``). The active sources are printed at start-up and
act once per step, after all the explicit stages, in this order: rain, then
Manning friction. They act on a dry domain too.

The bed-slope term is not one of these sources: it is part of each explicit
stage (:ref:`sec:HydrostaticReconstruction`).

.. _sec:Rain:

Rain
----

Rain is the same in every cell, with a rate :math:`r(t)` given by a table of
times :math:`t_0 < \dots < t_N` and rates :math:`r_0, \dots, r_N \ge 0`.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``flume.rain.type``
     - Rate
   * - ``constant``
     - :math:`r(t) = r_0` (``flume.rain.rate``).
   * - ``hyetograph``, ``step``
     - :math:`r_k` from :math:`t_k` to :math:`t_{k+1}`
       (``flume.rain.file``).
   * - ``hyetograph``, ``linear``
     - Linear in time between the entries.

Before :math:`t_0` the rate is :math:`r_0`; from :math:`t_N` on it is
:math:`r_N`.

**Step.** The step from :math:`t^n` to :math:`t^{n+1}` adds to every depth
the integral of the rate,

.. math::

   \Delta h = I(t^{n+1}) - I(t^{n}) , \qquad
   I(t) = \int_{t_0}^{t} r(\tau)\, d\tau ,

evaluated in closed form from the table. The discharges are not changed:
rain brings no horizontal momentum. Rain falls on every cell, wet or dry.

**Properties.**

* The depth added over a run is :math:`I(T) - I(0)` to round-off, whatever
  the steps: a change of rate inside a step is integrated exactly.
* Each cell keeps the part of the rain that the addition rounded away and
  adds it with the next step (a compensated sum). In single precision the
  rain of one step can be below the rounding of the depth; the sum keeps it.
* Rain is first order in time at both orders of the scheme.

The solver asks a forcing provider (``Source/FLUME_Forcing.H``) for the depth
of rain over an interval and for the largest rate over an interval
(:ref:`sec:RainStep`). Rain that varies in space comes from two further
providers (:ref:`sec:RainSpace`).

.. _sec:RainSpace:

Rain that varies in space
-------------------------

Both providers give each cell the integral of the rate over the cell and
the step, divided by the cell area. The compensated sum and the step bound
are those of :ref:`sec:Rain`.

**Gridded rain** (``gridded``). Entry :math:`k` of the index file is a time
:math:`T_k` and a grid of rates :math:`g_{k,m}` (m/s) on cells :math:`A_m`.
On a solver cell :math:`\Omega_{ij}` of area :math:`|\Omega_{ij}|`,

.. math::

   M_{k,ij} = \frac{1}{|\Omega_{ij}|} \sum_m |\Omega_{ij} \cap A_m|\, g_{k,m} ,
   \qquad
   \Delta h_{ij} = \sum_k w_k\, M_{k,ij} , \qquad
   w_k = W_k(t^{n+1}) - W_k(t^{n}) ,

where :math:`W_k(t)` is the integral up to :math:`t` of the weight of entry
:math:`k` in the rate: 1 from :math:`T_k` to :math:`T_{k+1}` (``step``), or
the hat function that is 1 at :math:`T_k` and 0 at :math:`T_{k\pm1}`
(``linear``); the first weight is 1 before :math:`T_0` and the last after
:math:`T_N`.

* Cell edges are formed in grid units by one formula, so neighbouring cells
  share each edge to the last bit and the rain on the domain is the
  integral of the grid over it, to round-off.
* An edge within 16 epsilons of a whole grid unit is set to it, so a grid
  cell beyond it is not read through a sliver of round-off.
* The weights are differences of :math:`W_k`, so the weights of all steps
  add up to those of the run: the volume delivered is the integral of the
  input, whatever the steps.
* The step bound is the largest mapped rate of the grids an interval uses.

**Moving cell** (``moving_cell``). With peak rate :math:`P`, standard
deviation :math:`\sigma`, centre :math:`\mathbf{c}_0` at :math:`t = 0` and
velocity :math:`\mathbf{v}`,

.. math::

   r(\mathbf{x}, t) = P \exp\left(-\frac{|\mathbf{x} - \mathbf{c}(t)|^2}{2\sigma^2}\right),
   \qquad
   \mathbf{c}(t) = \mathbf{c}_0 + \mathbf{v}\, t .

The mean over a cell :math:`[x_a, x_b] \times [y_a, y_b]` at one time is
exact:

.. math::

   \frac{P}{|\Omega_{ij}|}\, \frac{\pi \sigma^2}{2}\,
   \left[\mathrm{erf}\,\frac{x_b - c_x}{\sqrt{2}\sigma} - \mathrm{erf}\,\frac{x_a - c_x}{\sqrt{2}\sigma}\right]
   \left[\mathrm{erf}\,\frac{y_b - c_y}{\sqrt{2}\sigma} - \mathrm{erf}\,\frac{y_a - c_y}{\sqrt{2}\sigma}\right] .

* Each difference is formed from :math:`\mathrm{erfc}` on the side away
  from the centre, so a cell far from the centre keeps its digits.
* In time: the 4-point Gauss-Legendre rule on
  :math:`n = \max(1, \lceil 4|\mathbf{v}|\,\Delta t/\sigma \rceil)`
  sub-intervals, so the centre moves at most :math:`\sigma/4` in each. A
  cell at rest is exact with one.
* A cell of peak 0 adds nothing and is not integrated. A step that needs
  more than :math:`10^4` sub-intervals stops the run.
* Centre positions are formed in double precision relative to
  ``flume.prob_lo``. The step bound is :math:`P`.

.. _sec:Friction:

Manning friction
----------------

With the Manning coefficient :math:`n` (``flume.friction.manning`` or
``flume.friction.manning_file``) and the discharge
:math:`\mathbf{q} = (hu, hv)`:

.. math::

   \mathbf{S}_f = -\frac{g\, n^2\, |\mathbf{q}|\, \mathbf{q}}{h^{7/3}} .

In uniform flow down a slope :math:`S` it balances gravity, which gives
Manning's normal discharge :math:`q_n = h^{5/3} S^{1/2} / n`.

**Implicit update.** With :math:`\mathbf{q}^*` and :math:`h` from the
explicit terms, the friction step is backward Euler in the discharge at
fixed depth:

.. math::

   \mathbf{q} + \Delta t\, \frac{g\, n^2\, |\mathbf{q}|}{h^{7/3}}\, \mathbf{q}
      = \mathbf{q}^* ,
   \qquad
   \mathbf{q} = r\, \mathbf{q}^* , \qquad
   r = \frac{2}{1 + \sqrt{1 + 4 a}} , \qquad
   a = \frac{\Delta t\, g\, n^2 |\mathbf{q}^*|}{h^{7/3}} .

Since :math:`0 < r \le 1`, friction never reverses the flow and sets no
limit on the time step. For uniform flow on a plane the steady state of the
update is Manning's balance for any :math:`\Delta t`. A cell at or below
:math:`h_{\mathrm{dry}}` is not touched.

**At second order** friction acts once, on the mean of the two explicit
stages:

.. math::

   U^* = U^n + \Delta t\, L(U^n), \quad
   M = \tfrac{1}{2}\left(U^n + U^* + \Delta t\, L(U^*)\right), \quad
   U^{n+1} = F(M) ,

with :math:`L` the explicit terms and :math:`F` the friction step.

**Limitations.**

* Friction is first order in time at both orders.
* At first order the bed-slope term drives uniform flow on a plane with
  :math:`g S (h - S \Delta x / 2)` instead of :math:`g h S` (Delestre et
  al., 2012). The steady discharge is low by the fraction
  :math:`1 - \sqrt{1 - S \Delta x / (2 h)} \approx S \Delta x / (4 h)`.
* A film thinner than the drop of the bed over one cell
  (:math:`S \Delta x > h`) is driven with a fraction
  :math:`h / (2 S \Delta x)` of :math:`g h S` at first order. The
  second-order reconstruction is exact on a plane.
* In a stiff film at second order the second stage starts from a discharge
  that friction has not damped; the outflow limit keeps the depth
  non-negative.
