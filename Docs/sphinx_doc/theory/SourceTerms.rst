Source terms
============

A source term changes the state in a cell without a flux between cells. Each
source is a class with a name and one step of :math:`\Delta t`
(``Source/FLUME_Source.H``). The active sources are printed at start-up and
act once per step, after all the explicit stages, in this order: rain,
inflows inside the domain, infiltration, storm-drain inlets, then Manning
friction. They act on a dry domain too.

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
  A negative carry (rain rounded up) waits while the cell holds less than
  it, so rain never takes the depth below 0, also in a cell that another
  source (infiltration) has just emptied.
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

.. _sec:InflowSource:

Inflows inside the domain
-------------------------

An inflow (``flume.inflows.*``) adds the discharge :math:`Q(t)`,
m\ :sup:`3`/s, of a hydrograph file to the cell that holds a point, or to
the cells a line crosses: a river that enters away from the domain's sides,
a culvert outlet, a pumped discharge or a manhole that overflows. The file
and its interpolation are those of a discharge segment
(:ref:`sec:Inflow`).

**Cells of a line.** The line from :math:`(x_0, y_0)` to :math:`(x_1, y_1)`
is cut where it crosses the lines of the grid; each piece goes to the cell
that holds its middle. Cell :math:`c` takes the part

.. math::

   f_c = \frac{L_c}{L} , \qquad \sum_c f_c = 1 ,

with :math:`L_c` the length of the line inside it and :math:`L` the whole
length. A piece along a cell face, to rounding, goes to the cell above it
or to its right (as a point on a face does), and one on the upper edge of
the domain to the last cell. A piece shorter than :math:`10^{-12} L`, the
rounding where the line passes through a corner of the grid, is dropped,
so a cell the line only touches at a corner is not one of its cells. A point is one cell with :math:`f_c = 1`, the cell of a gauge
at that location (:ref:`sec:Gauges`).

**Step.** The step from :math:`t^n` to :math:`t^{n+1}` adds to the depth
of cell :math:`c`

.. math::

   \Delta h_c = \frac{f_c \left(I(t^{n+1}) - I(t^n)\right)}{\Delta x\, \Delta y} ,
   \qquad I(t) = \int_{t_0}^{t} Q(\tau)\, d\tau ,

with :math:`I` in closed form from the table, as for the rain. The
discharges are not changed: the water enters with no horizontal momentum,
and the bed and the pressure of the water around it set it moving.

**Properties.**

* The volume added over a run is :math:`I(T) - I(0)` to round-off, whatever
  the steps; the balance column ``inflow_<name>`` holds it
  (:doc:`WaterBalance`).
* Each cell keeps the part of the depth that the addition rounded away, as
  for the rain (a compensated sum), so a small discharge into deep water is
  not lost.
* The step bound of the rain covers the inflows (:ref:`sec:RainStep`).
* Inflows act after the rain and before the infiltration, so the soil and
  the drains see their water in the same step.

**Limitations.**

* No momentum: a jet that enters fast (a pipe outlet into a pool) has its
  speed set by the flow around it, not by the pipe. A river that enters
  through a side with its own momentum is a discharge segment
  (:ref:`sec:Inflow`). For the same reason a line across a stream that
  flows slows it there: the new water takes a share of the stream's
  momentum, as a drag of :math:`Q u / A` per unit area.
* The discharge is spread over the cells by length only, not by depth or
  conveyance.
* A cell of an inflow cannot be a building cell (``flume.buildings.file``);
  the run stops.
* The source is first order in time at both orders of the scheme, like the
  other sources.

.. _sec:Infiltration:

Infiltration
------------

Infiltration takes water from the surface into the soil
(``flume.infiltration.type``). Each step a loss model gives a depth
:math:`\Delta F` for each cell, which the cell loses, never more than the
water :math:`W` it holds after the rain of the step. The water that leaves
takes its share of the discharge: :math:`\mathbf{q}` is scaled by
:math:`(h - \Delta F)/h`, so the velocity does not change, and a cell left at
or below :math:`h_{\mathrm{dry}}` holds no discharge. The cumulative loss
:math:`F` of each cell is written to the plotfiles as ``infiltrated_depth``.
It adds the depth actually removed, :math:`h - h_{\mathrm{new}}` as stored
(the subtraction rounds), in a compensated sum like the rain's, so the depth
plus :math:`F` of each cell stays the water it got, to a few epsilons over
any number of steps. So the water balances: water on the surface plus :math:`F`, summed over
the cells, is the initial water plus the rain, less what leaves through the
sides. Building cells and cells with :math:`K_s = 0` or :math:`CN = 100`
infiltrate nothing.

**Green-Ampt** (``green_ampt``), as in Rossman and Huber (2016, section
4.4). With the saturated hydraulic conductivity :math:`K_s`, the suction
head at the wetting front :math:`\psi_s`, the initial moisture deficit
:math:`\theta_d` and :math:`M = \psi_s \theta_d`, the infiltration rate of a
saturated surface is

.. math::

   f = K_s \left(1 + \frac{M}{F}\right) .

With :math:`i = W/\Delta t`, the rate at which water is available in the
step:

* :math:`i \le K_s`: the cell loses :math:`W`.
* Otherwise the surface is saturated once :math:`F` reaches

  .. math::

     F_s = \frac{K_s M}{i - K_s}

  (their eq. 4-29). Below :math:`F_s` the cell loses :math:`W` if
  :math:`F + W < F_s`. Else :math:`F` rises to :math:`F_s` at the rate
  :math:`i`, after :math:`t_s = (F_s - F)/i`, and follows the integrated
  equation for the rest of the step (eq. 4-32):

  .. math::

     K_s \tau = (F_2 - F_1) - M \ln\frac{F_2 + M}{F_1 + M} ,

  with :math:`\tau = \Delta t - t_s` (or :math:`\Delta t` when already
  saturated), :math:`F_1` the cumulative loss at the start of
  :math:`\tau` and :math:`F_2` at its end.

**Why the integrated equation.** It is exact for any :math:`\tau`, so the
ponded infiltration does not depend on the step. :math:`F_2` is found by
Newton's method from the upper bound
:math:`F_1 + K_s\tau + \sqrt{2 M K_s \tau}`, in double precision, with a
series for :math:`u - \ln(1 + u)` when :math:`u = (F_2 - F_1)/(F_1 + M)` is
small. When the surface is saturated the rate is at most :math:`i`, so the
loss is at most :math:`W`.

Representative parameters by soil class are in Rossman and Huber (2016,
Table 4-7, after Rawls et al., 1983). Ponded water at :math:`t = 0` is taken
as rain of rate :math:`h/\Delta t` in the first step; the surface then
saturates when :math:`F` reaches :math:`F_s \approx K_s M \Delta t / h`,
after :math:`K_s M \Delta t^2 / h^2`, instead of at once. Exact ponding
reaches the same :math:`F_s` after half that time, so :math:`F` lags by
about :math:`K_s M \Delta t^2 / (2h^2)`: :math:`3.7 \times 10^{-6}` s in
the test of :ref:`sec:VerifyInfiltration`.

**Curve number** (``curve_number``), the NRCS runoff curve number used as a
loss on the rain (NRCS, 2004; Woodward et al., 2003). With :math:`P` the rain
fallen on the cell since the start, the potential retention :math:`S` and
the initial abstraction :math:`I_a = \lambda S`,

.. math::

   S = 25.4 \left(\frac{1000}{CN} - 10\right)\ \mathrm{mm} , \qquad
   Q = \frac{(P - I_a)^2}{P - I_a + S}\ \ (P > I_a), \quad Q = 0\ \ (P \le I_a),

and the cumulative loss is :math:`L(P) = P - Q`, computed as
:math:`I_a + (P - I_a) S / (P - I_a + S)` above :math:`I_a`. :math:`P` is a
compensated sum, like the rain's, so a run of many small steps keeps it in
single precision. Each step the
cell loses :math:`L(P_{n+1}) - L(P_n)`, which is at most the rain of the step
(:math:`dL/dP \le 1`). Water that runs onto a cell from its neighbours is
not infiltrated: the loss depends on the cell's own rain only.

**Which CN.** The handbook tables of curve numbers assume
:math:`\lambda = 0.2`. FLUME does not convert: the curve number given must
match ``flume.infiltration.cn.ia_ratio``. For :math:`\lambda = 0.05`,
Woodward et al. (2003, eqs. 8 and 9) relate the retentions (in inches) and
the curve numbers by

.. math::

   S_{0.05} = 1.33\, S_{0.2}^{1.15} , \qquad
   CN_{0.05} = \frac{100}{1.879\,(100/CN_{0.2} - 1)^{1.15} + 1} .

For example :math:`CN_{0.2} = 75` gives :math:`CN_{0.05} = 65.31` (their
Table 4).

**Limitations.**

* One storm: the soil does not recover between storms (Rossman and Huber,
  2016, section 4.4.2, not modelled). There a saturated soil takes
  :math:`75/\sqrt{K_s}` hours (:math:`K_s` in in/h) to recover fully: 106 h
  for :math:`K_s = 12.7` mm/h, much longer than a flood of a few hours.
* Green-Ampt leaves out the depth of the ponded water from the head that
  drives infiltration (Rossman and Huber, 2016, eq. 4-27, assume it small
  next to :math:`\psi_s`). Under deep ponding the loss is then low.
* The curve number does not infiltrate run-on water, and it has no rate:
  a cell's loss follows its rain, whatever the soil under ponded water.
* Infiltration is first order in time at both orders of the scheme, like
  the other sources.

.. _sec:Drains:

Storm-drain inlets
------------------

An inlet (``flume.drains.*``) takes water from the surface of the cell that
holds it. The water leaves the model: there is no pipe network and no
surcharge back to the surface. Two capacity laws are available, chosen for
each name with ``flume.drains.<name>.type``: ``weir_orifice`` (the default),
from the depth alone, and ``approach_flow``, from the speed of the water
over the inlet.

**Weir and orifice.** The inlet is a grate in a sag (Brown et al.,
2009, HEC-22, inlets in sag locations). At the depth :math:`d` of its cell it takes the
smaller of the weir and the orifice capacity,

.. math::

   Q_w = C_w P d^{3/2}, \qquad Q_o = C_o A \sqrt{2 g d}, \qquad
   Q = \min(Q_w, Q_o),

with :math:`P` the perimeter of the grate (m), :math:`A` its open area
(m\ :sup:`2`), :math:`g` the gravity (m/s\ :sup:`2`), :math:`C_w` the weir
coefficient (1.66 in SI units) and :math:`C_o` the orifice coefficient
(0.67), the values of HEC-22. The two rates are equal at

.. math::

   d^* = \frac{C_o A \sqrt{2 g}}{C_w P} ;

the weir rate is the smaller below :math:`d^*` and the orifice rate above
it. :math:`d^*` is printed for each inlet at start-up.

**Step.** On a cell of area :math:`|\Omega| = \Delta x \, \Delta y` the
depth follows

.. math::

   \frac{dh}{dt} = -\min\left(a h^{3/2},\, b h^{1/2}\right), \qquad
   a = \frac{C_w P}{|\Omega|}, \qquad b = \frac{C_o A \sqrt{2g}}{|\Omega|} ,

with :math:`d^* = b/a`. The step integrates it exactly from the depth
:math:`h_0` the cell holds after the rain and the infiltration, for the time
:math:`\tau = \Delta t`. Above :math:`d^*` (orifice),

.. math::

   \sqrt{h(t)} = \sqrt{h_0} - \tfrac{1}{2} b t ,

which reaches :math:`d^*` at :math:`t_1 = 2(\sqrt{h_0} - \sqrt{d^*})/b`. At or
below :math:`d^*` (weir), from the depth :math:`h_s` it starts at
(:math:`h_0`, or :math:`d^*` after :math:`t_1`),

.. math::

   h(t) = \frac{h_s}{\left(1 + \tfrac{1}{2} a t \sqrt{h_s}\right)^2} .

The cell loses :math:`h_0 - h(\tau)`, in double precision. The water that
leaves takes its share of the discharge, as for infiltration, and the depth
removed is added to the cell's drained depth (``drained_depth``) by a
compensated sum. The balance column ``drained`` is the drained depth summed
over the cells times the cell area (:doc:`WaterBalance`).

**Why the exact solution.** An explicit step :math:`Q(h_0) \Delta t` can
take more water than the cell holds: a 1 mm film on a 1 m cell with
:math:`P = 100` m would lose 13 mm in the first step of the test, 2.5 s, at
the weir rate. The exact
solution never falls below 0 (the weir branch never empties a cell), so the
loss is at most the water of the cell without a cap, and it does not depend
on the step.

**Approach flow.** With ``type = approach_flow`` the inlet takes

.. math::

   Q = k \, u \, A \, \mathrm{Fr}^{\alpha}, \qquad
   \mathrm{Fr} = \frac{u}{\sqrt{g h}} ,

with :math:`u = \sqrt{(hu)^2 + (hv)^2} / h` the speed of the water of the
cell (m/s), :math:`h` its depth (m), :math:`A` the area of the grate of the
inlet (m\ :sup:`2`; the whole grate, bars and openings, not its open
area), :math:`k` a coefficient and :math:`\alpha` the exponent of the
Froude number. :math:`k` and :math:`\alpha` have no unit.

It is the unified formula of Xia et al. (2022, equation 8), who define a
drainage velocity :math:`U = Q / A` and fit :math:`U / u = k \,
\mathrm{Fr}^{\alpha}` to 78 laboratory runs with one grate (0.75 m by
0.45 m at full size, at the centre of a flume, draining freely): the
defaults :math:`k = 0.302` and :math:`\alpha = -0.816` are their values.
Their runs cover depths of 2.8 to 18.3 cm and Froude numbers of 0.05 to
0.89, with :math:`u` and :math:`h` measured 1 m upstream of the inlet.
Annis et al. (2026, equation 3) use the formula for Brooklyn and calibrated
:math:`k = 0.226` for inlets of 0.4 m\ :sup:`2`.

FLUME takes :math:`\alpha` between -1 and 0, so the capacity does not
decrease with the speed (it goes as :math:`u^{1+\alpha}`) or with the depth
(as :math:`h^{-\alpha/2}`): at :math:`\alpha = -1` it does not depend on the
speed, and at :math:`\alpha = 0` not on the depth.

**Step of the approach-flow law.** The water an inlet takes leaves with its
share of the discharge, so the speed of the cell does not change while it
drains. With :math:`u` fixed over the step,

.. math::

   \frac{dh}{dt} = -c \, h^{p}, \qquad
   c = \frac{k A g^{p} u^{1+\alpha}}{|\Omega|}, \qquad
   p = -\frac{\alpha}{2} ,

and, with :math:`q = 1 - p` (between 1/2 and 1),

.. math::

   h(t)^{q} = h_0^{q} - q \, c \, t ,

which reaches 0 at :math:`t_e = h_0^{q} / (q c)`: unlike the weir, this law
empties a cell in a finite time. With :math:`z = q c \tau / h_0^{q}` the
step takes :math:`h(\tau) = h_0 (1 - z)^{1/q}`, computed as
:math:`h_0 \exp(\ln(1 - z) / q)` with ``log1p``, so a small loss is not
the difference of two close powers. It takes 0 when :math:`z \ge 1`, so the
loss is at most the water of the cell, without a cap. A cell at or below
``flume.dry_depth`` is left as it is. The speed is that of the cell after
the rain, the inflows and the infiltration of the step, and before the
friction.

*The last film.* When :math:`h(\tau)` is above 0 but at or below
``flume.dry_depth``, the step takes all the water of the cell. A film at
that depth has no discharge, so the law could not take it later, and
whether it stayed would depend on where the steps land. The water taken
this way is at most the dry depth per cell, once per emptying, and it is
counted in the drained depth and the balance.

*What it means in still water.* With :math:`\alpha > -1` the capacity is 0
at :math:`u = 0`: an inlet under ponded, still water takes nothing. Only
:math:`u = 0` gives 0: a pond that moves a little, also at round-off
speeds, drains a little, as :math:`u^{1+\alpha}`. A pond that starts at
rest on a level surface stays at rest and keeps its water. At
:math:`\alpha = -1` the speed drops out, :math:`Q = k A \sqrt{g h}`, an
orifice law with :math:`C_o = k / \sqrt{2}`. :math:`k` has no upper limit,
so at :math:`\alpha = -1` a :math:`k` above :math:`\sqrt{2}` is an orifice
coefficient above 1.

**Inlets from a file.** With ``flume.drains.<name>.file`` a name stands for
all the inlets of a list, with one type and one set of parameters. A cell
that holds :math:`n` of them drains at :math:`n` times the rate of one (the
coefficients :math:`a`, :math:`b` and :math:`c` are multiplied by
:math:`n`).

**Limitations.**

* :math:`d` and :math:`u` are the depth and the speed of the cell, not
  those over the grate: an inlet smaller than the cell sees the means of
  the cell.
* A cell holds the inlets of one name. The ``weir_orifice`` capacity is
  that of a grate in a sag (no approach flow, no clogging factor, no curb
  opening); the ``approach_flow`` capacity is 0 in still water unless
  :math:`\alpha = -1`. The two are not combined in one cell.
* The ``approach_flow`` law is empirical; its coefficient is site-specific
  and is meant to be calibrated (Annis et al., 2026). Xia et al. (2022)
  note that it depends on where the inlet stands and on where the speed
  and the depth are taken.
* The law is used outside the range of its fit. The fit covers Froude
  numbers of 0.05 to 0.89; FLUME applies it at every speed. The capacity
  falls slowly with the speed: with the default exponent, as
  :math:`u^{0.184}`, so water at 1 mm/s drains at 28 % of the rate at 1 m/s.
* Exponents below -1 are not allowed. Xia et al. (2022, figure 5) fit
  other grates with exponents of -1.15 to -1.25 (data of other authors). With such an exponent the capacity grows without
  limit as the speed goes to 0, so these fits cannot be given to FLUME as
  they are.
* The inlets of Xia et al. drain freely. A sewer under pressure takes less.
* The ``approach_flow`` capacity does not depend on the width of the cell.
  On cells not much wider than the inlet it can be more than the discharge
  that crosses the cell; the cell then empties and the inlet takes what
  arrives.
* The sewer has no limit: there is no pipe network, so an inlet takes its
  capacity whatever the state of the sewer below it.
* HEC-22 does not give the capacity between the depths where the weir and
  the orifice formulas hold; the smaller of the two is taken there.
* The water does not come back: an inlet always drains while its cell is
  wet.
* The source is first order in time at both orders of the scheme, like the
  other sources.

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
