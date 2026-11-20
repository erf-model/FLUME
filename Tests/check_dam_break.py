#!/usr/bin/env python3
"""Check a FLUME run that starts from a dam: mass, and the exact solution.

    check_dam_break.py <prefix> --time T --direction x|y --position X
        --surface-lo A --surface-hi B [--gravity G] [--mass-tol-eps K]
        [--time-rtol R]
        ( --l1-max E --plateau-rtol P --shock-cells C | --l1-max E --front-max F
        | --first-step-cfl N | --mass-only --min-change D | --initial-only )
        [--velocity-lo U V] [--velocity-hi U V]

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step, whose number is not pinned (it depends on the
precision of the build). Both are read (Tests/check_plotfile.py holds the
reader) and checked:

  * the first plotfile holds the dam: in a cell whose centre along the dam
    direction is below X the depth is max(A - bed, 0) and the discharges are
    that depth times --velocity-lo (default 0 0); in the other cells (centre
    at X or above) they come from B and --velocity-hi. The bed is the one in
    the plotfile. Within 8 epsilons of the larger of the surface, the bed and
    1 m (times 1 + the speed, for the discharges);
  * the time of the last plotfile is T, within R times T;
  * mass: the sum of the depth over all cells is the same in the two plotfiles
    within K machine epsilons (of the precision of the data) times the sum.
    The sums are exact sums of the stored numbers.

Without --mass-only the bed must be flat at 0 with still water on both sides
of the dam, and the last plotfile is compared with the exact solution
of a dam break over a wet bed (a rarefaction, a constant state and a shock). With xi = (x - X) / T measured from the
high water towards the low water, h1 the deeper depth, h0 the shallower,
c = sqrt(G h):

    xi < -c1               h = h1,                    u = 0
    -c1 < xi < um - cm     h = (2 c1 - xi)^2 / (9 G),  u = 2 (c1 + xi) / 3
    um - cm < xi < S       h = hm,                    u = um
    S < xi                 h = h0,                    u = 0

    um = 2 (c1 - cm) = (hm - h0) sqrt(G (hm + h0) / (2 hm h0)),  S = hm um / (hm - h0)

hm solves the middle equation (bisection here). Checked:

  * the state does not vary across the dam direction: depth and discharge
    along the dam differ from the first row by at most 64 epsilons of h1, and
    the discharge along the dam is at most that;
  * L1 error: the mean of |h - h_exact| over the cells, divided by h1 - h0,
    and the mean of |q - q_exact| divided by hm um, are each at most E;
  * plateau: in the middle 60 % of the constant state between the rarefaction
    and the shock, every depth is within P hm of hm and every velocity within
    P um of um;
  * shock: the depth crosses (hm + h0) / 2 within C cells of X + S T.

When the shallower side is dry (h0 = 0) the solution is Ritter's (Delestre et
al., 2013, section 4.1.2):
the rarefaction above reaches down to h = 0 at xi = 2 c1, and beyond it the bed
is dry. Then --front-max F is given in place of --plateau-rtol and
--shock-cells, and the checks are:

  * no depth is negative;
  * the state does not vary across the dam direction, as above;
  * L1 error: the mean of |h - h_exact| divided by h1, and the mean of
    |q - q_exact| divided by h1 c1, are each at most E;
  * beyond the exact front, X + 2 c1 T, no depth is above F h1. A first-order
    scheme carries a trace of water one cell per step, ahead of the front.

With --first-step-cfl N the run is the dry-bed dam break stopped after its
first step, taken at the Courant number N (flume.cfl); --time is not given.
At a face with still water h1 on one side and none on the other the flux uses
the dry-bed wave speeds, -c1 and 2 c1, and is 2 c1 h1 / 3 of mass and
G h1^2 / 3 of momentum. The step is N / (c1 / dx + c1 / dy): the water is
still, and both directions count. With n = N dy / (dx + dy) for a dam along x
(dx and dy exchanged along y), after the step, exactly,

    first dry cell     h = 2 n h1 / 3,        q = n c1 h1 / 3
    last wet cell      h = h1 (1 - 2 n / 3),  q = n c1 h1 / 6

and every other cell is as it started. Checked within 64 epsilons of h1 (depth)
and of c1 h1 (discharge), together with the time (within 8 epsilons, since
the step is computed in the precision of the build) and the mass.

With --initial-only the run takes no step: there is one plotfile, and only the
check of the dam in it is made. --time is not given.

With --mass-only the comparison with an exact solution is left out (the bed may be rough and
the water may move at the start), and the depth must have changed somewhere by
at least D between the two plotfiles: the run did move water.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_plotfile as cp  # noqa: E402


def read(plotfile):
    for f in (plotfile / "Header", plotfile / "Level_0" / "Cell_H"):
        if not f.is_file():
            sys.exit(f"FAIL: {f} does not exist")
    h = cp.read_header(plotfile / "Header")
    if h["vars"] != cp.VARS:
        sys.exit(f"FAIL: {plotfile}: variables {h['vars']}")
    _, fab_on_disk = cp.read_cell_h(plotfile / "Level_0" / "Cell_H")
    cells, real_size = cp.read_cell_data(plotfile / "Level_0", fab_on_disk)
    nx = h["domain_hi"][0] + 1
    ny = h["domain_hi"][1] + 1
    if len(cells) != nx * ny:
        sys.exit(f"FAIL: {plotfile}: {len(cells)} cells in the data files, expected {nx * ny}")
    return h, cells, real_size


def wet_bed_middle(h1, h0, g):
    """hm and um of the constant state between the rarefaction and the shock."""
    c1 = math.sqrt(g * h1)

    def f(hm):
        return 2.0 * (c1 - math.sqrt(g * hm)) - (hm - h0) * math.sqrt(0.5 * g * (hm + h0) / (hm * h0))
    lo, hi = h0, h1  # f(lo) > 0 > f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    hm = 0.5 * (lo + hi)
    return hm, 2.0 * (c1 - math.sqrt(g * hm))


def wet_bed(xi, h1, h0, hm, um, g):
    """Depth and velocity at xi = (x - X) / T; the deep water is at xi < 0."""
    c1 = math.sqrt(g * h1)
    cm = math.sqrt(g * hm)
    shock = hm * um / (hm - h0)
    if xi <= -c1:
        return h1, 0.0
    if xi <= um - cm:
        return (2.0 * c1 - xi) ** 2 / (9.0 * g), 2.0 * (c1 + xi) / 3.0
    if xi <= shock:
        return hm, um
    return h0, 0.0


def ritter(xi, h1, g):
    """Depth and velocity over a dry bed at xi = (x - X) / T; the water is at xi < 0."""
    c1 = math.sqrt(g * h1)
    if xi <= -c1:
        return h1, 0.0
    if xi <= 2.0 * c1:
        return (2.0 * c1 - xi) ** 2 / (9.0 * g), 2.0 * (c1 + xi) / 3.0
    return 0.0, 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("prefix")
    p.add_argument("--time", type=float)
    p.add_argument("--direction", choices=["x", "y"], required=True)
    p.add_argument("--position", type=float, required=True)
    p.add_argument("--surface-lo", type=float, required=True)
    p.add_argument("--surface-hi", type=float, required=True)
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--mass-tol-eps", type=float, default=16.0)
    p.add_argument("--time-rtol", type=float, default=1e-12)
    p.add_argument("--l1-max", type=float)
    p.add_argument("--plateau-rtol", type=float)
    p.add_argument("--shock-cells", type=float)
    p.add_argument("--front-max", type=float)
    p.add_argument("--first-step-cfl", type=float)
    p.add_argument("--mass-only", action="store_true")
    p.add_argument("--initial-only", action="store_true")
    p.add_argument("--velocity-lo", type=float, nargs=2, default=[0.0, 0.0])
    p.add_argument("--velocity-hi", type=float, nargs=2, default=[0.0, 0.0])
    p.add_argument("--min-change", type=float)
    a = p.parse_args()
    if a.initial_only:
        if a.time is not None or a.mass_only or a.first_step_cfl is not None:
            sys.exit("check_dam_break.py: --initial-only goes with no --time, --mass-only or --first-step-cfl")
    elif a.first_step_cfl is not None:
        if a.time is not None or min(a.surface_lo, a.surface_hi) != 0.0 or a.mass_only:
            sys.exit("check_dam_break.py: --first-step-cfl needs a dry side, and no --time or --mass-only")
    elif a.time is None:
        sys.exit("check_dam_break.py: --time is required")
    if a.mass_only:
        if a.min_change is None:
            sys.exit("check_dam_break.py: --mass-only needs --min-change")
    elif a.first_step_cfl is not None or a.initial_only:
        pass
    elif min(a.surface_lo, a.surface_hi) == 0.0:
        if None in (a.l1_max, a.front_max):
            sys.exit("check_dam_break.py: a dry side needs --l1-max and --front-max, or --mass-only")
    elif None in (a.l1_max, a.plateau_rtol, a.shock_cells):
        sys.exit("check_dam_break.py: give --l1-max, --plateau-rtol and --shock-cells, or --mass-only")

    errors = []
    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    n_expected = 1 if a.initial_only else 2
    if len(found) != n_expected or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and {n_expected - 1} later plotfile(s), found {found}")
    step = int(found[-1][len(a.prefix):])
    head0, first, real_size = read(Path(found[0]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23

    # The dam in the first plotfile.
    d_dam = 0 if a.direction == "x" else 1
    worst_init = 0.0
    for c, v in first.items():
        x = head0["prob_lo"][d_dam] + (c[d_dam] + 0.5) * head0["dx"][d_dam]
        low = x < a.position
        surface = a.surface_lo if low else a.surface_hi
        u = a.velocity_lo if low else a.velocity_hi
        h = max(surface - v[3], 0.0)
        scale = max(abs(surface), abs(v[3]), 1.0)
        worst_init = max(worst_init, abs(v[0] - h) / scale,
                         abs(v[1] - h * u[0]) / (scale * (1.0 + abs(u[0]))),
                         abs(v[2] - h * u[1]) / (scale * (1.0 + abs(u[1]))))
    if not worst_init <= 8.0 * eps:
        errors.append(f"the first plotfile is off the dam initial state by {worst_init / eps:.3g} epsilons,"
                      " at most 8")
    if a.initial_only:
        if errors:
            print(f"FAIL: {a.prefix}")
            for e in errors:
                print("  " + e)
            sys.exit(1)
        print(f"PASS: {found[0]}: the dam initial state within {worst_init / eps:.3g} epsilons")
        return
    head, last, _ = read(Path(found[1]))
    if a.first_step_cfl is not None:
        c1 = math.sqrt(a.gravity * max(a.surface_lo, a.surface_hi))
        a.time = a.first_step_cfl / (c1 / head["dx"][0] + c1 / head["dx"][1])
        if step != 1:
            errors.append(f"the last plotfile is of step {step}, expected 1")
        # The step is computed in the precision of the build.
        a.time_rtol = max(a.time_rtol, 8.0 * eps)
    if head["step"] != step:
        errors.append(f"step: got {head['step']}, expected {step}")
    if not abs(head["time"] - a.time) <= a.time_rtol * a.time:
        errors.append(f"time: got {head['time']!r}, expected {a.time!r}")
    for (i, j), v in last.items():
        if not all(math.isfinite(x) for x in v):
            errors.append(f"cell ({i}, {j}) holds a value that is not finite")
            break

    # Mass: exact sums of the stored depths.
    m0 = math.fsum(v[0] for v in first.values())
    m1 = math.fsum(v[0] for v in last.values())
    mass_tol = a.mass_tol_eps * eps * m0
    if not abs(m1 - m0) <= mass_tol:
        errors.append(f"mass: the depths sum to {m1!r}, at the start to {m0!r}: difference"
                      f" {m1 - m0:.3g}, tolerance {mass_tol:.3g}")
    report = f"mass change {(m1 - m0) / (eps * m0):.3g} epsilons of the mass"

    if a.first_step_cfl is not None:
        d = 0 if a.direction == "x" else 1
        dx = head["dx"][d]
        h1 = max(a.surface_lo, a.surface_hi)
        c1 = math.sqrt(a.gravity * h1)
        cfl = a.time * c1 / dx   # Courant number along the dam direction
        sign = 1 if a.surface_lo > a.surface_hi else -1   # the way the water runs
        m_dam = (a.position - head["prob_lo"][d]) / dx    # index of the first cell above the dam
        if abs(m_dam - round(m_dam)) > 1e-9:
            sys.exit("check_dam_break.py: --first-step-cfl needs the dam on a face")
        m_dam = round(m_dam)
        m_wet, m_dry = (m_dam - 1, m_dam) if sign > 0 else (m_dam, m_dam - 1)
        expected = {m_wet: (h1 * (1.0 - 2.0 * cfl / 3.0), sign * cfl * c1 * h1 / 6.0),
                    m_dry: (2.0 * cfl * h1 / 3.0, sign * cfl * c1 * h1 / 3.0)}
        worst_h = 0.0
        worst_q = 0.0
        for c, v in last.items():
            h_ex, q_ex = expected.get(c[d], (first[c][0], 0.0))
            worst_h = max(worst_h, abs(v[0] - h_ex) / h1)
            worst_q = max(worst_q, abs(v[1 + d] - q_ex) / (c1 * h1), abs(v[2 - d]) / (c1 * h1))
        if not (worst_h <= 64.0 * eps and worst_q <= 64.0 * eps):
            errors.append(f"after one step the depth is off by {worst_h / eps:.3g} epsilons of h1 and the"
                          f" discharge by {worst_q / eps:.3g} epsilons of c1 h1; at most 64."
                          f" Expected (h, q) in the cells at the dam: {expected}")
        report += (f"; first step: depth within {worst_h / eps:.3g} epsilons, discharge within"
                   f" {worst_q / eps:.3g} epsilons of the closed form")
    elif a.mass_only:
        change = max(abs(last[c][0] - first[c][0]) for c in first)
        if not change >= a.min_change:
            errors.append(f"the depth changed by at most {change:.3g} m, expected at least {a.min_change}")
        report += f"; largest depth change {change:.3g} m"
    else:
        d = 0 if a.direction == "x" else 1
        n = head["domain_hi"][d] + 1
        n_across = head["domain_hi"][1 - d] + 1
        dx = head["dx"][d]
        g = a.gravity
        deep_is_lo = a.surface_lo > a.surface_hi
        h1 = max(a.surface_lo, a.surface_hi)
        h0 = min(a.surface_lo, a.surface_hi)
        dry = h0 == 0.0
        c1 = math.sqrt(g * h1)
        # Over a dry bed the "constant state" is empty: hm = 0 at um = 2 c1.
        hm, um = (0.0, 2.0 * c1) if dry else wet_bed_middle(h1, h0, g)
        shock = um if dry else hm * um / (hm - h0)
        sign = 1.0 if deep_is_lo else -1.0

        def cell(m, t):
            return last[(m, t)] if d == 0 else last[(t, m)]

        # Uniform across the dam direction; no discharge along the dam.
        across_tol = 64.0 * eps * h1
        worst_across = 0.0
        for m in range(n):
            ref = cell(m, 0)
            for t in range(n_across):
                v = cell(m, t)
                worst_across = max(worst_across, abs(v[0] - ref[0]), abs(v[1 + d] - ref[1 + d]),
                                   abs(v[2 - d]))
        if not worst_across <= across_tol:
            errors.append(f"the state varies along the dam by {worst_across:.3g}, tolerance {across_tol:.3g}")

        l1_h = 0.0
        l1_q = 0.0
        worst_plateau = 0.0
        cm = math.sqrt(g * hm)
        width = shock - (um - cm)
        p_lo = (um - cm) + 0.2 * width
        p_hi = shock - 0.2 * width
        n_plateau = 0
        crossing = None
        prev = None
        lowest = math.inf
        ahead = 0.0
        for m in range(n):
            x = head["prob_lo"][d] + (m + 0.5) * dx
            xi = sign * (x - a.position) / a.time
            h_ex, u_ex = ritter(xi, h1, g) if dry else wet_bed(xi, h1, h0, hm, um, g)
            h = cell(m, 0)[0]
            q = sign * cell(m, 0)[1 + d]  # discharge towards the low water
            l1_h += abs(h - h_ex)
            l1_q += abs(q - h_ex * u_ex)
            if dry:
                lowest = min(lowest, h)
                if xi > shock:
                    ahead = max(ahead, h)
                continue
            if p_lo <= xi <= p_hi:
                n_plateau += 1
                worst_plateau = max(worst_plateau, abs(h - hm) / hm, abs(q / h - um) / um)
            # Where the depth crosses the middle of the shock, at xi > 0.
            level = 0.5 * (hm + h0)
            if prev is not None and xi > 0.0 and prev[1] > 0.0 and (prev[2] - level) * (h - level) <= 0.0 \
                    and prev[2] != h:
                crossing = prev[0] + (x - prev[0]) * (level - prev[2]) / (h - prev[2])
            prev = (x, xi, h)
        l1_h /= n * (h1 - h0)
        l1_q /= n * (h1 * c1 if dry else hm * um)
        if not (l1_h <= a.l1_max and l1_q <= a.l1_max):
            errors.append(f"L1 error: depth {l1_h:.4g}, discharge {l1_q:.4g}, at most {a.l1_max}")
        if dry:
            if not lowest >= 0.0:
                errors.append(f"the lowest depth is {lowest!r}: negative")
            if not ahead <= a.front_max * h1:
                errors.append(f"beyond the exact front the depth reaches {ahead:.3g} m, at most"
                              f" {a.front_max * h1:.3g} m")
            report += (f"; L1 error {l1_h:.4g} (depth), {l1_q:.4g} (discharge); lowest depth {lowest:.3g} m;"
                       f" at most {ahead:.3g} m beyond the exact front at"
                       f" {a.position + sign * shock * a.time:.4f} m; variation along the dam {worst_across:.3g}")
        elif n_plateau < 10:
            errors.append(f"only {n_plateau} cells in the constant state: the test needs at least 10")
        if not dry and not worst_plateau <= a.plateau_rtol:
            errors.append(f"constant state: depth or velocity off by {worst_plateau:.4g} of hm = {hm:.6g} m,"
                          f" um = {um:.6g} m/s; at most {a.plateau_rtol}")
        x_shock = a.position + sign * shock * a.time
        if dry:
            pass
        elif crossing is None:
            errors.append("the depth does not cross the middle of the shock")
        elif not abs(crossing - x_shock) <= a.shock_cells * dx:
            errors.append(f"shock at {crossing:.6g} m, expected {x_shock:.6g} m within {a.shock_cells} cells of {dx} m")
        if not dry:
            report += (f"; L1 error {l1_h:.4g} (depth), {l1_q:.4g} (discharge); constant state within"
                       f" {worst_plateau:.3g} over {n_plateau} cells; shock at"
                       f" {crossing if crossing is None else round(crossing, 4)} m (exact {x_shock:.4f} m);"
                       f" variation along the dam {worst_across:.3g}")

    if errors:
        print(f"FAIL: {a.prefix}")
        for e in errors:
            print("  " + e)
        print("  (" + report + ")")
        sys.exit(1)
    print(f"PASS: {found[1]}: {report}")


if __name__ == "__main__":
    main()
