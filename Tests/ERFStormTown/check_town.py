#!/usr/bin/env python3
"""Check the town case: the coupled run, its stand-alone run and the run
without buildings.

Usage (in the directory of the coupled run, as the check of Tests/run_case.py):

    check_town.py <town> <step> --standalone <command...>
                  --no-buildings <command...>

<town> is the directory prepare_town.py wrote (window.txt, rain.txt,
buildings.asc). The coupled run (inputs_town) took <step> ERF steps of
ERF_DT = 0.75 s, each of FLUME_PER_ERF = 3 FLUME steps of 0.25 s, and wrote
the ERF plotfiles plt<n> and the flood plotfiles flood<n> at n = <step> / 2
and <step>. The checks:

rain: the rain of the coupled water_balance.txt, at each of its lines, must
be the integral of the storm's grids over the square and over the window up
to the line's time (each rate held to the next entry, as rain.txt is read
with flume.rain.interpolation = step, FLUME's default), to 1e-12 of the total. FLUME maps each
grid to its cells by the mean over each cell (FLUME_RainField.H), so the rain
on the square is that of the 16 storm cells of 250 m it covers.

balance: the coupled water balance closes (Tests/check_balance.py, with the
factor BALANCE_FACTOR = 1024 in place of 64) and has FLUME_PER_ERF * <step>
steps, to the window's end. The factor: the square holds 4 x 4 storm cells
of 250 m, each 100 x 100 flood cells that get the same rain at every step,
and the bed is a plane, so the water is a sheet whose depth is the same in
many cells. FLUME sums the volumes and the rain of each step over the 40000
cells of a rank; the roundings of a sum of equal values add up in one
direction rather than at random, which the bound of check_balance.py does not
hold. The largest ratios measured are 247 with the buildings (in the first
minute; below 64 after 4 min) and 518 without them (the sheet stays even).
1024 is about twice the largest measured ratio, not a bound: it depends on
the cells of each sum (the ranks and boxes). The residual stays within
1.6e-9 m3 of 2.9e4 m3 and does not drift.

footprints: the buildings of ERF and of FLUME are the same. On the ground
level of plt<step>, the ERF cells with terrain_IB_mask = 1 (solid) must be
the ERF columns whose 4 x 4 flood cells are all building cells of flood<step>
(bed_elevation above the plane bed); a flood cell of a building must lie in
a column with terrain_IB_mask above 0 (the column or its rim, which
erf.if_snap_partial_cells makes half solid). In each solid column the height
of the solid, 10 m times the sum of terrain_IB_mask up the column, must be
within half an ERF cell (5 m) of the height of FLUME's building there.

immersed forcing: the air inside ERF's buildings is at rest: at plt<step>,
the largest speed in a solid cell (terrain_IB_mask = 1) whose six neighbours
are solid too (the ground counts as solid) must be below IB_FRACTION of the
largest speed of the domain. A solid cell next to the air is left out: its
plotted speed is the mean of its faces, and the faces it shares with the air
carry the air's speed (up to 1.97 m/s in one run). ERF sets the momentum of
the faces inside a building to zero after its projection, so this checks
that ERF's buildings are on and that freeze is active. It does not test the
drag coefficient (erf.if_Cd_momentum = 0 passes it; the README says so).

standalone: the command (FLUME alone on the same flood grid, the same
buildings and rain) is run in standalone/; its flood plotfiles at FLUME
steps FLUME_PER_ERF * <step> / 2 and FLUME_PER_ERF * <step> must be the
coupled flood<n> to the bit, and its water_balance.txt the coupled one, line
for line. FLUME uses its own bed inside ERF; the buildings ERF sees do not
reach it. The water balances agree byte for byte because ERF's flood grid and
the stand-alone run put the same box on each rank (four boxes of 200 x 200
cells on four ranks); another rank count or mapping changes the order of
the sums.

no buildings: the command (the same without flume.buildings.file) is run in
no_buildings/: its water balance must close (factor BALANCE_FACTOR), its
rain must be the coupled rain to 64 eps, and the buildings must change the flood: the largest depth
outside the buildings must be above CHANGE[0] = 2 times that without them,
and the outflow above CHANGE[1] = 1.05 times (4.0 and 1.11 in one run). The
README gives the numbers of one run.
"""
import math
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TESTS = HERE.parent
sys.path.insert(0, str(TESTS))
sys.path.insert(0, str(TESTS / "ERFReplay"))
sys.path.insert(0, str(HERE))
import check_balance  # noqa: E402
import check_replay  # noqa: E402
from check_balance import read_log  # noqa: E402
from check_plotfile import FAB_HEAD, read_cell_h, read_header  # noqa: E402
from prepare_town import ERF_DX, FLOOD_DX, Grid, read_index, square_mean  # noqa: E402

MM_PER_H = 1.0 / 3.6e6  # m/s per mm/h
ERF_DT = 0.75
FLUME_PER_ERF = 3
SLOPE_Y = 0.01  # flume.bed.slope of flume_town, rise over run in y
IB_FRACTION = 0.05
ROOF = 1.0e-3  # a flood cell is a building cell above this height, m (the lowest is 10 m)
CHANGE = (2.0, 1.05)  # the buildings at least double the street depth and add 5 % to the outflow
BALANCE_FACTOR = 1024.0
EPS = 2.0 ** -52


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def read_plotfile(path):
    """Header and {(i, j, k): [one value per component]} of level 0."""
    p = Path(path)
    h = read_header(p / "Header")
    _, fab_on_disk = read_cell_h(p / "Level_0" / "Cell_H")
    cells = {}
    for name, offset in fab_on_disk:
        with open(p / "Level_0" / name, "rb") as f:
            f.seek(offset)
            head = f.readline()
            m = FAB_HEAD.fullmatch(head)
            if not m:
                fail(f"{p / 'Level_0' / name}: cannot read the FAB header {head!r}")
            size = int(m.group(3))
            order = [int(x) for x in m.group(4).split()]
            if size not in (4, 8) or order not in (list(range(size, 0, -1)),
                                                   list(range(1, size + 1))):
                fail(f"{p / 'Level_0' / name}: unknown real format in {head!r}")
            lo = [int(m.group(k)) for k in (5, 6, 7)]
            hi = [int(m.group(k)) for k in (8, 9, 10)]
            ncomp = int(m.group(14))
            ni, nj, nk = (hi[d] - lo[d] + 1 for d in range(3))
            n = ni * nj * nk
            raw = f.read(n * ncomp * size)
            if len(raw) != n * ncomp * size:
                fail(f"{p / 'Level_0' / name}: the data of box {lo} {hi} are cut short")
            endian = "<" if order[0] == size else ">"
            v = struct.unpack(f"{endian}{n * ncomp}{'d' if size == 8 else 'f'}", raw)
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        cells[(lo[0] + i, lo[1] + j, lo[2] + k)] = [
                            v[c * n + (k * nj + j) * ni + i] for c in range(ncomp)]
    return h, cells


def read_window(town):
    s = (Path(town) / "window.txt").read_text().split()
    x0, y0, size, t0, window, depth = (float(v) for v in s)
    return x0, y0, size, t0, window, depth


class SquareRain:
    """The rain on the square up to a time, m3, from the index rain.txt."""

    def __init__(self, index, x0, y0, size):
        rows = read_index(index)
        self.times = [t for t, _ in rows]
        self.rates = [square_mean(Grid(p), x0, y0) * MM_PER_H for _, p in rows]
        self.area = size * size

    def volume(self, t):
        parts = []
        for k, tk in enumerate(self.times):
            if tk >= t:
                break
            end = min(self.times[k + 1], t) if k + 1 < len(self.times) else t
            parts.append(self.rates[k] * (end - tk))
        return math.fsum(parts) * self.area


def check_rain(town, x0, y0, size, window):
    rain = SquareRain(Path(town) / "rain.txt", x0, y0, size)
    _, names, rows = read_log("water_balance.txt")
    i_t, i_r = names.index("time"), names.index("rain")
    total = rain.volume(window)
    if total <= 0.0:
        fail("no rain in the window")
    worst = 0.0
    for r in rows:
        want = rain.volume(r[i_t])
        d = abs(r[i_r] - want)
        worst = max(worst, d)
        if d > 1e-12 * total:
            fail(f"the rain of the water balance at {r[i_t]:g} s is {r[i_r]!r} m3, not the "
                 f"storm's {want!r} m3 over the square")
    print(f"PASS: the rain of the coupled water balance is that of the storm's grids over the "
          f"square at each of its {len(rows)} lines: {rows[-1][i_r]!r} m3 "
          f"({rows[-1][i_r] / rain.area * 1e3:.4g} mm) at {rows[-1][i_t]:g} s, largest "
          f"difference {worst:.2e} m3, bound {1e-12 * total:.2e} m3")
    return rows[-1][i_r]


def check_coupled_balance(nstep, window):
    ok, msg, _ = check_balance.check(Path("water_balance.txt"), BALANCE_FACTOR)
    if not ok:
        fail(f"the coupled water balance: {msg}")
    _, names, rows = read_log("water_balance.txt")
    last = rows[-1]
    if int(last[0]) != FLUME_PER_ERF * nstep or last[names.index("time")] != window:
        fail(f"the coupled water balance ends at step {last[0]:.0f}, {last[names.index('time')]!r}"
             f" s, not step {FLUME_PER_ERF * nstep} at {window!r} s")
    print(f"PASS: the coupled water balance: {msg}; {FLUME_PER_ERF} FLUME steps per ERF step")


def flood_buildings(flood_cells, hdr):
    """FLUME's building height of each flood cell, m: bed_elevation less the plane bed."""
    ib = hdr["vars"].index("bed_elevation")
    out = {}
    for (i, j, _), v in flood_cells.items():
        y = (j + 0.5) * FLOOD_DX
        out[(i, j)] = v[ib] - SLOPE_Y * y
    return out


def check_footprints(nstep):
    he, erf = read_plotfile(f"plt{nstep:05d}")
    hf, fl = read_plotfile(f"flood{nstep:05d}")
    if "terrain_IB_mask" not in he["vars"]:
        fail(f"plt{nstep:05d} has no terrain_IB_mask: ERF has no immersed buildings "
             "(erf.buildings_type = ImmersedForcing, erf.plot_vars_1)")
    im = he["vars"].index("terrain_IB_mask")
    height = flood_buildings(fl, hf)
    r = int(round(ERF_DX / FLOOD_DX))
    ni = he["domain_hi"][0] - he["domain_lo"][0] + 1
    nj = he["domain_hi"][1] - he["domain_lo"][1] + 1
    nk = he["domain_hi"][2] - he["domain_lo"][2] + 1
    if (hf["domain_hi"][0] + 1, hf["domain_hi"][1] + 1) != (r * ni, r * nj):
        fail(f"the flood grid is not ERF's columns refined by {r}")
    solid = 0
    worst_h = 0.0
    for J in range(nj):
        for I in range(ni):
            m0 = erf[(I, J, 0)][im]
            sub = [height[(r * I + a, r * J + b)] for b in range(r) for a in range(r)]
            built = [h > ROOF for h in sub]
            if (m0 == 1.0) != all(built):
                fail(f"ERF column ({I}, {J}): terrain_IB_mask {m0!r} at the ground, but "
                     f"{sum(built)} of its {r * r} flood cells are building cells")
            if any(built) and not m0 > 0.0:
                fail(f"ERF column ({I}, {J}) is open (terrain_IB_mask 0) but holds "
                     f"{sum(built)} building cells of FLUME")
            if m0 == 1.0:
                solid += 1
                h_erf = ERF_DX * math.fsum(erf[(I, J, k)][im] for k in range(nk))
                h_fl = math.fsum(sub) / len(sub)
                worst_h = max(worst_h, abs(h_erf - h_fl))
                if abs(h_erf - h_fl) > 0.5 * ERF_DX:
                    fail(f"ERF column ({I}, {J}): the solid is {h_erf:g} m high, FLUME's "
                         f"building {h_fl:g} m")
    if solid == 0:
        fail("ERF has no solid column")
    print(f"PASS: ERF's {solid} solid columns are FLUME's building columns; their heights "
          f"agree within {worst_h:g} m (bound {0.5 * ERF_DX:g} m)")
    return he, erf


def check_forcing(he, erf):
    im = he["vars"].index("terrain_IB_mask")
    iu = [he["vars"].index(v) for v in ("x_velocity", "y_velocity", "z_velocity")]

    def speed(v):
        return math.sqrt(sum(v[c] ** 2 for c in iu))

    # A cell outside the domain counts as air (the x seam is not wrapped): this
    # can only leave a cell out of the check. No building reaches the domain's sides.
    def solid(i, j, k):
        if k < 0:
            return True
        v = erf.get((i, j, k))
        return v is not None and v[im] == 1.0

    inner = [speed(v) for (i, j, k), v in erf.items()
             if v[im] == 1.0 and all(solid(*n) for n in ((i + 1, j, k), (i - 1, j, k),
                                                         (i, j + 1, k), (i, j - 1, k),
                                                         (i, j, k + 1), (i, j, k - 1)))]
    edge = max(speed(v) for v in erf.values() if v[im] == 1.0)
    top = max(speed(v) for v in erf.values())
    if not inner:
        fail("ERF has no solid cell inside a building")
    inside = max(inner)
    if not inside < IB_FRACTION * top:
        fail(f"the largest speed inside ERF's buildings is {inside:.3g} m/s, not below "
             f"{IB_FRACTION:g} of the domain's largest, {top:.3g} m/s")
    print(f"PASS: the immersed forcing holds the air inside the buildings at rest: largest "
          f"speed {inside:.3g} m/s in {len(inner)} solid cells away from the air "
          f"({edge:.3g} m/s in a solid cell next to it), {top:.3g} m/s in the domain")


def flood2d(where, step):
    h, cells = check_replay.flood(where, step)
    return h, cells


def check_standalone(nstep, command):
    here = Path(".").resolve()
    check_replay.run(command, here / "standalone")
    for n in (nstep // 2, nstep):
        ha, a = flood2d(".", n)
        hb, b = flood2d("standalone", FLUME_PER_ERF * n)
        if ha["time"] != hb["time"]:
            fail(f"flood{n:05d} is at {ha['time']!r} s, the stand-alone plotfile at "
                 f"{hb['time']!r} s")
        if not a or a.keys() != b.keys():
            fail(f"flood{n:05d} and the stand-alone plotfile do not hold the same cells")
        diff = [k for k in a if a[k] != b[k]]
        if diff:
            k = diff[0]
            fail(f"flood{n:05d}: {len(diff)} cells differ from the stand-alone run, "
                 f"first {k}: {a[k]} against {b[k]}")
        print(f"PASS: flood{n:05d} is the stand-alone flood at step {FLUME_PER_ERF * n}, "
              f"{ha['time']:g} s, to the bit ({len(a)} cells)")
    la = Path("water_balance.txt").read_text()
    lb = (here / "standalone" / "water_balance.txt").read_text()
    if la != lb:
        fail("the coupled water_balance.txt is not that of the stand-alone run")
    print("PASS: the coupled water_balance.txt is that of the stand-alone run, byte for byte")


def outflow(where):
    _, names, rows = read_log(Path(where) / "water_balance.txt")
    i_t, i_o = names.index("time"), names.index("out_south")
    q = [((rows[k][i_t] + rows[k + 1][i_t]) / 2,
          (rows[k + 1][i_o] - rows[k][i_o]) / (rows[k + 1][i_t] - rows[k][i_t]))
         for k in range(len(rows) - 1)]
    return q, rows[-1][i_o], rows[-1][names.index("rain")]


def street_depth(where, step, buildings):
    _, cells = check_replay.flood(where, step)
    return max(v[0] for k, v in cells.items() if not buildings.get(k, 0.0) > ROOF)


def check_no_buildings(nstep, command, rain):
    here = Path(".").resolve()
    check_replay.run(command, here / "no_buildings")
    ok, msg, _ = check_balance.check(here / "no_buildings" / "water_balance.txt",
                                     BALANCE_FACTOR)
    if not ok:
        fail(f"no_buildings/: {msg}")
    print(f"PASS: no_buildings/: {msg}")
    qb, out_b, rain_b = outflow(".")
    qn, out_n, rain_n = outflow("no_buildings")
    if abs(rain_n - rain) > 64 * EPS * rain:
        fail(f"no_buildings/ gets {rain_n!r} m3 of rain, the coupled run {rain!r} m3")
    hf, fl = read_plotfile(f"flood{nstep:05d}")
    bld = flood_buildings(fl, hf)
    db = street_depth(".", nstep, bld)
    dn = street_depth("no_buildings", FLUME_PER_ERF * nstep, bld)
    tb, pb = max(qb, key=lambda p: p[1])
    tn, pn = max(qn, key=lambda p: p[1])
    print(f"with buildings:    out {out_b:.6g} m3 of {rain:.6g} m3, peak {pb:.4g} m3/s at "
          f"{tb:g} s, largest depth off the buildings {db:.4g} m")
    print(f"without buildings: out {out_n:.6g} m3 of {rain_n:.6g} m3, peak {pn:.4g} m3/s at "
          f"{tn:g} s, largest depth {dn:.4g} m")
    if not (db > CHANGE[0] * dn and out_b > CHANGE[1] * out_n):
        fail(f"the buildings do not change the flood: the largest street depth must be above "
             f"{CHANGE[0]:g} times, and the outflow above {CHANGE[1]:g} times, those without them")
    print(f"PASS: the buildings change the flood: street depth {db / dn:.3g} times (bound "
          f"{CHANGE[0]:g}), outflow {out_b / out_n:.4g} times (bound {CHANGE[1]:g})")


def main():
    args = sys.argv[1:]
    usage = ("usage: check_town.py <town> <step> --standalone <command...> "
             "--no-buildings <command...>")
    if "--standalone" not in args or "--no-buildings" not in args:
        fail(usage)
    i = args.index("--standalone")
    j = args.index("--no-buildings")
    if not (i == 2 and j > i + 1 and len(args) > j + 1):
        fail(usage)
    town = Path(args[0])
    nstep = int(args[1])
    if nstep < 2 or nstep % 2:
        fail("<step> must be even")
    x0, y0, size, _, window, _ = read_window(town)
    if nstep * ERF_DT != window:
        fail(f"{nstep} ERF steps of {ERF_DT:g} s are not the window of {window:g} s")
    rain = check_rain(town, x0, y0, size, window)
    check_coupled_balance(nstep, window)
    he, erf = check_footprints(nstep)
    check_forcing(he, erf)
    check_standalone(nstep, args[i + 1:j])
    check_no_buildings(nstep, args[j + 1:], rain)


if __name__ == "__main__":
    main()
