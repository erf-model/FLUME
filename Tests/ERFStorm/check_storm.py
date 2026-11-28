#!/usr/bin/env python3
"""Check the storm case: the coupled run, its replay and the hyetograph run.

Usage (in the directory of the coupled run, as the check of Tests/run_case.py):

    check_storm.py <index> <dt> <step> --replay <command...>
                   --hyetograph <command...>

The coupled run (inputs_storm) took the rain of ERF's Kessler microphysics
for <step> ERF steps of <dt> s, one FLUME step each, and wrote it to the
index <index> (flume.rain.replay_file). The checks:

rain: make_hyetograph.py writes hyetograph.txt, the catchment-mean rate of
each grid of the index. ERF's accumulated surface rain (precip_total_accum,
kg/m2 = mm, of each plt2d plotfile, every 60 steps) summed over the
catchment's 8 x 10 atmospheric columns of 2 km must equal the integral of
the hyetograph up to the plotfile's time times the catchment's area, to 1e-12
of the rain on the catchment at <step>. The hyetograph thus holds the rain
ERF made over the catchment, at the steps it made it.

replay: the stand-alone run with flume.rain.type = gridded, in replay/, must
give the coupled flood at <step> / 2 and <step> to round-off: the checks of
Tests/ERFReplay/check_replay.py (64 eps (h + R) for the depth, 64 eps of the
largest discharge, the rain of the water balance).

hyetograph: the same stand-alone run with flume.rain.type = hyetograph and
hyetograph.txt, in hyeto/, must exit with 0 and close its water balance
(Tests/check_balance.py, with the factor BALANCE_FACTOR = 256 in place of
64: the hyetograph wets every cell of the domain at once and, from about 18
to 36 min, all its rain infiltrates. The rain is summed over the 8192 cells of a
rank at each step and the infiltration from the cumulative depth of each
cell. The hyetograph makes these values equal in every cell, so the
roundings of each sum add up in one direction, up to about 1000 eps of the
volume, which the bound of check_balance.py does not hold. The largest ratio
measured is about 99). The hyetograph rains on the whole domain, so the
rain of its water balance, at each line (each step), must be the
hyetograph's depth up to then times the domain's area, to 1e-12 of the
total: the rates hold from one entry to the next (flume.rain.interpolation
= step, the default), and the catchment then gets the coupled rain, spread
evenly over it, at the steps ERF made it.

hydrographs: the discharge at the outlet is the change of out_outlet of
water_balance.txt over each step. The replay's peak must be at least 3 times
the hyetograph run's and come at least 1 h earlier, and both peaks must come
at least 1 h before the end of the runs. These bounds are loose on purpose:
the cell is convection, and its rain moves with the compiler and the
processor; the README gives the numbers of one run.
"""
import bisect
import math
import subprocess
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
from check_plotfile import read_cell_data, read_cell_h, read_header  # noqa: E402
from prepare_catchment import CATCHMENT  # noqa: E402

MM_PER_H = 1.0 / 3.6e6  # m/s per mm/h
PEAK_RATIO = 3.0
PEAK_LEAD = 3600.0
END_MARGIN = 3600.0
BALANCE_FACTOR = 256.0


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def read_hyetograph(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.split("#")[0].split()
        if s:
            rows.append((float(s[0]), float(s[1])))
    return rows


class Hyetograph:
    """The depth of rain of a hyetograph up to a time, m: each rate holds from
    its entry to the next, and the last from its entry on, as
    flume.rain.interpolation = step in FLUME_TimeTable; 0 before the first."""

    def __init__(self, rows):
        self.times = [t for t, _ in rows]
        self.rates = [r * MM_PER_H for _, r in rows]
        self.cum = [0.0]
        for k in range(len(rows) - 1):
            self.cum.append(self.cum[-1] + self.rates[k] * (self.times[k + 1] - self.times[k]))

    def depth(self, t):
        k = bisect.bisect_right(self.times, t) - 1
        return 0.0 if k < 0 else self.cum[k] + self.rates[k] * (t - self.times[k])


def erf_catchment_rain(plotfile, box):
    """ERF's accumulated rain over the columns inside box, m3, and the time."""
    p = Path(plotfile)
    h = read_header(p / "Header")
    _, fab = read_cell_h(p / "Level_0" / "Cell_H")
    cells, _ = read_cell_data(p / "Level_0", fab)
    if h["vars"] != ["precip_total_accum"]:
        fail(f"{plotfile} holds {h['vars']}, not precip_total_accum alone")
    lo, hi = h["prob_lo"], h["prob_hi"]
    n = [h["domain_hi"][d] - h["domain_lo"][d] + 1 for d in (0, 1)]
    dx = [(hi[d] - lo[d]) / n[d] for d in (0, 1)]
    parts = []
    count = 0
    for (i, j), v in cells.items():
        x = lo[0] + (i + 0.5) * dx[0]
        y = lo[1] + (j + 0.5) * dx[1]
        if box[0] <= x < box[1] and box[2] <= y < box[3]:
            parts.append(v[0] * 1e-3 * dx[0] * dx[1])
            count += 1
    return h["time"], math.fsum(parts), count, dx[0] * dx[1]


def check_rain(index, nstep):
    r = subprocess.run([sys.executable, str(HERE / "make_hyetograph.py"), str(index),
                        "hyetograph.txt"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    print(r.stdout, end="")
    if r.returncode != 0:
        fail("make_hyetograph.py exited with " + str(r.returncode))
    rows = Hyetograph(read_hyetograph("hyetograph.txt"))
    area = (CATCHMENT[1] - CATCHMENT[0]) * (CATCHMENT[3] - CATCHMENT[2])
    plotfiles = sorted(Path(".").glob("plt2d?????"))
    if not plotfiles or plotfiles[-1].name != f"plt2d{nstep:05d}":
        fail(f"the coupled run wrote no plt2d{nstep:05d}")
    results = []
    for p in plotfiles:
        t, erf, count, col = erf_catchment_rain(p, CATCHMENT)
        if count * col != area:
            fail(f"{p}: the catchment covers {count} columns of {col:g} m2, not {area:g} m2")
        results.append((p.name, t, erf, rows.depth(t) * area))
    total = results[-1][2]
    if total <= 0.0:
        fail("no rain fell on the catchment")
    bound = 1e-12 * total
    for name, t, erf, hy in results:
        if abs(erf - hy) > bound:
            fail(f"{name} ({t:g} s): ERF's rain on the catchment is {erf!r} m3, the "
                 f"hyetograph's {hy!r} m3: they differ by more than {bound:.3g} m3")
    worst = max(abs(erf - hy) for _, _, erf, hy in results)
    print(f"PASS: ERF's rain on the catchment ({area / 1e6:g} km2) is that of hyetograph.txt "
          f"at each of the {len(results)} plt2d plotfiles: {total!r} m3 "
          f"({total / area * 1e3:.4g} mm) at {results[-1][1]:g} s, largest difference "
          f"{worst:.2e} m3, bound {bound:.2e} m3")
    return rows, total, area


def outlet_hydrograph(where):
    _, names, rows = read_log(Path(where) / "water_balance.txt")
    t = [r[names.index("time")] for r in rows]
    v = [r[names.index("out_outlet")] for r in rows]
    q = [((t[k] + t[k + 1]) / 2, (v[k + 1] - v[k]) / (t[k + 1] - t[k]))
         for k in range(len(t) - 1)]
    return q, t[-1]


def check_hyetograph(command, rows, area):
    """The stand-alone run of hyetograph.txt in hyeto/."""
    here = Path(".").resolve()
    check_replay.run(command + [f"flume.rain.file={(here / 'hyetograph.txt')}"], here / "hyeto")
    ok, msg, _ = check_balance.check(here / "hyeto" / "water_balance.txt", BALANCE_FACTOR)
    if not ok:
        fail(f"hyeto/: {msg}")
    print(f"PASS: hyeto/: {msg}")
    h, _ = check_replay.flood("hyeto", 0)
    domain = (h["prob_hi"][0] - h["prob_lo"][0]) * (h["prob_hi"][1] - h["prob_lo"][1])
    _, names, log = read_log(here / "hyeto" / "water_balance.txt")
    i_t, i_r = names.index("time"), names.index("rain")
    worst = 0.0
    scale = rows.depth(log[-1][i_t]) * domain
    for line in log:
        t = line[i_t]
        want = rows.depth(t) * domain
        worst = max(worst, abs(line[i_r] - want))
        if abs(line[i_r] - want) > 1e-12 * scale:
            fail(f"hyeto/: the rain of the water balance at {t:g} s is {line[i_r]!r} m3, not "
                 f"the hyetograph's depth times the domain's {domain:g} m2, {want!r} m3")
    total = log[-1][i_r]
    print(f"PASS: hyeto/: the rain of the water balance is the hyetograph's depth times the "
          f"domain's area at each of its {len(log)} lines (largest difference {worst:.2e} m3 "
          f"of {total!r} m3); {scale * area / domain!r} m3 of it on the catchment")


def check_hydrographs():
    qr, end_r = outlet_hydrograph("replay")
    qh, end_h = outlet_hydrograph("hyeto")
    tr, pr = max(qr, key=lambda p: p[1])
    th, ph = max(qh, key=lambda p: p[1])
    print(f"replay:     peak {pr:.4g} m3/s at {tr:g} s (run to {end_r:g} s)")
    print(f"hyetograph: peak {ph:.4g} m3/s at {th:g} s (run to {end_h:g} s)")
    if ph <= 0.0 or pr < PEAK_RATIO * ph:
        fail(f"the replay's peak, {pr:.4g} m3/s, is not {PEAK_RATIO:g} times the "
             f"hyetograph run's, {ph:.4g} m3/s")
    if tr > th - PEAK_LEAD:
        fail(f"the replay peaks at {tr:g} s, not {PEAK_LEAD:g} s before the hyetograph run "
             f"({th:g} s)")
    for name, t, end in (("replay", tr, end_r), ("hyetograph", th, end_h)):
        if t > end - END_MARGIN:
            fail(f"the {name} run peaks at {t:g} s, within {END_MARGIN:g} s of its end {end:g} s")
    print(f"PASS: the replay's peak is {pr / ph:.3g} times the hyetograph run's and "
          f"{(th - tr) / 3600:.3g} h earlier; both peak at least {END_MARGIN / 3600:g} h "
          f"before the end")


def main():
    args = sys.argv[1:]
    usage = "usage: check_storm.py <index> <dt> <step> --replay <command...> --hyetograph <command...>"
    if "--replay" not in args or "--hyetograph" not in args:
        fail(usage)
    i = args.index("--replay")
    j = args.index("--hyetograph")
    if not (i == 3 and j > i + 1 and len(args) > j + 1):
        fail(usage)
    index = Path(args[0])
    dt = float(args[1])
    nstep = int(args[2])
    if not index.is_file():
        fail(f"the coupled run wrote no {index}")
    if not math.isfinite(dt) or dt <= 0 or nstep < 2 or nstep % 2:
        fail("<dt> must be > 0 and <step> even")
    rows, _, area = check_rain(index, nstep)
    r_max = check_replay.check_index(index, dt, nstep)
    check_replay.check_replay(index, nstep, args[i + 1:j], r_max)
    check_hyetograph(args[j + 1:], rows, area)
    check_hydrographs()


if __name__ == "__main__":
    main()
