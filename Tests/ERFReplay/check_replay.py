#!/usr/bin/env python3
"""Check the replay of the rain FLUME took from ERF (flume.rain.replay_file).

Usage (in the directory of the coupled run, as the check of Tests/run_case.py):

    check_replay.py <index> <dt> <step> --replay <command...>
                    --restart <checkpoint> <command...>

The coupled run here took the rain of ERF's microphysics for <step> ERF
steps of <dt> s, one FLUME step each, and wrote it to the index <index>
(flume.rain.replay_file). The checks:

index: every line is a time and a grid in the index's directory; the times
are multiples of <dt> and increase; the first line is the grid of no rain
(<stem>_zero.asc, all 0) at time 0 and the last the same grid at <step> * <dt>
(the end of the coupled rain); a grid of rain is followed by a line one step
later, a grid of no rain only by a grid of rain (dry steps add no line); no
grid of rain is all 0. The rain of the grids, rate * 1/3.6e6 * the time to the
next line * the cell area, summed, must be the rain of water_balance.txt at
<step> to 1e-12 of it plus what each cell may still hold in the carry of the
Rain source (FLUME_Rain.cpp: the part of a step's rain that the addition to
the depth rounds away waits for the next step), one eps of the largest depth
per cell times the cell area: the writer, without the reader.

replay: the command (the stand-alone run of the same flood, with
flume.rain.type = gridded) is run in replay/ with flume.rain.file = the
index; it must exit with 0. The depth of every cell of flood<n> for n = <step>
/ 2 and <step> must be that of the coupled run to 64 eps (h + R), with h the
coupled depth of the cell and R the largest rain a cell received in the
coupled run (the rate of each step is written with 17 digits, read back and
multiplied by the step: within 2 eps of the step's rain per step; a cell
whose depth the rain changes by less than one unit in its last place may
round to the neighbouring double; the flow carries these differences on, and
64 is a margin measured for this case, not a bound for every case); both
discharges to 64 eps of the largest discharge of either component; and the
rain of water_balance.txt to 64 eps of it. The checks assume one FLUME step
per ERF step and steps of <dt> that sum exactly in binary, as here.

recession: the same command in recession/ with flume.stop_time = 2 * <step> *
<dt> must reach that time, add no rain after the end of the coupled rain (the
rain column at step 2 * <step> is that of step <step>, to the bit), and move
water (some depth changes by 1e-6 m from flood<step>).

restart: the coupled run continues from <checkpoint> with the command. The
grids of the lines from the checkpoint's time on are removed first, and the
index is left whole: the restart must cut it back and write those grids
again. The index and every grid must then be those of the first run, to the
byte. The restart must print "restarted from <checkpoint>/Level_0/FLUME".

aborts: the restart in abort_index/ without the index, in abort_stale/ with
the index of the first run cut back to half the checkpoint's time (as a
restart from an older checkpoint that was stopped leaves it), and in
abort_held/ from a copy of the checkpoint without rain_replay_held, must stop
with the messages of flume.rain.replay_file (twice) and flume.restart.
"""
import filecmp
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TESTS))
from check_balance import read_log  # noqa: E402
from check_plotfile import read_cell_data, read_cell_h, read_header  # noqa: E402

EPS = 2.220446049250313e-16
MM_PER_H = 1.0 / 3.6e6


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def read_grid(path):
    """Header {key: value} and rows (top first) of an ESRI ASCII grid."""
    head = {}
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0][0].isalpha():
            head[parts[0].lower()] = float(parts[1])
        else:
            rows.append([float(v) for v in parts])
    return head, rows


def read_index(index):
    lines = []
    for n, line in enumerate(index.read_text().splitlines(), 1):
        parts = line.split()
        if len(parts) != 2:
            fail(f"{index}: line {n} is not a time and a grid: {line!r}")
        lines.append((float(parts[0]), parts[1]))
    return lines


def flood(where, step):
    p = Path(where) / f"flood{step:05d}"
    _, fab = read_cell_h(p / "Level_0" / "Cell_H")
    cells, _ = read_cell_data(p / "Level_0", fab)
    return read_header(p / "Header"), cells


def rain_at(where, step):
    _, names, rows = read_log(Path(where) / "water_balance.txt")
    for r in rows:
        if int(r[0]) == step:
            return r[names.index("rain")], r[names.index("time")]
    fail(f"{where}/water_balance.txt has no line at step {step}")


def run(command, where, expect_abort=None):
    where.mkdir(exist_ok=True)
    print(f"in {where.name}/: {' '.join(command)}", flush=True)
    r = subprocess.run(command, cwd=where, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    (where / "run.log").write_text(r.stdout)
    if expect_abort is None:
        if r.returncode != 0:
            print(r.stdout[-3000:])
            fail(f"the run in {where.name}/ exited with {r.returncode}")
    else:
        if r.returncode == 0:
            fail(f"the run in {where.name}/ did not stop")
        if not re.search(expect_abort, r.stdout):
            print(r.stdout[-3000:])
            fail(f"the run in {where.name}/ did not print {expect_abort!r}")
        print(f"PASS: {where.name}/ stopped with {expect_abort!r}")
    return r.stdout


def check_index(index, dt, nstep):
    lines = read_index(index)
    stem = index.stem
    zero = f"{stem}_zero.asc"
    t_end = nstep * dt
    if lines[0] != (0.0, zero) or lines[-1] != (t_end, zero):
        fail(f"{index}: the first line must be '0 {zero}' and the last '{t_end:g} {zero}'; "
             f"they are {lines[0]} and {lines[-1]}")
    for (t, g), (t2, g2) in zip(lines, lines[1:]):
        if not t2 > t:
            fail(f"{index}: the time {t2} does not follow {t}")
        if g != zero and t2 != t + dt:
            fail(f"{index}: the grid {g} at {t} s is not followed by a line at {t + dt} s")
        if g == zero and g2 == zero:
            fail(f"{index}: two grids of no rain follow each other at {t} and {t2} s")
    for t, _ in lines:
        if t / dt != round(t / dt):
            fail(f"{index}: the time {t} is not a multiple of {dt}")
    head, rows = read_grid(index.parent / zero)
    if any(v != 0.0 for r in rows for v in r):
        fail(f"{zero} is not all 0")
    area = head["cellsize"] ** 2
    parts = []
    largest = {}
    for (t, g), (t2, _) in zip(lines, lines[1:]):
        if g == zero:
            continue
        _, rows = read_grid(index.parent / g)
        if max(max(r) for r in rows) <= 0.0:
            fail(f"{g} holds no rain")
        for j, r in enumerate(rows):
            for i, v in enumerate(r):
                d = v * MM_PER_H * (t2 - t)
                parts.append(d * area)
                largest[(i, j)] = largest.get((i, j), 0.0) + d
    total = math.fsum(parts)
    rain, _ = rain_at(".", nstep)
    _, cells = flood(".", nstep)
    carry = len(cells) * EPS * max(v[0] for v in cells.values()) * area
    bound = 1e-12 * rain + carry
    ngrid = sum(1 for _, g in lines if g != zero)
    if abs(total - rain) > bound:
        fail(f"the rain of the grids of {index}, {total!r} m3, is not the rain of "
             f"water_balance.txt, {rain!r} m3, to {bound!r} m3")
    print(f"PASS: {index}: {len(lines)} lines, {ngrid} grids of rain, closed at {t_end:g} s; "
          f"their rain {total!r} m3 is that of water_balance.txt, {rain!r} m3 "
          f"(difference {abs(total - rain):.2e} m3, bound {bound:.2e} m3)")
    return max(largest.values())


def check_replay(index, nstep, command, r_max):
    here = Path(".").resolve()
    run(command + [f"flume.rain.file={index.resolve()}"], here / "replay")
    for step in (nstep // 2, nstep):
        ha, a = flood(".", step)
        hb, b = flood("replay", step)
        if ha["time"] != hb["time"]:
            fail(f"flood{step:05d}: time {hb['time']} in replay/, {ha['time']} here")
        worst = 0.0
        for k, va in a.items():
            vb = b[k]
            bound = 64 * EPS * (abs(va[0]) + r_max)
            worst = max(worst, abs(va[0] - vb[0]) / bound)
        qmax = max(max(abs(v[1]), abs(v[2])) for v in a.values())
        for comp in (1, 2):
            dq = max(abs(a[k][comp] - b[k][comp]) for k in a)
            if dq > 64 * EPS * qmax:
                fail(f"flood{step:05d}: {ha['vars'][comp]} differs by {dq!r}, above 64 eps of "
                     f"the largest discharge {qmax!r}")
        if worst > 1.0:
            fail(f"flood{step:05d}: a depth differs by {worst:.3g} times the bound 64 eps (h + R)")
        print(f"PASS: flood{step:05d} of the replay is the coupled flood: the largest "
              f"difference of depth is {worst:.2e} of the bound 64 eps (h + R), R = {r_max!r} m")
    ra, _ = rain_at(".", nstep)
    rb, _ = rain_at("replay", nstep)
    if abs(ra - rb) > 64 * EPS * ra:
        fail(f"the rain of the replay, {rb!r} m3, differs from the coupled {ra!r} m3")
    print(f"PASS: the rain of the replay, {rb!r} m3, is the coupled {ra!r} m3 "
          f"(difference {abs(ra - rb) / ra:.2e} of it)")


def check_recession(index, dt, nstep, command):
    here = Path(".").resolve()
    t_stop = 2 * nstep * dt
    run(command + [f"flume.rain.file={index.resolve()}", f"flume.stop_time={t_stop!r}"],
        here / "recession")
    r1, _ = rain_at("recession", nstep)
    r2, t2 = rain_at("recession", 2 * nstep)
    if t2 != t_stop:
        fail(f"the recession ended at {t2} s, not {t_stop} s")
    if r2 != r1:
        fail(f"the recession added rain after the coupled rain ended: {r1!r} m3 at step "
             f"{nstep}, {r2!r} m3 at step {2 * nstep}")
    _, a = flood("recession", nstep)
    _, b = flood("recession", 2 * nstep)
    moved = max(abs(a[k][0] - b[k][0]) for k in a)
    if moved < 1e-6:
        fail(f"no water moved in the recession (largest change of depth {moved!r} m)")
    print(f"PASS: the recession runs to {t_stop:g} s without rain after {nstep * dt:g} s "
          f"(rain {r2!r} m3) and the largest depth changes by {moved:.3g} m")


def check_restart(index, chk, command):
    here = Path(".").resolve()
    rain_dir = index.parent
    ref = here / "rain_first_run"
    if ref.exists():
        shutil.rmtree(ref)
    shutil.copytree(rain_dir, ref)
    t_chk = None
    for line in (Path(chk) / "Level_0" / "FLUME" / "Header").read_text().splitlines():
        if line.startswith("time "):
            t_chk = float(line.split()[1])
    if t_chk is None:
        fail(f"{chk}/Level_0/FLUME/Header holds no time")
    zero = f"{index.stem}_zero.asc"
    removed = [g for t, g in read_index(index) if t >= t_chk and g != zero]
    if not removed:
        fail(f"{index} has no grid at or after {t_chk} s: the restart would rewrite nothing")
    for g in removed:
        (rain_dir / g).unlink()
    print(f"removed {len(removed)} grids from {t_chk} s on; the index is left whole")

    # The aborts first, in their own directories, from the same checkpoint.
    a = here / "abort_index"
    a.mkdir(exist_ok=True)
    if not (a / chk).exists():
        shutil.copytree(chk, a / chk)
    run(command, a, expect_abort=r"flume\.rain\.replay_file: \S+: cannot be opened")
    c = here / "abort_stale"
    c.mkdir(exist_ok=True)
    if not (c / chk).exists():
        shutil.copytree(chk, c / chk)
    if (c / rain_dir.name).exists():
        shutil.rmtree(c / rain_dir.name)
    shutil.copytree(ref, c / rain_dir.name)
    t_cut = t_chk / 2
    kept = [f"{t!r} {g}\n" for t, g in read_index(index) if t < t_cut]
    (c / rain_dir.name / index.name).write_text("".join(kept) + f"{t_cut!r} {zero}\n")
    run(command, c, expect_abort=r"flume\.rain\.replay_file: \S+: does not hold what the run "
                                 r"wrote up to the checkpoint")
    b = here / "abort_held"
    b.mkdir(exist_ok=True)
    if (b / chk).exists():
        shutil.rmtree(b / chk)
    shutil.copytree(chk, b / chk)
    if (b / rain_dir.name).exists():
        shutil.rmtree(b / rain_dir.name)
    shutil.copytree(ref, b / rain_dir.name)
    h = b / chk / "Level_0" / "FLUME" / "Header"
    h.write_text("".join(s for s in h.read_text().splitlines(True)
                         if not s.startswith("multifab/rain_replay_held ")))
    run(command, b, expect_abort=r"flume\.restart: \S+: it was written without "
                                 r"flume\.rain\.replay_file")

    print(f"restart: {' '.join(command)}", flush=True)
    r = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    Path("restart.log").write_text(r.stdout)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        fail(f"the restart exited with {r.returncode}")
    if not re.search(rf"restarted from {re.escape(chk)}/Level_0/FLUME at step", r.stdout):
        fail(f"the restart did not print 'restarted from {chk}/Level_0/FLUME'")
    names = sorted(p.name for p in ref.iterdir())
    if sorted(p.name for p in rain_dir.iterdir()) != names:
        fail(f"after the restart {rain_dir} holds other files than after the first run")
    match, mismatch, errors = filecmp.cmpfiles(ref, rain_dir, names, shallow=False)
    if mismatch or errors:
        fail(f"after the restart these files differ from those of the first run: "
             f"{' '.join(mismatch + errors)}")
    print(f"PASS: after the restart from {chk} the index and the {len(match) - 1} grids are "
          f"those of the first run, to the byte")


def main():
    args = sys.argv[1:]
    if "--replay" not in args or "--restart" not in args:
        fail("usage: check_replay.py <index> <dt> <step> --replay <command...> "
             "--restart <checkpoint> <command...>")
    i = args.index("--replay")
    j = args.index("--restart")
    if not (i == 3 and j > i + 1 and len(args) > j + 2):
        fail("usage: check_replay.py <index> <dt> <step> --replay <command...> "
             "--restart <checkpoint> <command...>")
    index = Path(args[0])
    dt = float(args[1])
    nstep = int(args[2])
    replay = args[i + 1:j]
    chk = args[j + 1]
    restart = args[j + 2:]
    if not index.is_file():
        fail(f"the coupled run wrote no {index}")
    if not math.isfinite(dt) or dt <= 0 or nstep < 2 or nstep % 2:
        fail("<dt> must be > 0 and <step> even")
    r_max = check_index(index, dt, nstep)
    check_replay(index, nstep, replay, r_max)
    check_recession(index, dt, nstep, replay)
    check_restart(index, chk, restart)


if __name__ == "__main__":
    main()
