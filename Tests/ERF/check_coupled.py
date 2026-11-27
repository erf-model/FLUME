#!/usr/bin/env python3
"""Check the flood plotfiles and logs of a run against those of another run.

Usage (in the directory of the run, as the check of Tests/run_case.py):

    check_coupled.py <prefix> --dt <s> [--substeps]
    check_coupled.py <prefix> --dt <s> [--substeps] --reference <dir> --log <file>...
                     [--column-mean <ERF plotfile> <flood plotfile> <ratio>]
                     [--restart <checkpoint> <command...>]

ERF names a flood plotfile with its own step, and FLUME writes its own step
and time in the header. Each plotfile <prefix>NNNNN here must carry FLUME
step N and time N * dt exactly, with N the step in its name, and there must
be at least two. FLUME takes no step longer than flume.max_dt = dt, so a time
of N * dt after N steps means that every step was dt: in a stand-alone run
that no step was cut short by the CFL or rain limit, and in a run inside ERF
that each ERF step was one FLUME step. This is all that is checked without
--reference (the stand-alone run).

With --substeps, dt is ERF's step and FLUME takes its own, shorter steps
inside it: each plotfile must carry time N * dt exactly (FLUME landed on the
end of every ERF step) and a FLUME step above N (after step 0). It goes with
--reference and --restart too, for a run with substeps checked against
another such run.

With --reference, every plotfile here must also be in <dir>, and the other
way round, with the same value in every cell, the same time and the same step
(Tests/compare_plotfiles.py --cells); and each log file (water_balance.txt,
gauge_<name>.txt) must be byte for byte that of <dir>.

With --column-mean, each field flood_<name> of the ERF plotfile must hold,
in the lowest layer of each column, the mean of <name> (depth, discharge_x,
discharge_y, bed_elevation, surface_elevation) over the ratio x ratio flood
cells of the column in the flood plotfile, within 4 machine epsilons of the
largest magnitude of that field, and 0 in the layers above. At least two
fields must be there.

With --restart, the run here has written the ERF checkpoint <checkpoint>
(ending in its step n). The flood plotfiles after step n are removed, the
logs are kept (they run past step n: the restart must cut them back), and
the command, which continues from that checkpoint, is run here. It must exit
with 0 and print FLUME's "restarted from <checkpoint>/Level_0/FLUME at step
k", where k is the step in <checkpoint>/Level_0/FLUME/Header: n, or with
--substeps the number of FLUME's steps so far, more than n. Then the
plotfiles and logs are checked against <dir> again.
"""
import filecmp
import struct
import re
import shutil
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TESTS))
from check_plotfile import FAB_HEAD, read_cell_data, read_cell_h, read_header  # noqa: E402
from compare_plotfiles import plotfiles  # noqa: E402


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def check_steps(prefix, dt, substeps=False):
    names = plotfiles(".", prefix)
    if len(names) < 2:
        fail(f"{len(names)} plotfiles {prefix}NNNNN here, at least 2 needed")
    for name in names:
        n = int(name[len(prefix):])
        h = read_header(Path(name) / "Header")
        if substeps:
            if h["time"] != n * dt or (h["step"] <= n if n > 0 else h["step"] != 0):
                fail(f"{name}: FLUME step {h['step']} at {h['time']!r} s, expected more than {n} "
                     f"steps at {n * dt!r} s (ERF's step {dt} s in FLUME steps)")
            continue
        if h["step"] != n or h["time"] != n * dt:
            fail(f"{name}: FLUME step {h['step']} at {h['time']!r} s, expected step {n} at "
                 f"{n * dt!r} s (every step {dt} s)")
    if substeps:
        print(f"PASS: {len(names)} plotfiles, each at ERF time N * {dt} s after more than N FLUME steps")
    else:
        print(f"PASS: {len(names)} plotfiles, each at FLUME step N and time N * {dt} s")
    return names


def check_reference(prefix, ref, logs):
    r = subprocess.run([sys.executable, str(TESTS / "compare_plotfiles.py"), prefix, str(ref),
                        "--cells"])
    if r.returncode != 0:
        fail(f"the flood plotfiles differ from those of {ref}")
    for log in logs:
        if not Path(log).is_file() or not (ref / log).is_file():
            fail(f"{log}: missing here or in {ref}")
        if not filecmp.cmp(log, ref / log, shallow=False):
            fail(f"{log} differs from {ref / log}")
    print(f"PASS: the flood plotfiles and {', '.join(logs)} are those of {ref}")


def read_layers(level_dir):
    """Every cell {(i, j, k): [one value per component]} of a level, all layers."""
    _, fab_on_disk = read_cell_h(level_dir / "Cell_H")
    cells = {}
    for name, offset in fab_on_disk:
        with open(level_dir / name, "rb") as f:
            f.seek(offset)
            m = FAB_HEAD.fullmatch(f.readline())
            if not m:
                fail(f"{level_dir / name}: cannot read the FAB header")
            size = int(m.group(3))
            endian = "<" if int(m.group(4).split()[0]) == size else ">"
            lo = [int(m.group(k)) for k in (5, 6, 7)]
            hi = [int(m.group(k)) for k in (8, 9, 10)]
            ncomp = int(m.group(14))
            ni, nj, nk = (hi[d] - lo[d] + 1 for d in range(3))
            n = ni * nj * nk
            v = struct.unpack(f"{endian}{n * ncomp}{'d' if size == 8 else 'f'}",
                              f.read(n * ncomp * size))
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        idx = (k * nj + j) * ni + i
                        cells[(lo[0] + i, lo[1] + j, lo[2] + k)] = [v[c * n + idx] for c in range(ncomp)]
    return cells, size


def check_column_mean(erf_plot, flood_plot, ratio):
    he, hf = read_header(Path(erf_plot) / "Header"), read_header(Path(flood_plot) / "Header")
    erf_cells, size = read_layers(Path(erf_plot) / "Level_0")
    _, fab = read_cell_h(Path(flood_plot) / "Level_0" / "Cell_H")
    flood_cells, _ = read_cell_data(Path(flood_plot) / "Level_0", fab)
    eps = 2.0 ** -52 if size == 8 else 2.0 ** -23
    k0 = he["domain_lo"][2]
    fields = [v for v in he["vars"] if v.startswith("flood_")]
    if len(fields) < 2:
        fail(f"{erf_plot} holds the flood fields {fields}: at least 2 needed")
    for field in fields:
        c_e, c_f = he["vars"].index(field), hf["vars"].index(field[len("flood_"):])
        scale = max(abs(v[c_f]) for v in flood_cells.values())
        worst = 0.0
        for (i, j, k), v in erf_cells.items():
            if k != k0:
                if v[c_e] != 0.0:
                    fail(f"{erf_plot} {field} at {(i, j, k)}: {v[c_e]!r} above the lowest layer, not 0")
                continue
            mean = sum(flood_cells[(ratio * i + a, ratio * j + b)][c_f]
                       for a in range(ratio) for b in range(ratio)) / ratio ** 2
            worst = max(worst, abs(v[c_e] - mean))
        if worst > 4.0 * eps * scale:
            fail(f"{erf_plot} {field}: differs from the mean of {flood_plot} by {worst!r}, "
                 f"bound {4.0 * eps * scale!r}")
        print(f"PASS: {erf_plot} {field} is the mean over {ratio} x {ratio} cells of {flood_plot} "
              f"(largest difference {worst:.3g})")


def main():
    args = sys.argv[1:]
    restart = []
    if "--restart" in args:
        i = args.index("--restart")
        restart = args[i + 1:]
        args = args[:i]
    column = None
    if "--column-mean" in args:
        i = args.index("--column-mean")
        if i + 3 >= len(args):
            fail("--column-mean takes an ERF plotfile, a flood plotfile and the ratio")
        column = (args[i + 1], args[i + 2], int(args[i + 3]))
        del args[i:i + 4]
    logs = []
    if "--log" in args:
        i = args.index("--log")
        logs = args[i + 1:]
        args = args[:i]
    substeps = "--substeps" in args
    args = [a for a in args if a != "--substeps"]
    ref = None
    if "--reference" in args:
        i = args.index("--reference")
        if i + 1 >= len(args):
            fail("--reference takes a directory")
        ref = Path(args[i + 1])
        del args[i:i + 2]
    if len(args) != 3 or args[1] != "--dt" or bool(ref) != bool(logs) or (restart and not ref):
        fail("usage: check_coupled.py <prefix> --dt <s> [--substeps] [--reference <dir> --log <file>... "
             "[--restart <checkpoint> <command...>]]")
    prefix, dt = args[0], float(args[2])

    check_steps(prefix, dt, substeps)
    if ref is None:
        return
    if not ref.is_dir() or ref.resolve() == Path(".").resolve():
        fail(f"no reference directory {ref}, or it is this directory")
    check_reference(prefix, ref, logs)
    if column:
        check_column_mean(*column)
    if not restart:
        return

    chk, command = restart[0], restart[1:]
    m = re.search(r"(\d+)$", chk)
    if not m or not (Path(chk) / "Header").is_file() or not command:
        fail(f"--restart needs a checkpoint written here (no {chk}/Header) and a command")
    n = int(m.group(1))
    removed = [p for p in plotfiles(".", prefix) if int(p[len(prefix):]) > n]
    if not removed:
        fail(f"no flood plotfile after the checkpoint's step {n}: the restart would show nothing")
    for p in removed:
        shutil.rmtree(p)
    print(f"removed {' '.join(removed)}; restart: {' '.join(command)}", flush=True)
    run = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         encoding="utf-8", errors="backslashreplace")
    print(run.stdout, flush=True)
    if run.returncode != 0:
        fail(f"the restart exited with {run.returncode}")
    # FLUME's step in its checkpoint: n, or with substeps the number of
    # FLUME's steps in ERF's first n, which is more than n.
    header = Path(chk) / "Level_0" / "FLUME" / "Header"
    words = dict(line.split(None, 1) for line in header.read_text().splitlines() if " " in line) \
        if header.is_file() else {}
    if "step" not in words:
        fail(f"no 'step' line in {header}")
    k = int(words["step"])
    if (k <= n) if substeps else (k != n):
        fail(f"{header}: FLUME's step {k}, {'more than ' if substeps else ''}{n} expected")
    m = re.search(rf"restarted from {re.escape(chk)}/Level_0/FLUME at step (\d+)\b", run.stdout)
    if not m or int(m.group(1)) != k:
        fail(f"FLUME did not report continuing from {chk}/Level_0/FLUME at step {k}")
    check_steps(prefix, dt, substeps)
    check_reference(prefix, ref, logs)


if __name__ == "__main__":
    main()
