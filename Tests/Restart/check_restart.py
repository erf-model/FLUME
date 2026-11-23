#!/usr/bin/env python3
"""Continue a run from its checkpoint and compare it with an uninterrupted run.

Usage (in the directory of the first part of the run, as the check of
Tests/run_case.py):

    check_restart.py --reference <dir> --gauges <n> [--balance-rtol <r>]
                     [--rewound] [--cut-line] --restart <command...>

The first part of the run has already run in this directory and written the
checkpoint that the command's flume.restart=<name> names. The script then

1. with --rewound, checks that the first part logged past the checkpoint
   (so the restart has lines of water_balance.txt and of the gauge files to
   remove); with --cut-line, also adds to each log the start of a line of a
   later step, without its end of line, as a run stopped while writing it
   leaves (the restart must remove it);
2. runs the restart command here; it must exit with 0, print
   "restarted from <name> at step <n>", take step n + 1 first, and write no
   plotfile or checkpoint of step n or before;
3. checks that every plotfile of the reference directory is here, with the
   same step and time and every value of every cell (all variables) the same
   to the last bit, and that the last one comes from the restarted part;
4. checks that the reference has <n> gauge files and that each is byte for
   byte that of the reference;
5. checks water_balance.txt: byte for byte that of the reference with
   --balance-rtol 0 (the default; same rank count and boxes), or else the
   same steps and times, and each value within r times the largest magnitude
   of its column (the sums over the ranks are taken in another order).

    check_restart.py --expect <name>...

checks only that each of the named checkpoints (directories with a Header)
exists here: the check of the reference runs.

    check_restart.py --compare <dir>

checks only that the plotfiles here are those of <dir>, the same set, every
value of every cell to the last bit (step 3 above): the same run on other
ranks or boxes.
"""
import filecmp
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_plotfile as cp  # noqa: E402


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def split(argv):
    out = {"--reference": [], "--gauges": [], "--balance-rtol": [], "--restart": [], "--expect": [],
           "--compare": []}
    flags = {"--rewound": False, "--cut-line": False}
    key = None
    for a in argv:
        if key == "--restart":
            out[key].append(a)  # everything after --restart is the command
            continue
        if a in out:
            key = a
            continue
        if a in flags:
            flags[a] = True
            key = None
            continue
        if key is None:
            fail(f"check_restart.py: unexpected argument {a!r}")
        out[key].append(a)
    return out, flags


def plotfiles(where):
    return sorted(d.name for d in Path(where).iterdir() if d.is_dir() and re.fullmatch(r"plt\d{5,}", d.name))


def read_plotfile(p):
    h = cp.read_header(p / "Header")
    _, fab_on_disk = cp.read_cell_h(p / "Level_0" / "Cell_H")
    cells, _ = cp.read_cell_data(p / "Level_0", fab_on_disk)
    return h, cells


def data_lines(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.split())
    return rows


def log_steps(path):
    return [int(r[0]) for r in data_lines(path)]


def compare_plotfiles(ref, same_set):
    """Every plotfile of `ref` must be here (and no other with `same_set`),
    with the same variables, step and time, and every value of every cell
    the same to the last bit. Returns the names."""
    mine, theirs = plotfiles("."), plotfiles(ref)
    missing = sorted(set(theirs) - set(mine))
    if missing or len(theirs) < 2:
        fail(f"plotfiles of {ref} missing here: {missing} (reference has {theirs})")
    if same_set and mine != theirs:
        fail(f"plotfiles here {mine}, in {ref} {theirs}")
    n = 0
    for name in theirs:
        h_mine, c_mine = read_plotfile(Path(name))
        h_ref, c_ref = read_plotfile(ref / name)
        for k in ("vars", "time", "step"):
            if h_mine[k] != h_ref[k]:
                fail(f"{name}: {k} {h_mine[k]!r} here, {h_ref[k]!r} in {ref}")
        if c_mine.keys() != c_ref.keys():
            fail(f"{name} holds other cells here than in {ref}")
        for c, v in c_mine.items():
            if len(v) != len(h_ref["vars"]) or list(v) != list(c_ref[c]):
                fail(f"{name} cell {c}: {list(v)} here, {list(c_ref[c])} in {ref}")
            n += 1
    print(f"PASS: {len(theirs)} plotfiles, {n} cells, every value identical to those in {ref}")
    return theirs


def main():
    args, flags = split(sys.argv[1:])
    if args["--expect"]:
        for name in args["--expect"]:
            if not (Path(name) / "Header").is_file():
                fail(f"no checkpoint {name} (no {name}/Header)")
        print(f"PASS: checkpoints {' '.join(args['--expect'])} written")
        return
    if args["--compare"]:
        if len(args["--compare"]) != 1:
            fail("check_restart.py: --compare takes one directory")
        ref = Path(args["--compare"][0])
        if not ref.is_dir() or ref.resolve() == Path(".").resolve():
            fail(f"no reference directory {ref}, or it is this directory")
        compare_plotfiles(ref, same_set=True)
        return
    if len(args["--reference"]) != 1 or len(args["--gauges"]) != 1 or not args["--restart"]:
        fail("check_restart.py: --reference <dir>, --gauges <n> and --restart <command...>"
             " are required")
    ref = Path(args["--reference"][0])
    rtol = float(args["--balance-rtol"][0]) if args["--balance-rtol"] else 0.0
    if not ref.is_dir():
        fail(f"no reference directory {ref}")
    if ref.resolve() == Path(".").resolve():
        fail(f"the reference directory {ref} is this directory")

    names = [a.split("=", 1)[1] for a in args["--restart"] if a.startswith("flume.restart=")]
    if len(names) != 1:
        fail("the restart command needs exactly one flume.restart=<checkpoint>")
    chk = names[0]
    m = re.search(r"(\d+)$", chk)
    if not m or not (Path(chk) / "Header").is_file():
        fail(f"the first part wrote no checkpoint {chk}")
    chk_step = int(m.group(1))
    gauges = sorted(p.name for p in ref.glob("gauge_*.txt"))
    if len(gauges) != int(args["--gauges"][0]):
        fail(f"the reference has the gauge files {gauges}, not {args['--gauges'][0]}")
    first_last = max(log_steps("water_balance.txt"))
    if flags["--rewound"]:
        logs = ["water_balance.txt"] + gauges
        for f in logs:
            if max(log_steps(f)) <= chk_step:
                fail(f"the first part logged {f} only to step {max(log_steps(f))}: nothing to rewind"
                     f" past the checkpoint of step {chk_step}")
        print(f"first part: logs to step {first_last}, checkpoint at step {chk_step}")
    if flags["--cut-line"]:
        if not flags["--rewound"]:
            fail("--cut-line needs --rewound")
        for f in ["water_balance.txt"] + gauges:
            with open(f, "a") as fh:
                fh.write(str(chk_step)[:-1])  # e.g. "40" of a line "400 ..." cut short
        print(f"added a line cut short to {', '.join(['water_balance.txt'] + gauges)}")

    print("restart:", " ".join(args["--restart"]), flush=True)
    run = subprocess.run(args["--restart"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         encoding="utf-8", errors="backslashreplace")
    print(run.stdout, flush=True)
    if run.returncode != 0:
        fail(f"the restart exited with {run.returncode}")
    if not re.search(rf"restarted from {re.escape(chk)} at step {chk_step}\b", run.stdout):
        fail(f"the restart did not report continuing from {chk} at step {chk_step}")
    steps = [int(s) for s in re.findall(r"^Step (\d+):", run.stdout, re.M)]
    if not steps or steps[0] != chk_step + 1:
        fail(f"the restart's first step is {steps[:1]}, not {chk_step + 1}")
    written = [(n, int(s)) for n, s in re.findall(r"^Wrote (\S*?)(\d+) at time", run.stdout, re.M)]
    early = [f"{n}{s:05d}" for n, s in written if s <= chk_step]
    if early or not written:
        fail(f"the restart wrote {early or 'nothing'}: nothing of step {chk_step} or before")

    # Plotfiles: those of the reference, every cell to the last bit.
    theirs = compare_plotfiles(ref, same_set=False)
    last = read_plotfile(ref / theirs[-1])[0]["step"]
    if last <= max(chk_step, first_last):
        fail(f"the last plotfile, of step {last}, is not past the first part (step {first_last})")

    # Gauges: one cell each, no sum over the ranks.
    for g in gauges:
        if not Path(g).is_file() or not filecmp.cmp(g, ref / g, shallow=False):
            fail(f"{g} differs from its copy in {ref}")
    print(f"PASS: {len(gauges)} gauge file(s) identical to those in {ref}")

    # The water balance.
    log = "water_balance.txt"
    if rtol == 0.0:
        if not filecmp.cmp(log, ref / log, shallow=False):
            fail(f"{log} differs from its copy in {ref}")
        print(f"PASS: {log} identical to its copy in {ref}")
        return
    a, b = data_lines(log), data_lines(ref / log)
    if [r[:2] for r in a] != [r[:2] for r in b] or any(len(x) != len(y) for x, y in zip(a, b)):
        fail(f"{log}: other steps, times or columns than in {ref}")
    worst = 0.0
    ncol = len(b[0])
    # Columns 2 to the last but one are volumes; the last is the residual,
    # compared against the largest volume.
    scale = max(abs(float(r[c])) for r in b for c in range(2, ncol - 1))
    for c in range(2, ncol):
        col_scale = scale if c == ncol - 1 else max(abs(float(r[c])) for r in b)
        for x, y in zip(a, b):
            d = abs(float(x[c]) - float(y[c]))
            if d > rtol * col_scale:
                fail(f"{log} step {x[0]} column {c}: {x[c]} here, {y[c]} in {ref}")
            if col_scale > 0.0:
                worst = max(worst, d / col_scale)
    print(f"PASS: {log}: {len(a)} lines, largest difference {worst:.2e} of its column's scale"
          f" (allowed {rtol:g})")


if __name__ == "__main__":
    main()
