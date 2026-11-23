#!/usr/bin/env python3
"""Run one test case in its own directory, then check it.

    run_case.py --dir <dir> [--expect-output <regex>] [--reject-output <regex>]
                --run <command...> --check <command...>
    run_case.py --dir <dir> --expect-abort <regex> --run <command...>

The directory is removed and made again, so a test never reads the output of an
earlier run. With --check, the run must exit with 0 and the check command, run
in the same directory, must exit with 0. With --expect-abort, the run must exit
with a non-zero code and its output must match the regular expression. With
--expect-output, the output of a run that is checked must also match that
regular expression; with --reject-output it must not match that one. A run that is checked and leaves water_balance.txt must
also close its water balance to round-off (Tests/check_balance.py).
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_balance  # noqa: E402


def split_args(argv):
    """Split argv at the --dir, --expect-abort, --expect-output, --reject-output, --run and
    --check markers."""
    out = {"--dir": [], "--expect-abort": [], "--expect-output": [], "--reject-output": [],
           "--run": [], "--check": []}
    key = None
    for a in argv:
        if a in out and not (key == "--check" and a != "--run"):
            # Once --check starts, everything after it belongs to the check
            # command, which has options of its own.
            key = a
            continue
        if key is None:
            sys.exit(f"run_case.py: unexpected argument {a!r}")
        out[key].append(a)
    return out


def main():
    # The run's output is printed again below; a byte the console cannot show
    # must not fail the test.
    sys.stdout.reconfigure(errors="backslashreplace")
    args = split_args(sys.argv[1:])
    if len(args["--dir"]) != 1 or not args["--run"]:
        sys.exit("run_case.py: --dir <dir> and --run <command...> are required")
    expect = args["--expect-abort"]
    if len(expect) > 1 or (bool(expect) == bool(args["--check"])):
        sys.exit("run_case.py: give exactly one of --check <command...> and --expect-abort <regex>")
    if expect and not expect[0].strip():
        sys.exit("run_case.py: the --expect-abort expression is empty and would match any output")
    output = args["--expect-output"]
    if len(output) > 1 or (output and (expect or not output[0].strip())):
        sys.exit("run_case.py: --expect-output takes one expression that is not empty, and goes with --check")
    reject = args["--reject-output"]
    if len(reject) > 1 or (reject and (expect or not reject[0].strip())):
        sys.exit("run_case.py: --reject-output takes one expression that is not empty, and goes with --check")

    run_dir = Path(args["--dir"][0])
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    print("run:", " ".join(args["--run"]), flush=True)
    run = subprocess.run(args["--run"], cwd=run_dir, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, encoding="utf-8", errors="backslashreplace")
    print(run.stdout, flush=True)

    if expect:
        if run.returncode == 0:
            sys.exit("FAIL: the run was expected to abort but exited with 0")
        if not re.search(expect[0], run.stdout):
            sys.exit(f"FAIL: the run aborted without the expected message /{expect[0]}/")
        print(f"PASS: aborted (exit code {run.returncode}) with the expected message")
        return

    if run.returncode != 0:
        sys.exit(f"FAIL: the run exited with {run.returncode}")
    if output and not re.search(output[0], run.stdout):
        sys.exit(f"FAIL: the output of the run does not match /{output[0]}/")
    if reject and re.search(reject[0], run.stdout):
        sys.exit(f"FAIL: the output of the run matches /{reject[0]}/")
    print("check:", " ".join(args["--check"]), flush=True)
    check = subprocess.run(args["--check"], cwd=run_dir)
    if check.returncode != 0:
        sys.exit(f"FAIL: the check exited with {check.returncode}")

    # Every run that keeps a water balance (flume.log_int > 0, the default)
    # must close it to round-off at every line (Tests/check_balance.py).
    log = run_dir / "water_balance.txt"
    if log.exists():
        ok, msg, _ = check_balance.check(log)
        print(f"balance: {msg}", flush=True)
        if not ok:
            sys.exit(f"FAIL: water balance: {msg}")


if __name__ == "__main__":
    main()
