#!/usr/bin/env python3
"""Check a run with flume.log_int = 0 against the same run with the log.

    check_no_log.py <prefix> <reference directory>

The run must leave no water_balance.txt and no gauge file, and its plotfiles
<prefix>NNNNN must be byte for byte those of the reference run
(Tests/compare_plotfiles.py): the log reads the state and changes nothing.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    left = [p.name for p in Path(".").iterdir()
            if p.name == "water_balance.txt" or p.name.startswith("gauge_")]
    if left:
        sys.exit(f"FAIL: the run with flume.log_int = 0 wrote {sorted(left)}")
    ref = Path(sys.argv[2]) / "water_balance.txt"
    if not ref.exists():
        sys.exit(f"FAIL: the reference run wrote no {ref}")
    r = subprocess.run([sys.executable, str(HERE.parent / "compare_plotfiles.py"),
                        sys.argv[1], sys.argv[2]])
    if r.returncode != 0:
        sys.exit("FAIL: the plotfiles differ from those of the run with the log")
    print("PASS")


if __name__ == "__main__":
    main()
