#!/usr/bin/env python3
"""
Validate that the final column of gzipped TSV files only contains 0/1 labels.
Only the first --limit non-empty lines are checked (default 100000).
"""

import argparse
import gzip
from pathlib import Path
from typing import Tuple


def validate_file(path: Path, limit: int) -> Tuple[bool, int, str, int]:
    checked = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for lineno, raw in enumerate(handle, 1):
            stripped = raw.replace("\x00", "").strip()
            if not stripped:
                continue
            checked += 1
            label = stripped.split("\t")[-1].strip()
            if label not in {"0", "1"}:
                return False, checked, label, lineno
            if checked >= limit:
                break
    return True, checked, "", lineno if 'lineno' in locals() else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="Paths to gzipped TSV files to validate.")
    parser.add_argument("--limit", type=int, default=100_000, help="Lines to check per file (non-empty).")
    args = parser.parse_args()

    exit_code = 0
    for file_arg in args.files:
        path = Path(file_arg)
        if not path.exists():
            print(f"[ERROR] {path} does not exist")
            exit_code = 1
            continue
        ok, checked, label, lineno = validate_file(path, args.limit)
        if ok:
            print(f"[OK] {path} — checked {checked} non-empty lines (limit {args.limit})")
        else:
            print(f"[FAIL] {path} — invalid label '{label}' at logical row {checked} (file line {lineno})")
            exit_code = 2
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
