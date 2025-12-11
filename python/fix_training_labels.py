#!/usr/bin/env python3
import gzip
import os
from pathlib import Path

TRAIN_ROOT = Path(__file__).resolve().parents[1] / "training_datasets"


def normalize_label(raw: str) -> str:
    label = raw.replace("\x00", "").strip()
    if not label:
        normalize_label.missing += 1  # type: ignore[attr-defined]
        return "0"
    unique = set(label)
    if unique <= {"0"}:
        return "0"
    if unique <= {"1"}:
        return "1"
    raise ValueError(f"Unexpected label '{raw}' -> '{label}'")


def should_log(line_number: int) -> bool:
    if line_number == 1:
        return True
    if line_number <= 100_000:
        return line_number % 1_000 == 0
    return line_number % 25_000 == 0


def fix_file(path: Path, avg_lines: int) -> bool:
    changed = False
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(path, "rt") as src, gzip.open(tmp_path, "wt") as dst:
        for idx, line in enumerate(src, start=1):
            clean = line.rstrip("\n")
            parts = clean.split("\t")
            label = normalize_label(parts[-1])
            if parts[-1] != label:
                changed = True
                parts[-1] = label
                if should_log(idx):
                    remaining_pct = 100.0 if avg_lines <= 0 else max(avg_lines - idx, 0) / avg_lines * 100.0
                    print(
                        f"    changed line {idx:,} (~{remaining_pct:6.2f}% remaining assuming {avg_lines:,} rows) in {path}"
                    )
            dst.write("\t".join(parts))
            dst.write("\n")
    os.replace(tmp_path, path)
    return changed


def estimate_lines(path: Path, sample_rows: int = 100_000) -> int:
    total_bytes = path.stat().st_size
    read_rows = 0
    read_bytes = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            read_rows += 1
            read_bytes += len(line.encode("utf-8"))
            if read_rows >= sample_rows:
                break
    if read_rows == 0:
        return 0
    bytes_per_row = read_bytes / read_rows
    if bytes_per_row == 0:
        return 0
    estimate = int(total_bytes / bytes_per_row)
    return max(estimate, read_rows)


def main():
    files = sorted(TRAIN_ROOT.rglob("*.tsv.gz"))
    if not files:
        print(f"No gzipped TSVs found under {TRAIN_ROOT}")
        return
    fixed = 0
    normalize_label.missing = 0  # type: ignore[attr-defined]
    for path in files:
        avg_lines = estimate_lines(path)
        print(f"[start] {path} (estimated rows={avg_lines:,})")
        changed = fix_file(path, avg_lines)
        if changed:
            fixed += 1
            print(f"[done ] normalized labels in {path}")
        else:
            print(f"[done ] labels already clean in {path}")
    print(
        f"Done. Files updated: {fixed} / {len(files)}; empty labels coerced: {normalize_label.missing}"
    )


if __name__ == "__main__":
    main()
