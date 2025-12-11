#!/usr/bin/env python3
"""Compare two hparam_run_results.csv files and highlight ERR deltas."""
import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def canonical_param(row: dict) -> str:
    label = row.get("parameter_label", "") or ""
    name = row.get("parameter_name", "") or ""
    if label.startswith("baseline_"):
        return name or label.replace("baseline_", "", 1)
    if label.startswith("general_"):
        parts = label.split("_", 2)
        if len(parts) >= 2:
            return parts[1]
    primary = label.split("_", 1)[0]
    return primary or name


def load_rows(csv_path: Path) -> List[Tuple[Tuple[str, str, str, str], float, dict]]:
    rows: List[Tuple[Tuple[str, str, str, str], float, dict]] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            err = row.get("err")
            if not err:
                continue
            try:
                err_val = float(err)
            except ValueError:
                continue
            key = (
                row.get("dataset", ""),
                row.get("subset", ""),
                row.get("truth_mode", ""),
                canonical_param(row),
            )
            rows.append((key, err_val, row))
    return rows


def format_value(row: Optional[dict]) -> str:
    if not row:
        return "-"
    val = row.get("parameter_value") or row.get("parameter_label", "-")
    return val


def format_err(entry: Optional[Tuple[float, dict]]) -> str:
    if not entry:
        return "-"
    return f"{entry[0]:.6f}"


def main():
    parser = argparse.ArgumentParser(description="Compare ERR deltas across two hparam CSVs.")
    parser.add_argument("base", type=Path, help="Reference CSV (e.g., previous sweep).")
    parser.add_argument("candidate", type=Path, help="New CSV to compare against the base.")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit printed rows (0 = show all, default).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Optional path to write the comparison table as CSV.",
    )
    args = parser.parse_args()

    base_rows = load_rows(args.base)
    cand_rows = load_rows(args.candidate)
    base_map: Dict[Tuple[str, str, str, str], Tuple[float, dict]] = {
        key: (err, row) for key, err, row in base_rows
    }

    comparison = []
    for key, cand_err, cand_row in cand_rows:
        base_entry = base_map.get(key)
        delta = None
        if base_entry:
            delta = cand_err - base_entry[0]
        comparison.append(
            {
                "dataset": key[0],
                "subset": key[1],
                "truth_mode": key[2],
                "parameter_label": cand_row.get("parameter_label", ""),
                "parameter_name": canonical_param(cand_row),
                "base_value": format_value(base_entry[1] if base_entry else None),
                "base_err": base_entry[0] if base_entry else None,
                "candidate_value": format_value(cand_row),
                "candidate_err": cand_err,
                "delta_err": delta,
            }
        )

    # include baseline-only rows for visibility
    candidate_keys = {key for key, _, _ in cand_rows}
    for key, base_entry in base_map.items():
        if key in candidate_keys:
            continue
        comparison.append(
            {
                "dataset": key[0],
                "subset": key[1],
                "truth_mode": key[2],
                "parameter_label": base_entry[1].get("parameter_label", ""),
                "parameter_name": canonical_param(base_entry[1]),
                "base_value": format_value(base_entry[1]),
                "base_err": base_entry[0],
                "candidate_value": "",
                "candidate_err": None,
                "delta_err": None,
            }
        )

    comparison.sort(
        key=lambda row: (
            row["delta_err"] if row["delta_err"] is not None else float("inf"),
            row["dataset"],
            row["subset"],
            row["parameter_label"],
        )
    )

    header = (
        f"{'dataset':<18} {'subset':<14} {'truth':<15} {'param_label':<32} "
        f"{'base(val/err)':<24} {'cand(val/err)':<24} {'delta':<10}"
    )
    print(header)
    print("-" * len(header))

    rows_to_print = comparison if args.limit <= 0 else comparison[: args.limit]
    for row in rows_to_print:
        base_str = (
            f"{row['base_value']}/{row['base_err']:.6f}"
            if row["base_err"] is not None
            else "-"
        )
        if row["candidate_err"] is not None:
            cand_str = f"{row['candidate_value']}/{row['candidate_err']:.6f}"
        else:
            cand_str = "-"
        delta_str = f"{row['delta_err']:.6f}" if row["delta_err"] is not None else "-"
        print(
            f"{row['dataset']:<18} {row['subset']:<14} {row['truth_mode']:<15} {row['parameter_label']:<32} "
            f"{base_str:<24} {cand_str:<24} {delta_str:<10}"
        )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "dataset",
                    "subset",
                    "truth_mode",
                    "parameter_label",
                    "parameter_name",
                    "base_value",
                    "base_err",
                    "candidate_value",
                    "candidate_err",
                    "delta_err",
                ],
            )
            writer.writeheader()
            for row in comparison:
                writer.writerow(row)


if __name__ == "__main__":
    main()
