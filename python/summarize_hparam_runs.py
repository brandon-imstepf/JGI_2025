#!/usr/bin/env python3
"""
Summarize training logs into CSV tables for spreadsheet analysis.

This script supports two primary use cases:
  1. **Sweep analysis** – Parse `logs/hparam_sweep/<dataset_subset>/` to collate best ERR/WER/etc.
  2. **Baseline analysis** – With `--glob`, parse standalone logs (e.g., `logs/hparam_baseline/*.log`).

Data extracted per log:
  * Full command-line arguments (the script derives parameter names/values from the log header or file path).
  * The `Best: ERR=` / `WER=` / `FPR=` / `FNR=` / `CTF=` metrics, when present.
  * Dataset/subset names are inferred from the log location; `cohort` and `truth_mode` are derived from the dataset name.

Outputs:
  * `hparam_error_summary.csv` – best ERR per dataset/subset/parameter.
  * `hparam_run_results.csv` – one row per run (even failed ones).
  * `hparam_success_summary.csv` – counts of successful runs per dataset/subset.
  * Optional per-subset CSVs via `--per-subset`.

Usage examples:
  * Sweeps (default root = logs/hparam_sweep):
      `python3 python/summarize_hparam_runs.py --per-subset`
  * Baselines:
      `python3 python/summarize_hparam_runs.py --glob 'logs/hparam_baseline/*.log' --out-dir logs/hparam_baseline`
"""

import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

PROJECT_ROOT = Path("/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf")
LOG_ROOT = PROJECT_ROOT / "logs" / "hparam_sweep"
KNOWN_SUBSETS = ("full", "floats_only", "onehot_only", "floats_plus_onehot_no_middle")

BEST_PATTERN = re.compile(
    r"Best:\s+ERR=\s+(?P<err>[0-9.]+)\s+WER=\s+(?P<wer>[0-9.]+)\s+FPR=\s+(?P<fpr>[0-9.]+)\s+"
    r"FNR=\s+(?P<fnr>[0-9.]+)\s+CTF=\s+(?P<ctf>[0-9.\-]+)"
)
CMD_PATTERN = re.compile(r"([A-Za-z_]+)=([0-9.eE+\-]+)")
BASELINE_PATTERN = re.compile(
    r"baseline_(?P<dataset>.+?)_(?P<subset>full|floats_only|onehot_only|floats_plus_onehot_no_middle)"
    r"(?:_(?P<truth>ncbi|ncbi_intersect))?$"
)


def parse_command_line(line: str) -> Dict[str, str]:
    return {k: v for k, v in CMD_PATTERN.findall(line)}


def split_combo(name: str) -> Tuple[str, str]:
    for subset in KNOWN_SUBSETS:
        suffix = f"_{subset}"
        if name.endswith(suffix):
            dataset = name[: -len(suffix)]
            if dataset.endswith("_"):
                dataset = dataset[:-1]
            return dataset, subset
    return name, "unknown"


def classify_dataset(dataset: str) -> Tuple[str, str]:
    ds_lower = dataset.lower()
    if "training" in ds_lower:
        cohort = "october"
    elif "sorted" in ds_lower or "10-23" in ds_lower:
        cohort = "july"
    else:
        cohort = "unknown"
    truth_mode = "ncbi_intersect" if "intersect" in ds_lower else "ncbi"
    return cohort, truth_mode


def extract_param_value(label: str, cmd_map: Dict[str, str]) -> Tuple[str, str]:
    parts = label.split("_")
    prefix = parts[0]
    suffix = parts[-1]
    core = "_".join(parts[1:-1]) if len(parts) > 2 else (parts[1] if len(parts) > 1 else parts[0])

    if prefix == "triage" or core == "triage":
        ptriage = cmd_map.get("ptriage")
        ntriage = cmd_map.get("ntriage")
        val = f"ptriage={ptriage},ntriage={ntriage}" if ptriage and ntriage else suffix
        return "triage", val

    for pre in ("general_", "mse_", "legacy_", "wbce_", "tversky_"):
        if core.startswith(pre):
            core = core[len(pre) :]
    if not core and len(parts) > 1:
        core = parts[1]

    if core in cmd_map:
        return core, cmd_map[core]
    return core or parts[0], suffix


def parse_log(path: Path) -> Dict:
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    cmd_map = parse_command_line(lines[0]) if lines else {}
    match = BEST_PATTERN.search(text)
    status = "SUCCESS" if match else "FAIL"
    metrics = {k: float(v) for k, v in match.groupdict().items()} if match else {}
    dataset, subset = split_combo(path.parent.name)
    cohort, truth_mode = classify_dataset(dataset)
    baseline_match = BASELINE_PATTERN.match(path.stem)
    if baseline_match:
        dataset = baseline_match.group("dataset")
        subset = baseline_match.group("subset")
        truth_mode = baseline_match.group("truth") or truth_mode
        cohort, _ = classify_dataset(dataset)
    elif subset == "unknown":
        for candidate in KNOWN_SUBSETS:
            if candidate in path.stem:
                subset = candidate
                break
    param_name, param_value = extract_param_value(path.stem, cmd_map)
    return {
        "dataset": dataset,
        "subset": subset,
        "cohort": cohort,
        "truth_mode": truth_mode,
        "parameter": path.stem,
        "parameter_name": param_name,
        "parameter_value": param_value,
        "status": status,
        **metrics,
    }


def collect_records(dirs: Iterable[Path]) -> List[Dict]:
    records: List[Dict] = []
    for combo_dir in dirs:
        if not combo_dir.is_dir():
            continue
        for log_path in combo_dir.glob("*.log"):
            records.append(parse_log(log_path))
    return records


def summarize(records: Sequence[Dict]) -> Tuple[List[Dict], List[Tuple[str, str, int]]]:
    best_rows: Dict[Tuple[str, str, str, str], Dict] = {}
    combo_counts: Dict[Tuple[str, str], int] = {}
    run_rows: List[Dict] = []
    for rec in records:
        key_combo = (rec["dataset"], rec["subset"])
        status = rec["status"]
        if status == "SUCCESS":
            combo_counts[key_combo] = combo_counts.get(key_combo, 0) + 1
            key = (rec["dataset"], rec["subset"], rec["parameter_name"], rec["parameter_value"])
            current = best_rows.get(key)
            if current is None or rec["err"] < current["err"]:
                best_rows[key] = {
                    "dataset": rec["dataset"],
                    "subset": rec["subset"],
                    "cohort": rec["cohort"],
                    "truth_mode": rec["truth_mode"],
                    "parameter_name": rec["parameter_name"],
                    "parameter_value": rec["parameter_value"],
                    "err": rec["err"],
                    "wer": rec.get("wer"),
                    "fpr": rec.get("fpr"),
                    "fnr": rec.get("fnr"),
                    "ctf": rec.get("ctf"),
                }
            run_rows.append(
                {
                    "dataset": rec["dataset"],
                    "subset": rec["subset"],
                    "cohort": rec["cohort"],
                    "truth_mode": rec["truth_mode"],
                    "parameter_label": rec["parameter"],
                    "parameter_name": rec["parameter_name"],
                    "parameter_value": rec["parameter_value"],
                    "err": rec["err"],
                    "wer": rec.get("wer"),
                    "fpr": rec.get("fpr"),
                    "fnr": rec.get("fnr"),
                    "ctf": rec.get("ctf"),
                }
            )
        else:
            run_rows.append(
                {
                    "dataset": rec["dataset"],
                    "subset": rec["subset"],
                    "cohort": rec["cohort"],
                    "truth_mode": rec["truth_mode"],
                    "parameter_label": rec["parameter"],
                    "parameter_name": rec["parameter_name"],
                    "parameter_value": rec["parameter_value"],
                    "err": "",
                    "wer": "",
                    "fpr": "",
                    "fnr": "",
                    "ctf": "",
                }
            )
    combo_summary = [(ds, subset, count) for (ds, subset), count in combo_counts.items()]
    combo_summary.sort(key=lambda x: (-x[2], x[0], x[1]))
    best_rows_list = list(best_rows.values())
    best_rows_list.sort(key=lambda x: (x["dataset"], x["subset"], x["parameter_name"], x["parameter_value"]))
    run_rows.sort(key=lambda x: (x["dataset"], x["subset"], x["parameter_label"]))
    return best_rows_list, combo_summary, run_rows


def write_csv(path: Path, rows: List[Dict], headers: List[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_combo_csv(path: Path, rows: List[Tuple[str, str, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["dataset", "subset", "success_count"])
        for dataset, subset, count in rows:
            writer.writerow([dataset, subset, count])


def main():
    parser = argparse.ArgumentParser(description="Summarize hparam sweep logs into CSV.")
    parser.add_argument("--log-root", default=str(LOG_ROOT), help="Root containing log subfolders (default: logs/hparam_sweep)")
    parser.add_argument("--log-dirs", nargs="*", help="Explicit log directories to process (overrides --log-root).")
    parser.add_argument("--glob", help="Glob pattern for log files (e.g., 'logs/hparam_baseline/*.log').")
    parser.add_argument("--out-dir", default=str(LOG_ROOT), help="Directory to write CSVs (default: log root)")
    parser.add_argument("--per-subset", action="store_true", help="Also write per-subset CSVs")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.glob:
        log_paths = [Path(p) for p in Path().glob(args.glob)]
        if not log_paths:
            print(f"[INFO] Glob '{args.glob}' matched no log files.")
            return
        records = [parse_log(p) for p in log_paths]
        combos = []
    elif args.log_dirs:
        combos = [Path(p) for p in args.log_dirs]
        records = collect_records(combos)
    else:
        log_root = Path(args.log_root)
        combos = sorted(p for p in log_root.iterdir() if p.is_dir())
        if not combos:
            print(f"[INFO] No log directories found under {log_root}")
            return
        records = collect_records(combos)
    if not records:
        print("[INFO] No log files found.")
        return

    best_rows, combo_summary, run_rows = summarize(records)
    headers = ["dataset", "subset", "cohort", "truth_mode", "parameter_name", "parameter_value", "err", "wer", "fpr", "fnr", "ctf"]
    best_path = out_dir / "hparam_error_summary.csv"
    write_csv(best_path, best_rows, headers)
    combo_path = out_dir / "hparam_success_summary.csv"
    write_combo_csv(combo_path, combo_summary)
    print(f"Wrote {best_path}")
    print(f"Wrote {combo_path}")
    run_path = out_dir / "hparam_run_results.csv"
    run_headers = ["dataset", "subset", "cohort", "truth_mode", "parameter_label", "parameter_name", "parameter_value", "err", "wer", "fpr", "fnr", "ctf"]
    write_csv(run_path, run_rows, run_headers)
    print(f"Wrote {run_path}")

    if args.per_subset:
        by_subset = defaultdict(list)
        for row in best_rows:
            by_subset[row["subset"]].append(row)
        for subset, rows in by_subset.items():
            subset_path = out_dir / f"hparam_error_summary_{subset}.csv"
            write_csv(subset_path, rows, headers)
            print(f"Wrote {subset_path}")


if __name__ == "__main__":
    main()
