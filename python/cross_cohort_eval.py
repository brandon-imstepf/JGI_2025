#!/usr/bin/env python3
"""
Cross-cohort evaluator for BBNet baselines.

Reads the baseline summary CSV, selects the top-N (lowest ERR) networks,
and evaluates each model on the opposite cohort's dataset using train.sh
in evaluate-only mode. Results are written to a consolidated CSV along
with per-model log files.
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Dict, List
import re

PROJECT_ROOT = Path("/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf")
TRAIN_SH = PROJECT_ROOT / "bbmap" / "train.sh"
DEFAULT_NET_ROOT = PROJECT_ROOT / "networks" / "hparam_baseline_nets"
RESULTS_CSV = PROJECT_ROOT / "logs" / "hparam_baseline" / "cross_cohort_eval.csv"
LOG_DIR = PROJECT_ROOT / "logs" / "hparam_baseline" / "cross_cohort_logs"

DATASET_INFO: Dict[str, Dict[str, Dict[str, str]]] = {
    "july": {
        "ncbi": {
            "dir": "NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055",
            "subsets": "NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets",
            "prefix": "NCBI_Datasets_10-23_SORTED.ncbi_truth",
        },
        "ncbi_intersect": {
            "dir": "NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052",
            "subsets": "NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052_subsets",
            "prefix": "NCBI_Datasets_10-23_SORTED.intersect_truth",
        },
    },
    "october": {
        "ncbi": {
            "dir": "NCBI_Datasets_Training_ncbi_truth_slurm_19888056",
            "subsets": "NCBI_Datasets_Training_ncbi_truth_slurm_19888056_subsets",
            "prefix": "NCBI_Datasets_Training.ncbi_truth",
        },
        "ncbi_intersect": {
            "dir": "NCBI_Datasets_Training_intersect_truth_slurm_19888054",
            "subsets": "NCBI_Datasets_Training_intersect_truth_slurm_19888054_subsets",
            "prefix": "NCBI_Datasets_Training.intersect_truth",
        },
    },
}
class ModelRecord(object):
    def __init__(
        self,
        dataset,
        subset,
        truth_mode,
        err,
        parameter_label,
        parameter_name,
        parameter_value,
        source_cohort,
        net_root,
    ):
        self.dataset = dataset
        self.subset = subset
        self.truth_mode = truth_mode
        self.err = err
        self.parameter_label = parameter_label
        self.parameter_name = parameter_name
        self.parameter_value = parameter_value
        self.source_cohort = source_cohort
        self.net_root = net_root

    def model_path(self) -> Path:
        return self.net_root / self.parameter_name / (self.parameter_label + ".bbnet")


def infer_cohort(dataset_name: str) -> str:
    dlower = dataset_name.lower()
    if "training" in dlower:
        return "october"
    if "sorted" in dlower or "10-23" in dlower:
        return "july"
    raise ValueError("Cannot infer cohort from dataset name '%s'" % dataset_name)


def load_top_models(csv_path: Path, top_n: int, net_root: Path) -> List[ModelRecord]:
    rows: List[ModelRecord] = []
    with csv_path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            err_str = (row.get("err", "") or "").strip()
            if not err_str:
                continue
            rows.append(
                ModelRecord(
                    dataset=row["dataset"],
                    subset=row["subset"],
                    truth_mode=row["truth_mode"],
                    err=float(err_str),
                    parameter_label=row["parameter_label"],
                    parameter_name=row["parameter_name"],
                    parameter_value=row["parameter_value"],
                    source_cohort=infer_cohort(row["dataset"]),
                    net_root=net_root,
                )
            )
    rows.sort(key=lambda r: r.err)
    return rows[:top_n]


def build_subset_file(cohort: str, truth_mode: str, subset: str) -> Path:
    info = DATASET_INFO[cohort][truth_mode]
    subset_dir = PROJECT_ROOT / "training_datasets" / info["subsets"] / subset
    filename = "%s.all.subset_%s.tsv.gz" % (info["prefix"], subset)
    path = subset_dir / filename
    if not path.exists():
        raise FileNotFoundError("Missing subset file: %s" % path)
    return path



def run_evaluation(record: ModelRecord, target_cohort: str, log_dir: Path):
    model_path = record.model_path()
    if not model_path.exists():
        print("[WARN] Model not found: %s" % model_path, file=sys.stderr)
        return None
    target_train = build_subset_file(target_cohort, record.truth_mode, record.subset)
    cmd = [
        str(TRAIN_SH),
        "evaluate=t",
        "netin=%s" % model_path,
        "in=%s" % target_train,
        "validate=%s" % target_train,
    ]
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / ("%s__to_%s.log" % (record.parameter_label, target_cohort))
    print(
        "[INFO] Evaluating %s on %s (%s, %s)"
        % (record.parameter_label, target_cohort, record.truth_mode, record.subset)
    )
    result = subprocess.run(
        cmd, cwd=PROJECT_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    log_path.write_text(combined)
    if result.returncode != 0:
        print("[WARN] Evaluation failed for %s. See %s" % (model_path, log_path), file=sys.stderr)
        return None
    best_line = ""
    for line in combined.splitlines():
        if line.startswith("Best:"):
            best_line = line
            break
    if not best_line:
        print("[WARN] No 'Best:' line found for %s. See %s" % (model_path, log_path), file=sys.stderr)
        return None
    metrics: Dict[str, str] = {}
    pattern = re.compile(r'(ERR|WER|FPR|FNR|CTF)=\s*([0-9.\-]+)')
    for key, value in pattern.findall(best_line):
        metrics[key.lower()] = value
    if "err" not in metrics:
        print("[WARN] Unable to parse metrics from %s." % model_path, file=sys.stderr)
        return None
    return {
        "dataset": record.dataset,
        "subset": record.subset,
        "truth_mode": record.truth_mode,
        "source_cohort": record.source_cohort,
        "target_cohort": target_cohort,
        "model": str(model_path.relative_to(PROJECT_ROOT)),
        "eval_err": metrics.get("err", ""),
        "eval_wer": metrics.get("wer", ""),
        "eval_fpr": metrics.get("fpr", ""),
        "eval_fnr": metrics.get("fnr", ""),
        "eval_ctf": metrics.get("ctf", ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate top-N baseline networks on the opposite cohort."
    )
    parser.add_argument(
        "--summary",
        default=str(PROJECT_ROOT / "logs" / "hparam_baseline" / "hparam_run_results.csv"),
        help="CSV with baseline run results (default: logs/hparam_baseline/hparam_run_results.csv)",
    )
    parser.add_argument("--top", type=int, default=10, help="Number of models to evaluate (default: 10)")
    parser.add_argument("--out", default=str(RESULTS_CSV), help="Output CSV path")
    parser.add_argument("--net-root", default=str(DEFAULT_NET_ROOT), help="Directory containing trained networks")
    args = parser.parse_args()

    summary_csv = Path(args.summary)
    if not summary_csv.exists():
        parser.error("Summary CSV not found: %s" % summary_csv)
    net_root = Path(args.net_root)
    if not net_root.is_absolute():
        net_root = PROJECT_ROOT / net_root
    models = load_top_models(summary_csv, args.top, net_root)
    if not models:
        parser.error("No models found in the summary CSV.")

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    results: List[Dict[str, str]] = []
    for record in models:
        target = "july" if record.source_cohort == "october" else "october"
        try:
            info = run_evaluation(record, target, LOG_DIR)
        except FileNotFoundError as exc:
            print("[WARN] %s" % exc, file=sys.stderr)
            continue
        if info:
            results.append(info)

    if not results:
        print("[WARN] No evaluations succeeded; nothing to write.", file=sys.stderr)
        return
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset",
                "subset",
                "truth_mode",
                "source_cohort",
                "target_cohort",
                "model",
                "eval_err",
                "eval_wer",
                "eval_fpr",
                "eval_fnr",
                "eval_ctf",
            ],
        )
        writer.writeheader()
        writer.writerows(results)
    print("[INFO] Wrote %s (%d rows)" % (out_path, len(results)))


if __name__ == "__main__":
    main()
