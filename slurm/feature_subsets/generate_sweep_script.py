#!/usr/bin/env python3
import gzip
import shlex
from pathlib import Path
from typing import List

PROJECT_ROOT = Path("/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf")
TRAIN_SCRIPT = PROJECT_ROOT / "bbmap" / "train.sh"
TRAINING_ROOT = PROJECT_ROOT / "training_datasets"
LOG_ROOT = PROJECT_ROOT / "logs" / "hparam_sweep"
NET_ROOT = PROJECT_ROOT / "comparison_models" / "hparam_sweep_nets"
SCRIPT_DIR = Path(__file__).parent

DATASET_BASES = [
    {
        "name": "ncbi_sorted",
        "train_dir": "NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets",
        "prefix": "NCBI_Datasets_10-23_SORTED.ncbi_truth",
    },
    {
        "name": "ncbi_sorted_intersect",
        "train_dir": "NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052_subsets",
        "prefix": "NCBI_Datasets_10-23_SORTED.intersect_truth",
    },
    {
        "name": "ncbi_training",
        "train_dir": "NCBI_Datasets_Training_ncbi_truth_slurm_19888056_subsets",
        "prefix": "NCBI_Datasets_Training.ncbi_truth",
    },
    {
        "name": "ncbi_training_intersect",
        "train_dir": "NCBI_Datasets_Training_intersect_truth_slurm_19888054_subsets",
        "prefix": "NCBI_Datasets_Training.intersect_truth",
    },
]

SUBSETS = ["full", "floats_only", "onehot_only", "floats_plus_onehot_no_middle"]

VALIDATION_CONFIG = {
    "dir": "Validation_NCBI_10_ncbi_truth_20251120_122824_subsets",
    "prefix": "Validation_NCBI_10.ncbi_truth",
}

FALLBACK_INPUTS = 356

GENERAL_PARAMS = {
    "alpha": 0.08,
    "balance": 0.2,
    "density": 1.0,
    "fpb": 0.08,
    "setsize": 60000,
    "vfraction": 0.1,
    "cutoffbackprop": 0.5,
    "spread": 0.05,
    "pem": 1.0,
    "nem": 1.0,
    "fpem": 10.5,
    "fnem": 10.5,
    "epem": 0.2,
    "enem": 0.2,
    "epm": 0.2,
    "anneal": 0.003,
    "annealprob": 0.225,
}

ACTIVATION_PARAMS = {
    "sig": 0.6,
    "tanh": 0.4,
    "rslog": 0.02,
    "msig": 0.02,
}

LEGACY_LOSS_PARAMS = {
    "positiveerrormult": 1.0,
    "negativeerrormult": 1.0,
    "falsepositiveerrormult": 10.5,
    "falsenegativeerrormult": 10.5,
    "excesspositiveerrormult": 0.2,
    "excessnegativeerrormult": 0.2,
}

LEGACY_OFFSET_PARAMS = {
    "falsepositeverrorincr": 0.05,
    "falsenegativeerrorincr": 0.05,
}

TRIAGE_VALUES = [0.008, 0.01, 0.012]

WBCE_PARAMS = {
    "lossposweight": 1.0,
    "lossnegweight": 1.0,
}

TVERSKY_PARAMS = {
    "lossalpha": 0.3,
    "lossbeta": 0.7,
}

SBATCH_ARGS = [
    "-A",
    "grp-org-gt-ga",
    "-q",
    "jgi_normal",
    "-t",
    "08:00:00",
    "-c",
    "8",
    "--mem=64G",
]


def infer_input_count(tsv_path: Path, fallback: int) -> int:
    if not tsv_path.exists():
        return fallback
    opener = gzip.open if tsv_path.suffix.endswith(".gz") else open
    with opener(tsv_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#dims"):
                parts = line.replace("\t", " ").split()
                if len(parts) >= 2:
                    try:
                        return int(float(parts[1]))
                    except ValueError:
                        pass
            else:
                return max(1, line.count("\t"))
    return fallback


def format_value(value: float) -> str:
    if abs(value) < 1e-12:
        return "0"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    magnitude = abs(value)
    if magnitude >= 1000:
        return f"{value:.0f}"
    if magnitude >= 100:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    if magnitude >= 10:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    if magnitude >= 1:
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return f"{value:.4f}".rstrip("0").rstrip(".")


def format_arg(arg: str) -> str:
    if "=" not in arg:
        return arg
    key, val = arg.split("=", 1)
    try:
        num = float(val)
        return f"{key}={format_value(num)}"
    except ValueError:
        return arg


def encode_value_token(token: str) -> str:
    return (
        token.replace("-", "neg")
        .replace(".", "p")
        .replace("=", "")
        .replace(",", "_")
        .replace("/", "_")
    )


def format_command(base_cmd: str, log_dir: Path, extra_args, log_name: str) -> str:
    pieces = [base_cmd.format(log_name=log_name)]
    pieces.extend(format_arg(arg) for arg in extra_args)
    cmd = " ".join(pieces)
    cmd += f" > {log_dir / (log_name + '.log')} 2>&1"
    return cmd


def build_base_cmd(train_path: Path, validate_path: Path, net_dir: Path, input_count: int) -> str:
    mindims = f"{input_count},128,96,64,32,1"
    maxdims = f"{input_count},256,192,128,64,1"
    return (
        f"{TRAIN_SCRIPT} -Xmx64g seed=12345 "
        f"in={train_path} "
        f"validate={validate_path} "
        f"mindims={mindims} "
        f"maxdims={maxdims} "
        f"nets=1 cycles=1 batches=200000 "
        f"out={net_dir}/{{log_name}}.bbnet"
    )


def build_job_lines(base_cmd: str, log_dir: Path) -> List[str]:
    lines: List[str] = []

    def add_job(name: str, extra):
        cmd = format_command(base_cmd, log_dir, extra, name)
        sbatch_cmd = " ".join(
            [
                "sbatch",
                *SBATCH_ARGS,
                f"--output={log_dir / (name + '.slurm.out')}",
                f"--error={log_dir / (name + '.slurm.err')}",
                f"--wrap={shlex.quote(cmd)}",
            ]
        )
        lines.append(sbatch_cmd + "\n")

    for param, default in GENERAL_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            if param == "density":
                value = min(1.0, value)
            value_label = encode_value_token(format_value(value))
            add_job(
                f"{param}_{value_label}",
                [
                    f"{param}={format_value(value)}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
            )

    activation_total = sum(ACTIVATION_PARAMS.values())
    for target, base_value in ACTIVATION_PARAMS.items():
        remaining_total = activation_total - base_value
        if remaining_total <= 0:
            continue
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            new_target = base_value * delta
            remaining_target = activation_total - new_target
            scale = remaining_target / remaining_total
            values = {}
            for name, val in ACTIVATION_PARAMS.items():
                values[name] = new_target if name == target else val * scale
            add_job(
                "activation_"
                + target
                + "_"
                + "_".join(f"{k}{encode_value_token(format_value(v))}" for k, v in values.items()),
                [
                    *(f"{k}={format_value(v)}" for k, v in values.items()),
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
            )

    for param, default in LEGACY_LOSS_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            add_job(
                f"{param}_{encode_value_token(format_value(value))}",
                [
                    f"{param}={format_value(value)}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
            )

    for val in TRIAGE_VALUES:
        encoded = encode_value_token(format_value(val))
        add_job(
            f"triage_ptr{encoded}_ntr{encoded}",
            [
                "loss=mse",
                "losslegacyweighting=t",
                f"ptriage={format_value(val)}",
                f"ntriage={format_value(val)}",
            ],
        )

    for param, step in LEGACY_OFFSET_PARAMS.items():
        for direction, tag in ((1, "high"), (-1, "low")):
            value = direction * step
            add_job(
                f"{param}_{encode_value_token(format_value(value))}",
                [
                    f"{param}={format_value(value)}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
            )

    for param, default in WBCE_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            add_job(
                f"{param}_{encode_value_token(format_value(value))}",
                [
                    "loss=wbce",
                    f"{param}={format_value(value)}",
                    "losslegacyweighting=f",
                ],
            )

    for param, default in TVERSKY_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            add_job(
                f"{param}_{encode_value_token(format_value(value))}",
                [
                    "loss=tversky",
                    f"{param}={format_value(value)}",
                    "losslegacyweighting=f",
                ],
            )

    return lines


def collect_dataset_configs():
    configs = []
    missing = []
    validation_dir = TRAINING_ROOT / VALIDATION_CONFIG["dir"]
    for base in DATASET_BASES:
        train_root = TRAINING_ROOT / base["train_dir"]
        for subset in SUBSETS:
            train_file = train_root / subset / f"{base['prefix']}.all.subset_{subset}.tsv.gz"
            if not train_file.exists():
                missing.append(f"{base['name']}:{subset} (training missing)")
                continue
            validate_file = (
                validation_dir / subset / f"{VALIDATION_CONFIG['prefix']}.all.subset_{subset}.tsv.gz"
            )
            if not validate_file.exists():
                missing.append(f"{base['name']}:{subset} (validation missing)")
                continue
            input_count = infer_input_count(train_file, FALLBACK_INPUTS)
            configs.append(
                {
                    "name": f"{base['name']}_{subset}",
                    "train": train_file,
                    "validate": validate_file,
                    "inputs": input_count,
                }
            )
    return configs, missing


def main():
    configs, missing = collect_dataset_configs()
    if not configs:
        print("No dataset/subset combinations are ready. Nothing to do.")
        if missing:
            print("Missing combos:")
            for item in missing:
                print(f"  - {item}")
        return

    for cfg in configs:
        log_dir = LOG_ROOT / cfg["name"]
        net_dir = NET_ROOT / cfg["name"]
        base_cmd = build_base_cmd(cfg["train"], cfg["validate"], net_dir, cfg["inputs"])
        lines = build_job_lines(base_cmd, log_dir)
        script_path = SCRIPT_DIR / f"run_hparam_sweep_{cfg['name']}.sh"
        with open(script_path, "w") as out:
            out.write("#!/bin/bash\n")
            out.write(f"mkdir -p {log_dir}\n")
            out.write(f"mkdir -p {net_dir}\n")
            for line in lines:
                out.write(line)
                out.write("sleep 2\n")
        script_path.chmod(0o755)

    print(f"Generated {len(configs)} sweep scripts in {SCRIPT_DIR}")
    if missing:
        print("Skipped combinations (files not ready):")
        for item in missing:
            print(f"  - {item}")


if __name__ == "__main__":
    main()
