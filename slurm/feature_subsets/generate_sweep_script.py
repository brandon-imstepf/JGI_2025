#!/usr/bin/env python3
import shlex
from pathlib import Path

PROJECT_ROOT = Path("/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf")
TRAIN_SCRIPT = PROJECT_ROOT / "bbmap" / "train.sh"
LOG_DIR = PROJECT_ROOT / "logs" / "hparam_sweep"
NET_DIR = PROJECT_ROOT / "comparison_models" / "hparam_sweep_nets"
OUTPUT_SCRIPT = Path(__file__).with_name("run_hparam_sweep.sh")

TRAIN_INPUT = (
    PROJECT_ROOT
    / "training_datasets"
    / "NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets"
    / "full"
    / "NCBI_Datasets_10-23_SORTED.ncbi_truth.all.subset_full.tsv.gz"
)
VALIDATE_INPUT = (
    PROJECT_ROOT
    / "training_datasets"
    / "Validation_NCBI_10_ncbi_truth_20251120_122824_subsets"
    / "full"
    / "Validation_NCBI_10.ncbi_truth.all.subset_full.tsv.gz"
)

BASE_CMD = (
    f"{TRAIN_SCRIPT} -Xmx64g seed=12345 "
    f"in={TRAIN_INPUT} "
    f"validate={VALIDATE_INPUT} "
    f"mindims=356,128,96,64,32,1 "
    f"maxdims=356,256,192,128,64,1 "
    f"nets=2 cycles=1 batches=200000 "
    f"out={NET_DIR}/{{log_name}}.bbnet"
)

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
    "falsepositeverrorincr": 0.0,
    "falsenegativeerrorincr": 0.0,
}

LEGACY_OFFSET_PARAMS = {
    "falsepositeverrorincr": 0.05,
    "falsenegativeerrorincr": 0.05,
}

WBCE_PARAMS = {
    "lossposweight": 1.0,
    "lossnegweight": 1.0,
}

TVERSKY_PARAMS = {
    "lossalpha": 0.3,
    "lossbeta": 0.7,
}

TRIAGE_VALUES = [0.008, 0.01, 0.012]

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


def format_command(extra_args, log_name):
    pieces = [BASE_CMD.format(log_name=log_name)]
    pieces.extend(extra_args)
    cmd = " ".join(pieces)
    cmd += f" > {LOG_DIR / (log_name + '.log')} 2>&1"
    return cmd


def main():
    lines = []
    sbatch_prefix = " ".join(["sbatch", *SBATCH_ARGS])

    # General trainer sweeps (+/-20%)
    for param, default in GENERAL_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            if param == "density":
                value = min(1.0, value)
            name = f"general_{param}_{tag}"
            cmd = format_command(
                [
                    f"{param}={value}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # Activation sweeps with renormalization
    activation_total = sum(ACTIVATION_PARAMS.values())
    for target, base_value in ACTIVATION_PARAMS.items():
        remaining_total = activation_total - base_value
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            new_target = base_value * delta
            values = {}
            if remaining_total <= 0:
                continue
            remaining_target = activation_total - new_target
            scale = remaining_target / remaining_total
            for name, val in ACTIVATION_PARAMS.items():
                if name == target:
                    values[name] = new_target
                else:
                    values[name] = val * scale
            name = f"activation_{target}_{tag}"
            cmd = format_command(
                [
                    *(f"{k}={v}" for k, v in values.items()),
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # Legacy sweep (+/-20%)
    for param, default in LEGACY_LOSS_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            name = f"legacy_{param}_{tag}"
            cmd = format_command(
                [
                    f"{param}={value}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # Triage sweep (specific values)
    for val in TRIAGE_VALUES:
        name = f"legacy_triage_{val:.3f}"
        cmd = format_command(
            [
                "loss=mse",
                "losslegacyweighting=t",
                f"ptriage={val}",
                f"ntriage={val}",
            ],
            name,
        )
        lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # Legacy offset sweep (zero baselines with +/- absolute adjustments)
    for param, step in LEGACY_OFFSET_PARAMS.items():
        for direction, tag in ((1, "high"), (-1, "low")):
            value = direction * step
            name = f"legacy_{param}_{tag}"
            cmd = format_command(
                [
                    f"{param}={value}",
                    "loss=mse",
                    "losslegacyweighting=t",
                    "ptriage=0.01",
                    "ntriage=0.01",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # WBCE sweep
    for param, default in WBCE_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            name = f"wbce_{param}_{tag}"
            cmd = format_command(
                [
                    "loss=wbce",
                    f"{param}={value}",
                    "losslegacyweighting=f",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    # Tversky sweep (alpha/beta only)
    for param, default in TVERSKY_PARAMS.items():
        for delta, tag in ((1.2, "high"), (0.8, "low")):
            value = default * delta
            name = f"tversky_{param}_{tag}"
            cmd = format_command(
                [
                    "loss=tversky",
                    f"{param}={value}",
                    "losslegacyweighting=f",
                ],
                name,
            )
            lines.append(f"{sbatch_prefix} --wrap={shlex.quote(cmd)}\n")

    with open(OUTPUT_SCRIPT, "w") as out:
        out.write("#!/bin/bash\n")
        out.write(f"mkdir -p {LOG_DIR}\n")
        out.write(f"mkdir -p {NET_DIR}\n")
        for line in lines:
            out.write(line)


if __name__ == "__main__":
    main()
