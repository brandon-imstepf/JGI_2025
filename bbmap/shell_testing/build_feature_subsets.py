#!/usr/bin/env python3
import argparse
import gzip
import os
from pathlib import Path
from typing import List, Tuple

FLOAT_COUNT = 8

SUBSET_NAMES = ("full", "floats_only", "onehot_only", "floats_plus_onehot_no_middle")


def infer_total_columns(sample_path: Path) -> int:
    with gzip.open(sample_path, "rt") as fh:
        first_line = fh.readline().rstrip("\n")
    return len(first_line.split("\t"))


def compute_masks(total_cols: int) -> dict:
    label_idx = total_cols - 1
    onehot_total = max(0, label_idx - FLOAT_COUNT)
    segment = onehot_total // 3
    remainder = onehot_total - segment * 3
    start_block = FLOAT_COUNT
    first_end = start_block + segment
    last_start = start_block + segment + remainder + segment
    masks = {}
    for name in SUBSET_NAMES:
        if name == "full":
            indices = list(range(total_cols))
        elif name == "floats_only":
            indices = list(range(FLOAT_COUNT))
            indices.append(label_idx)
        elif name == "onehot_only":
            indices = list(range(start_block, label_idx))
            indices.append(label_idx)
        elif name == "floats_plus_onehot_no_middle":
            indices = list(range(FLOAT_COUNT))
            indices.extend(range(start_block, first_end))
            indices.extend(range(last_start, label_idx))
            indices.append(label_idx)
        else:
            indices = list(range(total_cols))
        masks[name] = indices
    return masks


def sanitize_label(value: str) -> str:
    stripped = value.replace("\x00", "").strip()
    if not stripped:
        return stripped
    for ch in reversed(stripped):
        if ch in ("0", "1"):
            return ch
    return stripped


def stream_subsets(fin_path: Path, masks: dict, out_dir: Path, label_idx: int):
    out_dir.mkdir(parents=True, exist_ok=True)
    writers = {}
    try:
        for subset, mask in masks.items():
            out_path = out_dir / subset / fin_path.name.replace(".tsv.gz", f".subset_{subset}.tsv.gz")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            writers[subset] = (gzip.open(out_path, "wt"), mask, max(mask))
        with gzip.open(fin_path, "rt") as fin:
            for line in fin:
                clean = line.replace("\x00", "").rstrip("\r\n")
                parts = clean.split("\t")
                plen = len(parts)
                if plen <= label_idx:
                    continue
                label_value = sanitize_label(parts[label_idx])
                if label_value not in {"0", "1"}:
                    continue
                parts[label_idx] = label_value
                for subset, (writer, mask, max_idx) in writers.items():
                    if plen <= max_idx:
                        continue
                    subset_vals = [parts[i] for i in mask]
                    writer.write("\t".join(subset_vals))
                    writer.write("\n")
    finally:
        for writer, _, _ in writers.values():
            writer.close()


def build_subsets(input_dir: Path, output_dir: Path):
    input_files = {
        "all": next(input_dir.glob("*.all.tsv.gz")),
        "positives": next(input_dir.glob("*.positives.tsv.gz")),
        "negatives": next(input_dir.glob("*.negatives.tsv.gz")),
    }
    total_cols = infer_total_columns(input_files["all"])
    masks = compute_masks(total_cols)
    label_idx = total_cols - 1
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, src in input_files.items():
        stream_subsets(src, masks, output_dir, label_idx)


def main():
    parser = argparse.ArgumentParser(description="Build feature subsets for training datasets.")
    parser.add_argument("--input-dir", required=True, help="Path to original dataset folder.")
    parser.add_argument("--output-dir", required=True, help="Path for subset outputs.")
    args = parser.parse_args()
    build_subsets(Path(args.input_dir), Path(args.output_dir))


if __name__ == "__main__":
    main()
