#!/usr/bin/env python3
"""
Chunk-wise shuffler for large TSV files.

The script randomizes lines by shuffling fixed-size chunks in memory,
writing each shuffled chunk to a temporary file, shuffling chunk order,
then streaming the shuffled chunks to the final output. This avoids
loading the entire dataset into RAM while still providing uniform
randomization at the chunk granularity.
"""

import argparse
import os
import random
import shutil
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description="Shuffle large TSV files without loading everything into RAM.")
    parser.add_argument("--input", "-i", required=True, help="Path to the input TSV file.")
    parser.add_argument("--output", "-o", required=True, help="Path to the shuffled output TSV file.")
    parser.add_argument("--chunk", "-c", type=int, default=500_000,
                        help="Number of lines to shuffle per chunk (default: 500000).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--tmpdir", default=None, help="Optional temporary directory for chunk files.")
    return parser.parse_args()


def write_chunk(lines, tmpdir):
    random.shuffle(lines)
    tmp = tempfile.NamedTemporaryFile(delete=False, dir=tmpdir, mode="w", encoding="utf-8")
    tmp.writelines(lines)
    tmp.flush()
    tmp.close()
    return tmp.name


def shuffle_file(in_path, out_path, chunk_size, seed=None, tmpdir=None):
    if seed is not None:
        random.seed(seed)

    chunk_files = []
    lines = []

    with open(in_path, "r", encoding="utf-8") as infile:
        for line in infile:
            lines.append(line)
            if len(lines) >= chunk_size:
                chunk_files.append(write_chunk(lines, tmpdir))
                lines = []

    if lines:
        chunk_files.append(write_chunk(lines, tmpdir))

    random.shuffle(chunk_files)

    with open(out_path, "w", encoding="utf-8") as outfile:
        for chunk in chunk_files:
            with open(chunk, "r", encoding="utf-8") as cf:
                shutil.copyfileobj(cf, outfile)
            os.remove(chunk)


def main():
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"Input file not found: {args.input}")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    shuffle_file(args.input, args.output, args.chunk, args.seed, args.tmpdir)


if __name__ == "__main__":
    main()
