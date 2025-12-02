#!/usr/bin/env python3
"""
Find the first false-positive CDS in a CallGenes output GFF by comparing
against a reference GFF (default: CDS entries only).  Prints both the
contig/start/stop/strand tuple and a ready-to-copy oracle_debug argument.
"""
import argparse
import gzip
import os
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_truth(path, feature="CDS"):
    truth = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            if feature and parts[2] != feature:
                continue
            contig = parts[0].split()[0]
            start = int(parts[3])
            stop = int(parts[4])
            strand = "+" if parts[6] == "+" else "-"
            truth.add((contig, start, stop, strand))
    return truth


def find_first_fp(pred_path, truth_set, feature="CDS"):
    with open_text(pred_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            if feature and parts[2] != feature:
                continue
            contig = parts[0].split()[0]
            start = int(parts[3])
            stop = int(parts[4])
            strand = "+" if parts[6] == "+" else "-"
            quad = (contig, start, stop, strand)
            if quad not in truth_set:
                return quad
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Identify the first false-positive feature in a CallGenes GFF."
    )
    parser.add_argument("pred_gff", help="Predicted GFF (CallGenes output)")
    parser.add_argument("truth_gff", help="Reference GFF (truth set)")
    parser.add_argument(
        "--feature",
        default="CDS",
        help="Feature type to compare (default: CDS)",
    )
    args = parser.parse_args()

    truth = parse_truth(args.truth_gff, args.feature)
    if not truth:
        sys.exit("No truth entries found; check --feature or the truth GFF.")

    fp = find_first_fp(args.pred_gff, truth, args.feature)
    if not fp:
        print("No false positives detected for feature", args.feature)
        return

    contig, start, stop, strand = fp
    print(f"First false positive: {contig}:{start}-{stop}{strand}")
    print(
        "oracle_debug argument:",
        f'oracle_debug={contig}:{start}-{stop}:{strand}',
    )


if __name__ == "__main__":
    main()
