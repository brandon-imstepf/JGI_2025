#!/usr/bin/env python3
"""
List the first N false-positive features from a CallGenes GFF by comparing to
an authoritative truth set.  Helps build an oracle-debug catalog quickly.
"""

import argparse
import gzip
from typing import Iterable, List, Optional, Sequence, Tuple, Set


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_truth(path: str, feature: str) -> Set[Tuple[str, int, int, str]]:
    truth: Set[Tuple[str, int, int, str]] = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
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


def iterate_predictions(
    path: str, feature: str
) -> Iterable[Tuple[str, int, int, str, float]]:
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
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
            score = float(parts[5]) if parts[5] not in (".", "") else float("nan")
            yield (contig, start, stop, strand, score)


def collect_false_positives(
    preds: Iterable[Tuple[str, int, int, str, float]],
    truth: Set[Tuple[str, int, int, str]],
    limit: int,
) -> List[Tuple[str, int, int, str, float]]:
    results: List[Tuple[str, int, int, str, float]] = []
    for quad in preds:
        contig, start, stop, strand, score = quad
        if (contig, start, stop, strand) in truth:
            continue
        results.append(quad)
        if 0 < limit == len(results):
            break
    return results


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="List the first N false-positive features in a CallGenes output GFF."
    )
    parser.add_argument("pred_gff", help="Predicted GFF (CallGenes output)")
    parser.add_argument("truth_gff", help="Reference GFF (truth set)")
    parser.add_argument(
        "--feature",
        default="CDS",
        help="Feature type to compare (default: CDS)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Number of false positives to list (<=0 means no limit)",
    )
    args = parser.parse_args(argv)

    truth = load_truth(args.truth_gff, args.feature)
    if not truth:
        raise SystemExit("No truth entries found; check --feature or the truth GFF.")

    preds = iterate_predictions(args.pred_gff, args.feature)
    fps = collect_false_positives(preds, truth, args.limit)

    if not fps:
        print("No false positives detected for feature", args.feature)
        return 0

    print("Idx\tContig\tStart\tStop\tStrand\tScore\toracle_debug")
    for idx, (contig, start, stop, strand, score) in enumerate(fps, start=1):
        debug = f"{contig}:{start}-{stop}:{strand}"
        pretty_score = f"{score:.4f}" if score == score else "nan"
        print(f"{idx}\t{contig}\t{start}\t{stop}\t{strand}\t{pretty_score}\t{debug}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
