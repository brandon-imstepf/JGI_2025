# FastRandom RNG Fix & Recompile (2025-11-17)

## Goal
- Diagnose and fix the negative-layer assertion inside `ml.Trainer` runs that start from `cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig_tol3`, then ensure the tree builds with the corrected RNG.

## Commands
```bash
# Patch FastRandom then rebuild CallGenes/ML classes
./compile.sh
```

## Findings
- The crash originated in `CellNet.pickEdges` because `FastRandom.nextInt(int bound)` was casting `(nextLong() >>> 1)` to `int`, which can wrap around to a negative value when the upper long bit is set. Those negative draws fed `rand` and `start` inside `pickEdges`, tripping the assertion seen in the Slurm log.
- Replaced `nextInt(int bound)` with a thin wrapper around `nextLong(long bound)` so bounded draws remain non-negative regardless of the parent RNG state (`bbmap/current/shared/FastRandom.java:95-101`).
- Recompiled via `./compile.sh` to refresh `ml.Trainer` and dependent classes; build completed successfully (only the usual incubating-module warning remained).

## Next Steps
1. Relaunch `cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig_tol3` to confirm the RNG change removes the assertion and training progresses.
2. If the run succeeds, capture the new `.bbnet` path + metrics in a follow-up report so future agents can trace when the fix landed.
