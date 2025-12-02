# Session Summary – Oracle Tracing & Setup (2025-11-18)

## Goal
- Enable reproducible tracing for CallGenes false positives (CDS and RNA) so we can diagnose why the neural-net multiplier fails to suppress them.

## Work Completed
1. **Helper script placement** – Moved the “find first false positive” utility into the top-level `python/` folder (INIT.md updated) so non-core scripts stay outside `bbmap_fresh/`.
2. **False-positive discovery** – Ran `python3 python/find_first_false_positive.py …` on `oracle_test_11-18.gff` to capture the first CDS and rRNA cases:
   - CDS FP: `NC_013946.1:1964-2947:+`
   - rRNA FP: `NC_013946.1:256481-257974:+`
3. **Oracle tracing runs** – Replayed Mruber with `oracle=truth` and `oracle_debug` for each FP, logging to:
   - `logs/2025-11-18_oracle_trace_cds.log`
   - `logs/2025-11-18_oracle_trace.log` (rRNA)
   Both traces confirm an oracle score of 0 despite high heuristic scores.
4. **Temporary instrumentation** – Briefly forced non-CDS targets through `CallGenesHelper.modifyOrfScoreWithNeuralNetwork` when they matched `oracle_debug`, captured the rRNA trace, then removed the hook so production code remains unchanged.
5. **Catalog snapshot** – Quick script (untracked) counted ~600 CDS false positives and 13 RNA false positives in the `oracle=nn` run; top entries recorded in the shell history for follow-up.

## Next Steps
1. Expand the FP catalog (multiple CDS + RNA coordinates) and rerun truth-mode traces for each, dropping the logs under `logs/`.
2. Extract feature vectors and raw NN outputs for those ORFs and neighboring true positives to see whether the NN saturates or heuristics overpower the multiplier.
3. Summarize the patterns and propose training/scaling/multiplier adjustments, then rerun `run_benchmarking_all.sh` and refresh the notebook.

Artifacts from this session (GFFs, logs) live under `/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/{oracle_*,logs/}` and can be used to resume the analysis.
