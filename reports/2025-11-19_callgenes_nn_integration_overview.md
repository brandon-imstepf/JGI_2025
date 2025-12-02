# CallGenes NN Integration Overview (2025-11-19)

## Helper Architecture & Hooks
- **Multiplier Registry** – `bbmap/current/prok/CallGenesHelper.java:60-139` defines the configurable knobs: `nnCutoff`, `nnStrength`, and the `nnmultiplier_mode` enum (`ratio`, `ratio_pow2`, `ratio_pow3`, `advisory`, `sigmoid`, `exponential`, `direct_only`, `centered`). This block also includes the temporary verbose trace for `NZ_CAAJWH010000040.1:1050-1550-`, which writes detailed logs to `/clusterfs/.../logs/2025-11-19_nn_fp_verbose.log`.
- **CLI Parsing** – `bbmap/current/prok/CallGenesHelper.java:170-205` augments the parser with `nnmultiplier_mode`, `nn_sigmoid_slope`, and `nn_exp_alpha`, so any `callgenes.sh` invocation (or `run_callgenes_modes.sh` via `NN_MULT_MODE`) can toggle the strategy without code edits.

## Advisory Score Flow
- **Feature Generation & NN Call** – Inside `modifyOrfScoreWithNeuralNetwork` (`bbmap/current/prok/CallGenesHelper.java:937-1053`) each ORF’s 356-dim vector is fed to a thread-local `CellNet` (`applyInput` + `feedForward`). The raw advisory score is logged when the verbose hook is armed, and written to `scores.csv` along with the ORF ID, original score, modified score, and TP/FP label.
- **Multiplier Application** – After the advisory is computed, the helper chooses a scaling path:
  - `DIRECT_ONLY`: `modifiedScore = advisoryScore` (heuristic ignored).
  - `CENTERED`: `mult = (advisory − cutoff)/cutoff`, then `modifiedScore = originalScore * (mult * nnStrength)` (can flip scores negative).
  - All other modes: `desired = computeDesiredMultiplier(advisory)`, `finalMult = ((desired − 1) * nnStrength) + 1`, then `modifiedScore = originalScore * finalMult`.
  `computeDesiredMultiplier` (lines 1112-1155) implements the ratio/power/advisory/sigmoid/exponential transforms, so experimenting with penalty curves just means extending that function.
- **Truth-Oriented Guards** – When `oracle=truth`, the block just above (lines 1013-1024) adds zero-scored ORFs to `threadTruthRejects`, forcing `applyTruthRejections` (lines 1072-1095) to drop them before the DP trace.

## CallGenes Integration Points
- **Helper Initialization** – `bbmap/current/prok/CallGenes.java:30-124` instantiates `CallGenesHelper`, calls `helper.initialize`, and lets it create the log directory/scores path when `log=t truegenes=...` is set. That’s why every `run_callgenes_modes.sh` run now emits per-mode `log.txt` and `scores.csv` files.
- **Thread Usage** – `CallGenes.java:779-818` shows each worker thread calling `CallGenes.this.helper.getThreadLocalCopy()` and wiring the score writer via `helper.setScoreWriter(bswScore)`. From there every ORF goes through `CallGenesHelper.modifyOrfScoreWithNeuralNetwork` before dynamic programming chooses the final gene set.
- **DP Trace Context** – The per-list loop (`CallGenes.java:828-838`) funnels into the helper, so any change in `CallGenesHelper` (reject thresholds, multiplier math) directly affects which ORFs survive into the DP path.

## Diagnostics & Visualization
- **Verbose Log** – `/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/2025-11-19_nn_fp_verbose.log` captures the exact advisory/multiplier math for the hard-coded FP (`NZ_CAAJWH010000040.1:1050-1550-`). This is useful when calibrating thresholds.
- **Notebook Plots** – `python/data_analysis.ipynb` has three sections:
  1. Benchmarking tool comparison (`benchmarking_results_all/summary_results.csv`).
  2. Multiplier sweep aggregation (`callgenes_mode_results_*`), plotting stop precision/recall per mode/pass.
  3. Score multiplier histograms that ingest every `scores.csv` and show average/median `Modified/Original` ratios for TP vs FP.

With these references you should have clickable context for both the helper and the integration, plus plots/logs to evaluate any future tweaks.
