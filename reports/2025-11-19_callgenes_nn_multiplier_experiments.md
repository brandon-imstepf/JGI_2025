# CallGenes NN Multiplier Experiments & Logging (2025-11-19)

## Goals
- Capture full traces for known CallGenes false positives so we can see how the NN modifies heuristics in real time.
- Explore alternative multiplier strategies (ratio powers, direct advisory, sigmoid, exponential, NN-only) to see if any improve Validation_NCBI_10 benchmarks without hard gates.
- Ensure `run_callgenes_modes.sh` always logs scores so diagnostics are reproducible, and surface the results in `python/data_analysis.ipynb` for quick visualization.

## Code / Tool Changes
1. **Verbose FP instrumentation** – `bbmap_fresh/current/prok/CallGenesHelper.java` now contains a temporary hard-coded hook for `NZ_CAAJWH010000040.1:1050-1550-`. Every time that ORF is processed the helper logs entry conditions, advisory score, multiplier mode, and final score to `/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/2025-11-19_nn_fp_verbose.log`. This showed the NN outputs (~0.51) only shave ~15% off a 300+ heuristic score, explaining why the ORF survives.
2. **Logging defaults** – `bbmap_fresh/current/run_callgenes_modes.sh` now always passes `truegenes=... log=t`, so every run leaves a `log.txt` and `scores.csv` per genome/mode under `callgenes_mode_results_<mode>/<genome>.<mode>/`. Each multiplier sweep therefore captures its own score histograms.
3. **Multiplier framework** – Added `nnmultiplier_mode` (ratio, ratio_pow2, ratio_pow3, advisory, sigmoid, exponential, direct_only) with optional `nn_sigmoid_slope` and `nn_exp_alpha`. The new modes let us map the NN advisory into aggressive penalties/boosts or even replace the heuristic score entirely (`direct_only`). `run_callgenes_modes.sh` accepts `NN_MULT_MODE` and writes to `callgenes_mode_results_<mode>/` so runs don’t overwrite each other.
4. **Notebook recreation** – Rebuilt `python/data_analysis.ipynb` after it was accidentally zeroed. The new version has two sections: (a) aggregate plots for `run_benchmarking_all.sh` outputs, and (b) a CallGenes multiplier comparison that reads every `callgenes_mode_results_*` folder and plots stop precision/recall bars.

## Benchmarks
Ran `bbmap_fresh/current/run_callgenes_modes.sh` on Validation_NCBI_10 with `NN_MULT_MODE` set to each strategy (ratio, ratio_pow2, ratio_pow3, advisory, sigmoid, exponential, direct_only). Key two-pass NN totals:

| Mode          | TP_Stop | FP_Stop | FN_Stop | TP_Start | FP_Start | FN_Start |
|---------------|--------:|--------:|--------:|---------:|---------:|---------:|
| ratio         | 40 363  | 1 293   | 3 276   | 34 831   | 6 825    | 8 808    |
| ratio_pow2    | 39 607  | 1 550   | 4 032   | 33 895   | 7 262    | 9 744    |
| ratio_pow3    | 39 136  | 1 575   | 4 503   | 33 312   | 7 399    | 10 327   |
| advisory      | 40 383  | 1 964   | 3 256   | 34 368   | 7 979    | 9 271    |
| sigmoid       | 39 593  | 2 233   | 4 046   | 33 455   | 8 371    | 10 184   |
| exponential   | 40 363  | 1 293   | 3 276   | 34 831   | 6 825    | 8 808    |
| direct_only   | 40 363  | 1 293   | 3 276   | 34 831   | 6 825    | 8 808    |

None of the new scaling modes improved on the legacy ratio behavior; most increase FN counts or FP counts. Direct-only and exponential modes produced identical metrics to ratio, confirming that simply scaling/replacing the heuristic doesn’t eliminate known FPs when the NN output stays near 0.5.

## Remaining Gaps / Next Ideas
- The NN outputs for Validation_NCBI_10 are too tightly clustered around the cutoff. To make the NN genuinely decisive we either need a hard reject threshold or a retrained/calibrated network (post-RNG fix and label tolerance) with stronger separation between positives and negatives.
- All multiplier modes now have score logs, so any future tweak can be evaluated quickly by re-running the sweep and inspecting `callgenes_mode_results_<mode>/summary_results_callgenes.csv` alongside the notebook plots.

All artifacts (logs, score CSVs, FP catalogs) live under `/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/` for future agents.
