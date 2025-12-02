# CallGenes Centered Multiplier Experiment (2025-11-19)

## Summary
To test whether giving the neural network full authority over the heuristic score would suppress false positives, I added a `centered` multiplier mode inside `CallGenesHelper`. This mode maps the advisory probability directly onto a signed multiplier: `mult = (advisory − cutoff) / cutoff`, then scales the heuristic score by `mult * nnStrength`. Advisories below the cutoff yield negative multipliers, so their scores flip sign. After recompiling, I reran the entire Validation_NCBI_10 benchmark suite (`run_callgenes_modes.sh`) to compare the `centered` mode against the legacy ratio behavior.

## Implementation Notes
- `nnmultiplier_mode=centered` is parsed via `nnmultiplier_mode=<value>` and produces `mult = (advisory − cutoff) / cutoff` (with a safety guard when the cutoff is zero). The final score becomes `originalScore * (mult * nnStrength)`.
- Every other mode continues to funnel through the existing `((desired − 1) * nnStrength) + 1` blend, so only the `centered` mode produces negative scores when the NN is pessimistic.
- `run_callgenes_modes.sh` was invoked with `NN_MULT_MODE=centered`, generating an isolated results folder `callgenes_mode_results_centered/` so metrics can be compared alongside the earlier runs.

## Validation_NCBI_10 Results (two-pass runs)
Aggregating the ten genomes:

| Mode           | TP_Stop | FP_Stop | FN_Stop | TP_Start | FP_Start | FN_Start |
|----------------|--------:|--------:|--------:|---------:|---------:|---------:|
| ratio          | 40 363  | 1 293   | 3 276   | 34 831   | 6 825    | 8 808    |
| centered       | 40 363  | 1 293   | 3 276   | 34 831   | 6 825    | 8 808    |

The centered multiplier matched the legacy ratio totals exactly; neither precision nor recall changed. Per-genome numbers confirm the same pattern: false positives still carry enough residual score to survive the DP trace even when the NN assigns a negative multiplier.

## Conclusion & Next Steps
The `centered` mapping proved that merely flipping the score sign isn’t enough—the DP trace still considers those ORFs unless we explicitly reject negative scores or recalibrate the NN to produce stronger signal. To make real progress we likely need:
1. A hard rejection threshold (`advisory < t ⇒ orf.orfScore = -1`), or
2. A retrained/calibrated network whose outputs span the full `[0,1]` range on real CallGenes candidates.

The new mode and benchmark results remain in the repo for future comparison:
- Code changes: `bbmap_fresh/current/prok/CallGenesHelper.java` (centered mode + logging)
- Benchmark outputs: `callgenes_mode_results_centered/summary_results_callgenes.csv`

Further experiments (thresholding, recalibration, retraining) can now reuse the same benchmarking script and notebook visualizations.
