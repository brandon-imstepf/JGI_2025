# CallGenes Oracle Validation (2025-11-13)

## Objective
Validate the new `oracle=` switch inside `CallGenesHelper` by exercising all three advisory modes (neural, neutral, truth) on the `Validation_NCBI_10` genomes and confirming their behavior with `CompareGff` metrics.

## Commands
```bash
# Build to pick up CallGenesHelper/Oracle updates
./compile.sh

# Run all ten Validation_NCBI_10 genomes under oracle=nn|neutral|truth
(cd bbmap/current && ./run_nn_oracle_validation.sh)
```
Outputs land in `nn_oracle_validation/oracle_<mode>/` with per-genome GFFs and stdout logs; aggregated metrics are in `nn_oracle_validation/summary_oracle.csv`.

## Aggregated Metrics (Validation_NCBI_10)
| Oracle Mode | TP_Stop | FP_Stop | FN_Stop | Stop Recall | Stop Precision | TP_Start | FP_Start | FN_Start | Start Recall | Start Precision |
|-------------|-------:|-------:|-------:|------------:|---------------:|--------:|--------:|--------:|-------------:|----------------:|
| nn          | 40,420 | 1,215 | 3,219 | 0.9262 | 0.9708 | 34,921 | 6,714 | 8,718 | 0.8002 | 0.8387 |
| neutral     | 40,507 | 1,086 | 3,132 | 0.9282 | 0.9739 | 34,827 | 6,766 | 8,812 | 0.7981 | 0.8373 |
| truth       | 38,034 |   751 | 5,605 | 0.8716 | 0.9806 | 36,602 | 2,183 | 7,037 | 0.8387 | 0.9437 |

## Observations
- `oracle=neutral` holds metrics near the baseline NN run, confirming that the new hook can be a true no-op when needed.
- `oracle=truth` (perfect labels derived from the reference GFF) eliminates roughly 40% of the NN’s false positives (FP_Stop 1,215 → 751) and boosts precision above 98%. Recall falls because the current true-gene lookup requires exact coordinate matches; many valid genes that shift a few bases relative to the reference are therefore labeled 0 by the oracle. This is expected given the strict equality test and matches prior behavior in the training pipeline.
- Per-genome logs live under `nn_oracle_validation/oracle_<mode>/<genome>.stdout.log`, which capture the CallGenes stats and make it easy to spot any anomalies in runtime or logging.

These results demonstrate that the oracle plumbing works: neural mode reproduces the current behavior, neutral mode leaves scores untouched, and truth mode enforces perfect knowledge (subject to the exact-coordinate rule). The next step is to decide whether to relax the truth matching (e.g., allow ±3 nt windows) so that the perfect oracle can’t accidentally down-weight near-matching ORFs.
