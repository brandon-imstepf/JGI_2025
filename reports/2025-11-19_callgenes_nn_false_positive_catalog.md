# CallGenes NN False-Positive Catalog & Score Audit (2025-11-19)

## Goal
- Reproduce the current CallGenes vs CallGenes+NN performance on `Validation_NCBI_10` with the NN prefilter disabled to confirm the regression persists.
- Build a reusable catalog of neural false positives along with truth-oracle traces so we can quickly retest future fixes.
- Inspect the neural score distribution (original vs modified scores, advisory multipliers) to understand why high-scoring false positives survive the multiplier.

## Commands & Artifacts
```bash
# Compile BBTools/CallGenes (run from repo root)
./compile.sh

# CallGenes-only benchmark sweep (Runs 1/2 pass ± NN w/ prefilter_cutoff=0)
(cd bbmap_fresh/current && ./run_callgenes_modes.sh)
# -> callgenes_mode_results/summary_results_callgenes.csv

# Generate GFFs + FP tables for several genomes (NN, 2-pass)
(cd bbmap_fresh/current && ./callgenes.sh …)  # invoked via helper loop
python3 python/list_false_positives.py <pred.gff> <truth.gff.gz> --limit 25 \
    > oracle_fp_catalog/<genome>.nn_fp_cds.tsv  # (and --feature rRNA)

# Truth-mode oracle traces for top-5 CDS false positives per genome
(cd bbmap_fresh/current && ./callgenes.sh … oracle=truth \
    oracle_debug=<contig:start-stop:strand> oracle_debug_log=<logs/...>)

# Score CSV for detailed multiplier stats (requires log=t truegenes=…)
(cd bbmap_fresh/current && ./callgenes.sh … log=t truegenes=…)
```
Artifacts:
- `callgenes_mode_results/summary_results_callgenes.csv` – latest benchmark table.
- `oracle_fp_catalog/` – per-genome neural GFFs, FP tables (CDS & rRNA), score logs.
- `logs/2025-11-19_oracle_trace_*` – 20 oracle traces (5 CDS FPs × 4 genomes) confirming each ORF has oracle score 0 despite high heuristic scores.
- `python/list_false_positives.py` – new helper to dump the first N FPs with ready-to-use `oracle_debug` values.

## Findings
### 1. NN still regresses recall & precision (prefilter disabled)
Summing the ten `Validation_NCBI_10` genomes:

| Mode | TP_Start | TP_Stop | FP_Start | FP_Stop | FN_Start | FN_Stop |
|------|---------:|--------:|---------:|--------:|---------:|--------:|
| single_pass | 34 827 | 40 507 | 6 766 | 1 086 | 8 812 | 3 132 |
| single_pass_nn | 34 146 | 39 997 | 7 126 | 1 275 | 9 493 | 3 642 |
| two_pass | 35 762 | 41 030 | 6 364 | 1 096 | 7 877 | 2 609 |
| two_pass_nn | 34 831 | 40 363 | 6 825 | 1 293 | 8 808 | 3 276 |

Even without the NN-specific prefilter, the neural modes lose ~900 TP starts (−2.6 %) and add ~667 FN stops (+26 %) while also increasing FP stops (+197). This mirrors the November 13 report, so the integration change alone did not fix the regression.

### 2. Neural outputs cluster below the cutoff for both TPs and FPs
Using `log=t truegenes=…` we captured `scores.csv` for `GCF_901004845.1_7614_7_26_genomic` (2-pass NN run). The NN advisory score is inferred from the modified/original ratio (`advisory = (modified/original) × nnCutoff`, with `nnCutoff = 0.6`):

| Status | Count | Mean Advisory | P10 | Median | P90 | Mean Multiplier |
|--------|------:|--------------:|----:|-------:|----:|----------------:|
| TP | 3 303 | 0.57 | 0.52 | 0.57 | 0.61 | 0.95 |
| FP | 3 888 | 0.48 | 0.44 | 0.48 | 0.52 | 0.80 |

- Only 598 of 3 303 true positives (18 %) have advisory ≥ 0.6; 91 % of TPs are *below* the cutoff, so their scores are reduced.
- False positives do sit lower on average (0.48 vs 0.57) but the overlap is huge; the median multipliers differ by only ~0.15. A 15–20 % penalty is not enough to dislodge high-score FPs (e.g., the `NZ_CAAJWH010000040.1:1050-1550-` FP drops from 303 → 259 yet still wins the DP trace).
- Because `nnStrength` is capped at ≤ 1, we cannot amplify penalties beyond `advisory / cutoff`. The current NN must therefore output values near zero to actually zero out an ORF, but the runtime distribution never approaches that.

### 3. Confirmed FP catalog for regression testing
`python/list_false_positives.py` + oracle runs produced ready-to-replay coordinates for four genomes (top 25 CDS + rRNA). Examples (all logged under `logs/2025-11-19_*`):

| Genome | Feature | Example | Heuristic Score | Oracle |
|--------|---------|---------|----------------:|--------|
| GCF_015023715.1 | CDS | `NZ_JACVDX010000003.1:3-1601:+` | 408 | `oracle_truth = 0 → rejecting` |
| GCF_015023775.1 | CDS | `NZ_JACVDR010000003.1:818-2548:+` | 1 208 | `oracle_truth = 0 → rejecting` |
| GCF_031010175.1 | CDS | `NZ_JAGTGN010000001.1:96994-97440:+` | 222 | `oracle_truth = 0 → rejecting` |
| GCF_901004845.1 | CDS | `NZ_CAAJWH010000040.1:1050-1550:-` | 303 | `oracle_truth = 0 → rejecting` |

Each coordinate now has an accompanying oracle debug log that shows the truth score is 0 while the heuristic score stays high, matching the earlier Mruber case.

### 4. Training job aborted after RNG bug
`slurm/cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-12_new.err` ends with the pre–FastRandom fix assertion and a SLURM timeout (job 19698473). That means the `.bbnet` we are using comes from a run that never completed an evolutionary cycle after 400k batches. Combined with the observed score distribution drift, this strongly suggests we are using a stale/under-trained network that was never recalibrated after the RNG fix + label-tolerance changes.

## Next Steps
1. **Regenerate training data & retrain with the fixed RNG.** Rerun `run_builder_tol3.slurm` (prefilter disabled) followed by `train_total_cycles_8dim_6l_small_11-12.slurm` now that `FastRandom` is fixed, then benchmark the resulting `.bbnet`.
2. **Calibrate neural outputs before scaling.** Options include:
   - Relaxing the cutoff (e.g., calibrate via Validation_NCBI_10 so positives sit ≥0.6 while negatives fall below), or
   - Allowing `nnStrength > 1` / introducing a nonlinear penalty so a 0.48 advisory actually cuts the score in half instead of 20 %.
3. **Leverage the FP catalog.** Re-run `python/list_false_positives.py` and the oracle traces after each tweak to ensure those known-bad coordinates finally get suppressed.
4. **Automate score logging in benchmarks.** When `log=t truegenes=…` is set we get rich `scores.csv` files—hook that into the benchmarking scripts so we can monitor advisory distributions per genome.
