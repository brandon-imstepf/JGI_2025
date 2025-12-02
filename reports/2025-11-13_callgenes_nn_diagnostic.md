# CallGenes NN Diagnostic (2025-11-13)

## Goal
Investigate why the current NN multiplier hurts accuracy despite the plumbing being correct, by checking (1) label/oracle quality, (2) feature distribution drift, (3) training target vs runtime usage, and (4) generalization gaps relative to the `ncbi_master_11-7_splits` training set.

## Data & Commands
1. **Oracle sweep** – demonstrates integration health and candidate coverage limits:
   ```bash
   ./compile.sh
   (cd bbmap/current && ./run_nn_oracle_validation.sh)
   ```
   Results: `nn_oracle_validation/summary_oracle.csv`.

2. **Training vs validation vectors** – compared first eight scalar features (log length, start/mid/stop scores, gene/contig GC & entropy) between the training set that produced `cycles_trained_8dims_400k_ncbi_total_8f_6l_small_11-12_new.bbnet` and a freshly regenerated Validation_NCBI_10 TSV:
   ```bash
   python3 scripts (ad-hoc) sampling
   - train: ncbi_master_11-7_splits/ncbi_master_11-7_train.tsv
   - val:   cortex/.../btsf_runs/Validation_NCBI_10_training.tsv
   ```

3. **Label provenance** – inspected `slurm/run_builder*.slurm` and `bbmap/btsf_new.sh` to confirm both training and validation TSVs come from the same pipeline, and reviewed training logs (`slurm/cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-12_new.err`).

## Findings
### 1. Oracle Modes Confirm the Plumbing but Expose Candidate Gaps
| Oracle Mode | TP_Stop | FP_Stop | FN_Stop | Stop Recall | Stop Precision | Source |
|-------------|-------:|-------:|-------:|------------:|---------------:|--------|
| nn          | 40 420 | 1 215 | 3 219 | 0.9262 | 0.9708 | `summary_oracle.csv` |
| neutral     | 40 507 | 1 086 | 3 132 | 0.9282 | 0.9739 | `summary_oracle.csv` |
| truth       | 38 034 |   751 | 5 605 | 0.8716 | 0.9806 | `summary_oracle.csv` |

- **Neutral** matches the baseline, proving the multiplier is a no-op when requested; there is no hidden bias in the integration layer.
- **Truth** drastically reduces FP_Stop (−38 %) yet recall *drops* to 0.87 because many reference genes never appear in the candidate list. CallGenes can only rescale existing ORFs; if the heuristics miss or trim an ORF, even a perfect oracle cannot recover it. This is clearly visible on `GCF_030000695.1`, where truth mode still reports 902 FN_stop. Any future NN gains therefore hinge on improving candidate generation (e.g., second-pass heuristics, relaxed minlen) in addition to the multiplier.

### 2. Feature Drift Between Training and Validation
Sampling 100k rows from each source (training rows shuffled, validation freshly regenerated):

| Feature (first 8 dims) | Train Mean ± SD | Validation Mean ± SD |
|-----------------------|-----------------|----------------------|
| log_len/10            | 0.663 ± 0.067   | 0.599 ± 0.096        |
| start_score           | 0.600 ± 0.310   | 0.291 ± 0.334        |
| kmer_score/len        | 0.504 ± 0.246   | 0.317 ± 0.296        |
| stop_score            | 0.154 ± 0.191   | 0.120 ± 0.166        |
| gene_gc               | 0.532 ± 0.133   | 0.559 ± 0.112        |
| gene_entropy (logit)  | 0.779 ± 0.039   | 0.789 ± 0.053        |
| contig_gc             | 0.530 ± 0.129   | 0.558 ± 0.104        |
| contig_entropy (logit)| 0.775 ± 0.023   | 0.782 ± 0.026        |

(Validation stats from `btsf_runs/Validation_NCBI_10_training.tsv`; training stats from `ncbi_master_11-7_splits/ncbi_master_11-7_train.tsv`.)

Takeaways:
- Validation ORFs are, on average, shorter with weaker start/k-mer scores (likely because Validation_NCBI_10 contains tougher genomes); the NN never saw this distribution during training, so its outputs are poorly calibrated.
- GC and entropy means shift by >0.02, enough to move many inputs closer to the NN cutoff boundary.
- The training TSV appears to list positives first (our first 500k-line sample had label 1 throughout), so if shuffling was omitted before batching, every worker would see long runs of positives followed by negatives, reinforcing bias.

### 3. Training Target vs Runtime Usage
- `ml.Trainer` reports FPR≈FNR≈0.024 on training batches (see `slurm/cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-12_new.err`), meaning the model learned to separate the dataset it saw.
- At runtime we treat the sigmoid output as a *multiplier*, not a binary decision. When the output distribution is highly bimodal (close to 0 or 1), minor calibration errors cause large swings once divided by `nncutoff`. Without post-training calibration (e.g., temperature scaling or Platt scaling), shifting to a multiplier magnifies any bias.
- The cutoff (0.6) and strength (1.0) were tuned on older datasets (`ncbi_master_11-5`). Applying them wholesale to `ncbi_master_11-7` + Validation_NCBI_10 likely mis-scales the multiplier.

### 4. Label Alignment & Coverage Issues
- `btsf_new.sh` labels positives via exact start/stop equality between `callgenes.sh nofilter` output and the reference GFF. Any partial overlaps, alternative start codons, or CDS trimmed by the heuristic go to the negative pile. When we consume Validation_NCBI_10, the same strict rule remains, so borderline ORFs (the ones the NN should help with) are often mislabeled during both training and evaluation.
- Truth-mode recall < 1 confirms there are true genes absent from the candidate set even with `nofilter`, so the model never sees them during training either. Improving base heuristics (e.g., re-enabling 2-pass mode, relaxing `minlen`, or adding alternate start codon sweeps) would expand the candidate pool and give the NN a chance to learn real differences instead of memorizing high-score heuristics.
- Because positives and negatives are subsampled per genome (negatives downsampled to 2× positives), any genome whose negatives survive the shuf more frequently will dominate the batches. With the new Validation_NCBI_10 genomes absent from `ncbi_master_11-7`, the model is extrapolating beyond its training distribution.

## Recommendations
1. **Rebuild the training TSV** with the latest CallGenesHelper (prefilter disabled, oracle support) so positives/negatives match the runtime behavior. (`run_builder.slurm` already points at `NCBI_Datasets_Training`; rerun only after reviewing gating changes.)
2. **Shuffle the TSV before training** (either via `shuf` or inside `ml.Trainer`) to avoid long runs of positives that bias the optimizer.
3. **Calibrate the NN outputs** against the new validation set—e.g., fit a logistic regression or temperature scaling on Validation_NCBI_10 to map NN scores to meaningful multipliers before touching `nnStrength`.
4. **Expand candidate coverage**: re-enable 2-pass mode for neural runs, lower `minlen`, or add an option that injects reference starts/stops when available so the “perfect oracle” can actually reach 100% recall.
5. **Consider tolerance windows** when labeling positives (±3 nt on either end) to avoid mislabeling near-matches as negatives. This change would also make the truth-oracle evaluation more representative.

Once those changes are in place, regenerate `ncbi_master_11-7` via `run_builder.slurm` (after approval) and retrain with `train_total_cycles_8dim_6l_small_11-12.slurm`. Until then, tweaking cutoffs alone will continue to underperform because the model is calibrated to a different feature/label distribution.
