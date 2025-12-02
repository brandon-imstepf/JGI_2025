# Neural-Net Prefilter Investigation (2025-11-13)

- **Context**: Benchmarking in `benchmarking_results_all/summary_results.csv` showed CallGenes+NN runs consistently had lower recall (higher FN) than baseline despite similar FP counts.
- **Cause**: When `net=` is supplied, `CallGenesHelper` (see `bbmap/current/prok/CallGenesHelper.java:823-835`) applies a hard prefilter: any ORF with `orfScore / length < prefilterCutoff` (default `0.1`) is forced to `orf.orfScore = -1` before the NN multiplier executes. Baseline runs do not perform this filter, so NN mode silently discards marginal ORFs instead of merely rescaling them.
- **Evidence**:
  * Genome `GCF_901004845.1_7614_7_26_genomic` (single-pass) baseline vs NN:
    - Baseline (`…baseline.gff`): `TP_Stop=3233`, `FN_Stop=459`.
    - NN with `prefilter_cutoff=0.1` (`…nn_pref01.gff`): `TP_Stop=3179`, `FN_Stop=513` (−54 TP, +54 FN).
    - NN with `prefilter_cutoff=0.0` (`…nn_pref00.gff`): `TP_Stop=3233`, `FN_Stop=459` (matches baseline while retaining NN scoring for reweighting).
  * Similar FN inflation is visible across all validation genomes (see the new CallGenes comparison section in `python/data_analysis.ipynb` §6).
- **Recommendation**:
  1. Disable the NN-specific prefilter during benchmarking by passing `prefilter_cutoff=0` (or `nnallorfs=true`) so the NN operates purely as the expert multiplier.
  2. If a prefilter is needed for speed, lower it substantially and add logging to quantify how many ORFs it rejects.
  3. Re-run `bbmap/current/run_benchmarking_all.sh` after adjusting the parameters, then refresh the analysis notebook to confirm the NN no longer regresses recall.

Files produced during this analysis:
- `benchmarking_results_all/GCF_901004845.1_7614_7_26_genomic.baseline.gff`
- `benchmarking_results_all/GCF_901004845.1_7614_7_26_genomic.nn_pref01.gff`
- `benchmarking_results_all/GCF_901004845.1_7614_7_26_genomic.nn_pref00.gff`

Use `java -ea -cp bbmap/current gff.CompareGff …` to replicate the metrics above.
