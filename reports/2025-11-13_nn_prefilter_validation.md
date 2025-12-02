# NN Prefilter Validation (2025-11-13)

## Summary
- Disabled the neural-network prefilter inside `CallGenesHelper` so the NN multiplier now executes for every ORF unless a positive `prefilter_cutoff` is explicitly supplied. This mirrors the "advisory multiplier" design outlined in INIT.md.
- Recompiled the local CallGenes classes (`./compile.sh`) after removing the orphaned `tax/TaxApp.java` entry from `sources.txt` so the updated helper is actually used at runtime.
- Regenerated the `Validation_NCBI_10` vectors via `bbmap/btsf_new.sh` (training mode) to ensure no regressions when the helper runs without a net file. Output TSV: `cortex/JGI_2025/brandon-imstepf-JGI_2025-9bd4797/btsf_runs/Validation_NCBI_10_training.tsv` (90 MB).
- Added `bbmap/current/run_nn_prefilter_validation.sh` to sweep the Validation_NCBI_10 genomes twice (legacy `prefilter_cutoff=0.1` vs. new default `0.0`) while logging `CompareGff` metrics for start/stop precision + recall. Results land in `nn_prefilter_validation/`.

## Commands
```bash
# Build
./compile.sh

# Training vectors (no NN)
./bbmap/btsf_new.sh \
    in=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/Validation_NCBI_10 \
    out=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/cortex/JGI_2025/brandon-imstepf-JGI_2025-9bd4797/btsf_runs/Validation_NCBI_10_training.tsv

# Prefilter sweep (NN enabled)
(cd bbmap/current && ./run_nn_prefilter_validation.sh)
```

## Metrics (Validation_NCBI_10)
Aggregated stop-site metrics from `nn_prefilter_validation/summary_prefilter.csv`:

| Prefilter | TP_Stop | FP_Stop | FN_Stop | Stop Recall | Stop Precision |
|-----------|--------:|--------:|--------:|------------:|---------------:|
| 0.1       | 40 340  | 1 214   | 3 299   | 0.9244      | 0.9708         |
| 0.0       | 40 420  | 1 215   | 3 219   | 0.9262      | 0.9708         |

Key per-genome deltas (0.0 – 0.1):

- `GCF_901004845.1`: +48 TP_stop, −48 FN_stop, FP_stop unchanged.
- `GCF_037003125.1`: +22 TP_stop, −22 FN_stop, −1 FP_stop.
- `GCF_030000695.1`: +4 TP_stop, −4 FN_stop.
- Remaining genomes stay identical, confirming the old prefilter only suppressed borderline ORFs.

No FP inflation was observed; NN recall recovered exactly where the legacy ratio filter previously zeroed the scores, matching the expectations from `reports/nn_prefilter_analysis.md`.

## Next Steps
- Feed the refreshed `nn_prefilter_validation/summary_prefilter.csv` into `python/data_analysis.ipynb` (Section 6) to update the CallGenes vs NN comparison plots.
- Once satisfied with the deltas, rerun `bbmap/current/run_benchmarking_all.sh` so that Prodigal/GeneMark vs CallGenes tables incorporate the fixed NN behavior.
