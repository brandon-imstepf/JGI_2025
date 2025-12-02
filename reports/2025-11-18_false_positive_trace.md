# False Positive Trace – Mruber Genome (2025-11-18)

## Goal
- Identify a concrete false-positive ORF from an `oracle=nn` run, drive it through the truth-mode instrumentation, and capture the full decision trail so we can diagnose why the NN pipeline keeps it.

## Commands
```bash
# 1. Baseline NN run to harvest the predicted GFF
./bbmap_fresh/callgenes.sh \
    in=sh_test/Mruber.fna.gz \
    outgff=oracle_test_11-18.gff \
    truegenes=sh_test/Mruber.gff.gz \
    nn=slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet \
    nncutoff=0.6 strength=1.0 \
    oracle=nn log=t

# 2. Ask the helper script for the first rRNA false positive
python3 python/find_first_false_positive.py \
    oracle_test_11-18.gff \
    sh_test/Mruber.gff.gz \
    --feature rRNA

# 3. Truth-mode replay with oracle_debug logging enabled
./bbmap_fresh/callgenes.sh \
    in=sh_test/Mruber.fna.gz \
    outgff=oracle_debug_trace.gff \
    truegenes=sh_test/Mruber.gff.gz \
    nn=slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet \
    nncutoff=0.6 strength=1.0 \
    oracle=truth \
    oracle_debug=NC_013946.1:256481-257974:+ \
    oracle_debug_log=logs/2025-11-18_oracle_trace.log \
    log=t

# 4. Inspect the dedicated trace log
tail -n 40 logs/2025-11-18_oracle_trace.log
```

## Findings
- `python/find_first_false_positive.py` flagged `NC_013946.1:256481-257974+` (type `rRNA`) as the first false positive in the `oracle=nn` output and printed the corresponding `oracle_debug=NC_013946.1:256481-257974:+` switch.
- Re-running under `oracle=truth` with that debug target produced the following entries in `logs/2025-11-18_oracle_trace.log`:
  ```
  [ORACLE_DEBUG] forcing oracle/NN scoring for non-CDS target NC_013946.1:256481-257974+ type=rRNA
  [ORACLE_DEBUG] oracle truth score for NC_013946.1:256481-257974+ = 0.0 originalScore=1581.8936
  [ORACLE_DEBUG] rejecting NC_013946.1:256481-257974+ advisory=0.0 originalScore=1581.8936
  ```
  This confirms the oracle returns 0 (no matching truth entry) even though the heuristics assign a very high score, so this is a genuine FP candidate for the NN to suppress.
- After rejection, the ORF still appears in `oracle_test_11-18.gff` when `oracle=nn`, so we know the NN multiplier is not penalizing this region strongly enough; we now have a reproducible trace for further calibration work.

## Next Steps
1. Repeat the same workflow for a few more FP coordinates (covering both CDS and rRNA/tRNA) to build a short catalog of problematic cases.
2. Compare the feature vectors and raw NN outputs for those ORFs against nearby true positives to see whether the NN is saturating or whether the heuristics are overpowering the multiplier.
3. Once we understand the pattern, adjust training data / scaling and re-run the benchmarking suite.
