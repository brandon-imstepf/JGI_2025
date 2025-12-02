# CallGenes NN Tolerance & Training Updates (2025-11-13)

## Summary
We added a tolerance-aware truth lookup for CallGenes so that both the oracle workflow and the BTSF training pipeline can treat ORFs that differ only slightly from the NCBI reference (e.g., alternate start codons) as positives. This keeps Brian’s original scripts intact and confines new behavior to CallGenesHelper and associated Slurm wrappers. We also introduced a chunk-wise TSV shuffler and wired it into the existing training Slurm jobs so future neural nets see randomized batches without manual preprocessing.

## Key Changes
1. **CallGenesHelper updates** (`bbmap/current/prok/CallGenesHelper.java`)
   - Canonicalizes contig IDs and indexes all reference CDS entries once per run.
   - Adds `label_tolerance=<int>` (default 3 nt) so truth checks allow ±tolerance on both start and stop, applied everywhere we compare against the reference (oracle, logging, FN accounting, vector labels).
   - Exposes the tolerance-sensitive lookup to `CallGenesOracle`, so `oracle=truth` now honors the same notion of “match” used during training.

2. **CallGenesOracle** (`bbmap/current/prok/CallGenesOracle.java`)
   - Accepts a simple `TruthLookup` callback instead of owning the set, letting CallGenesHelper supply the canonical/tolerance-aware predicate.

3. **Training data shuffle** (`python/shuffle_tsv.py` + `slurm/train_total_cycles_*.slurm`)
   - New chunk shuffler randomizes huge TSVs by shuffling fixed-size batches to temporary files, keeping memory usage low.
   - Training Slurm scripts now call the shuffler before `train.sh`, ensuring positives/negatives are mixed prior to each run without altering Brian’s core tooling.

4. **Dataset builder variants** (`slurm/run_builder_tol3.slurm`, `slurm/run_builder_orig_tol3.slurm`)
   - Copies of the existing builders that pass `label_tolerance=3` to `btsf_new.sh`, ready for launch once the code is approved. Original scripts remain unchanged for comparison/backups.

## Label Tolerance Concept
- **What it does**: Accepts ORFs whose start and/or stop differ by up to `label_tolerance` nucleotides from the reference when determining whether they are “true genes.” The comparison is strand-aware and uses canonicalized contig IDs.
- **Default**: `label_tolerance = 3`, which covers single-codon shifts caused by alternate start sites or small alignment differences. Setting `label_tolerance=0` restores the old exact-match behavior.
- **CLI usage**: Available anywhere you run CallGenesHelper (e.g., `callgenes.sh … label_tolerance=5`). The value propagates to the oracle, the BTSF pipeline (via the new builder scripts), and reporting/logging.
- **Why it matters**: Validation_NCBI_10 (my baseline) contains many CDS that differ by a codon or two relative to the training references. With zero tolerance those ORFs were labeled as negatives, so the NN learned to down-rank exactly the cases we need to recover. Tolerance labeling ensures near-misses stay positives during training/inference, enabling the oracle to reach higher recall and letting the NN focus on genuine false positives.

## Next Steps
1. Review the code and the new `label_tolerance` flag. If the default of 3 nt looks good, regenerate the master TSVs via `run_builder_tol3.slurm` (or the `_orig` variant) before retraining.
2. Run the updated training Slurm jobs (once approved) so the NN consumes the shuffled TSVs.
3. Re-benchmark (`run_nn_oracle_validation.sh`, `run_benchmarking_all.sh`) to quantify improvements before merging upstream.
