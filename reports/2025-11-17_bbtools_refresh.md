# BBTools Refresh & File Migration (2025-11-17)

## Goal
- Stand up a clean `bbmap_fresh/` clone of the current BBTools repo, then transplant our CallGenes integration points and auxiliary tooling so future work happens against a modern baseline.

## Actions
1. **Core CallGenes files copied**
   - `bbmap/current/prok/CallGenes.java`
   - `bbmap/current/prok/CallGenesHelper.java`
   - `bbmap/current/prok/CallGenesOracle.java`
   - `bbmap/current/prok/GffToTsv.java`

2. **NN training/runtime assets moved**
   - Entire `bbmap/current/ml/` directory (custom BBNet/CNN trainers plus helpers) synced into `bbmap_fresh/current/ml/`.
   - Slurm job scripts copied via `rsync -a bbmap/current/slurm/ bbmap_fresh/current/slurm/`.

3. **Benchmarking/oracle wrappers restored**
   - `run_benchmarking_all.sh`, `run_benchmarking_single.sh`, `run_callgenes_modes.sh`
   - `run_nn_oracle_validation.sh`, `run_nn_prefilter_validation.sh`
   - `run_validation_cutoff.sh`, `run_validation_strength.sh`

All copies were one-way overwrites (`cp`/`rsync`) from the legacy `bbmap/` tree into `bbmap_fresh/`, so the new checkout now contains our CallGenes helper stack, NN trainers, Slurm launchers, and validation scripts exactly as they existed before.

## Next Steps
1. Review the migrated files inside `bbmap_fresh/` (especially CallGenesHelper and the run\_*.sh scripts) to ensure there were no upstream API changes we need to reconcile.
2. Update `README_git.md` with the new workflow (archive → clone → copy helper files) so future refreshes follow the same playbook.
