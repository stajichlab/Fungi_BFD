# Functional pipeline production rollout — SwissProt + MEROPS pilot → full run

| | |
|---|---|
| ID | none (the registry's T-022 is GeneMark ET; this plan was labelled T-022 by mistake) |
| Status | pending (2026-08-10) |
| Priority | high |
| Registry | no entry |

## Goal
Roll the newly containerized, SLURM-validated SwissProt + MEROPS modules out to
the full 10,943-genome BFD set, staged as decided on 2026-08-10:

1. **Pilot**: `run_merops` + `run_swissprot` only, over the first few hundred
   samples — measure real scale timing and resource behavior.
2. **Gate** on the pilot, then the full functional run
   (`nextflow/run_functional.sh`, which enables CAZy, MEROPS, SwissProt, SignalP,
   IDP/AIUPred, WoLF-PSORT, BUSCO genome+pep).

## Prior evidence (already committed)
- Two real-SLURM smokes, single real genome (`Moesziomyces_antarcticus_T-34`):
  - `smoke_swissprot/` — diamond engine, green (`tables_diamond/` snapshot kept).
  - `smoke_blastp/` — blastp 2.16.0 engine, green.
- Commits: `912054a` (modules/bin/script containerization), `5c1c5c1`
  (run_functional.sh wiring), `5f8a107` (schema), `b6ecc70` (learnings).
- Benchmark (8 threads, `-k`/`-max_target_seqs 20`, `-e 1e-5`):
  diamond 2.2.5 `--sensitive` **2m19s** (68,238 rows) vs blastp 2.16.0
  `-seg yes -soft_masking true` **~56 min** (80,502 rows). **Diamond stays default.**

## Pilot run details
- Launch from the project root (where production `samples.csv` lives), mirroring
  the smoke scope with everything but MEROPS/SwissProt off:
  ```
  sbatch nextflow/run_functional.sh \
      --n_test 200 \
      --run_setup true \
      --run_pfam false --run_cazy false \
      --run_merops true --run_swissprot true \
      --run_signalp false --run_tmhmm false --run_targetp false \
      --run_idp false --run_wolfpsort false --run_predgpi false \
      --run_busco_genome false --run_busco_pep false
  ```
- `run_build_duckdb` stays off — `build_BFD_duckDB.sh` globs the **full** table
  set (`asm_stats.parquet`, ...) and fails on a partial run (expected).
- All inputs/workdir on shared `/bigdata` (never `/scratch` — node-local NVMe).

## Pilot metrics to capture (the gate)
- Per-genome `RUN_SWISSPROT`/`RUN_MEROPS` wall-time distribution at scale
  (vs 2m19s single-genome benchmark), SLURM concurrency on `epyc`, queueing.
- `tables/{merops,swissprot,swissprot_annot}.parquet` row counts scale with
  sample count; `swissprot_acc→accession` join resolution stays ~100%.
- Note: to A/B engines at scale, use separate outdir trees — `storeDir` on
  `RUN_SWISSPROT` means a changed `swissprot_search` won't re-run if the
  per-species output already exists (learned 2026-08-10).

## Post-pilot
- Full functional run via `nextflow/run_functional.sh` (all tools above) across
  10,943 genomes; `-resume` against the pilot work dir reuses completed tasks.
- Cross-cutting: `asm_stats.parquet` (T-001) is a separate pipeline, not part of
  the functional run.
