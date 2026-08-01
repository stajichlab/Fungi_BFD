# PFAM/HMMER hmmsearch performance A/B (T-015)

## Purpose

Answer T-015's two open questions with real timing data instead of guessing:
does copying the Pfam-A HMM database to node-local scratch reduce PFAM scan
wall time, and does enabling real MPI parallelism (`pfam_tasks>1`, currently
unused in production) meaningfully speed up `hmmsearch`. See
`todo/pfam_hmmer_performance.md` for the full writeup and recommendation; this
folder holds the scripts and raw output that data is based on.

## Status

**Status**: complete

## Datasets

- Real production Pfam-A HMM database (`db-pfam` module,
  `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`)
- Real protein set: `input/pep/Malassezia_brasiliensis_CBS_14135.proteins.fa`
  (3786 proteins — chosen as a fast, real, already-annotated genome; noted in
  `todo/pfam_hmmer_performance.md` as a caveat since it's smaller than typical)

## Key Findings

- **Scratch-copy the DB: no benefit.** Shared-storage vs. scratch-copy wall
  time is statistically identical (~1% noise) — this task is compute-bound,
  not I/O-bound. `hmmsearch` also doesn't need the pressed `.h3*` index files,
  only the raw `Pfam-A.hmm` (that's an `hmmscan` requirement).
- **`RUN_PFAM`'s existing MPI code path is broken as written** — a plain
  multi-task `srun` fails on this cluster (`CPU binding outside of job step
  allocation`) unless `--cpu-bind=none` is added. Confirmed via an independent
  compiled MPI hello-world before trusting any `hmmsearch --mpi` timing. Dormant
  today only because `pfam_tasks`/`pfam_nodes` both default to 1.
- **MPI task scaling is real but small and non-monotonic**, peaking at
  `pfam_tasks=4` (+15% vs. single-task baseline) and getting *worse* beyond
  that (`pfam_tasks=16` is slower than no MPI at all). `pfam_tasks=2` is a
  specific trap (HMMER's rank 0 dispatches rather than searches, so `-n 2` is
  really 1 worker plus coordination overhead).

Full numbers and recommendation: `../../todo/pfam_hmmer_performance.md`.
General HPCC/MPI mechanism (not PFAM-specific): `.living/learnings.md`
(2026-08-01 entry) and `~/.claude/skills/nextflow-hpcc/SKILL.md`.

## Open Questions

- Scaling curve should be re-checked against a larger/more typical genome's
  protein set before picking a production `pfam_tasks` default — 3786 proteins
  is on the small side.

## Reproducibility

Each test is a standalone `sbatch` script (not a `run.sh` pipeline, since these
are one-off timing comparisons, not a repeatable derived-output pipeline):

```bash
sbatch analysis/pfam_hmmsearch_perf/ab_scratch_vs_shared.sh   # Q1: scratch vs shared
sbatch analysis/pfam_hmmsearch_perf/mpi_sanity/test_cpubind.sh # MPI launch diagnosis
sbatch analysis/pfam_hmmsearch_perf/ab_mpi_tasks_fixed.sh      # Q2: pfam_tasks=2,4 (fixed)
sbatch analysis/pfam_hmmsearch_perf/ab_mpi_scaling.sh          # Q2: pfam_tasks=8,16
```

`ab_mpi_tasks.sh` (no `_fixed` suffix) is kept as-is — it's the *broken*
invocation (matching `RUN_PFAM`'s current code, no `--cpu-bind=none`) that
demonstrated the bug; its failure is itself part of the finding.

## Outputs

| File | Description |
|------|-------------|
| `testA_shared.timing.txt`, `testA2_shared.timing.txt` | `hmmsearch` vs. shared-storage DB (`/usr/bin/time -v`) |
| `testB_scratch.timing.txt` | `hmmsearch` vs. scratch-copied DB |
| `testC_mpi4_fixed.timing.txt`, `testD_mpi2_fixed.timing.txt` | `hmmsearch --mpi`, `-n 4`/`-n 2`, with the `--cpu-bind=none` fix |
| `testE_mpi8.timing.txt`, `testE_mpi16.timing.txt` | `hmmsearch --mpi`, `-n 8`/`-n 16` |
| `*.domtblout`, `*.tblout` | Real HMMER output per test, used for the domain-count sanity check (all identical: 7819) |
| `mpi_sanity/hello.c` | Minimal MPI hello-world used to isolate the `--cpu-bind` launch bug from HMMER's own MPI code |
| `mpi_sanity/test_cpubind.sh`, `mpi_sanity/test_mpi_launch.sh` | The diagnostic jobs that found/confirmed the `--cpu-bind=none` fix |
