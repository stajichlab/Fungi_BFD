# PFAM/HMMER scan performance — scratch-copy DB, hmmscan vs hmmsearch, chunking

**Status**: future-investigation, not yet started.
**Added**: 2026-08-01, during T-014 (genome_stats/function storage reorg) work.
**Relationship to T-014**: separate concern — T-014 changed *where* PFAM's merged
output lands (bucketed storeDir, then Parquet staging via `MERGE_PFAM`), not how
the PFAM scan itself runs. This item is about the scan's own wall time, raised as
an adjacent open question while touching `MERGE_PFAM`/`pfamtbl_to_long.py` during
that work, not a blocker for or dependency of T-014.

## Open questions

1. **Scratch-copy the Pfam-A HMM database?** Would copying the Pfam-A HMM DB to
   node-local scratch (instead of reading it from shared/network storage per-task)
   meaningfully reduce PFAM scan wall time at this cluster's scale? Don't assume —
   profile actual I/O wait vs. compute time in a real PFAM task first.
2. **`hmmscan` vs `hmmsearch`?** `hmmscan` (query=protein, target=HMM DB) vs
   `hmmsearch` (query=HMM DB, target=protein) — general HMMER guidance favors
   `hmmsearch` for large HMM databases like Pfam-A, since it avoids the
   per-target DB reformatting/indexing overhead `hmmscan` incurs. Needs
   verification against this pipeline's *actual current implementation*
   (`RUN_PFAM` module) and real timing data, not assumed from general guidance.
3. **Chunking?** Is splitting the protein set into batches for parallelism worth
   the added complexity, or is plain `hmmsearch` + existing MPI parallelism
   already fast enough without it? This is the explicit open decision point —
   profile first, then decide.

## Recommended next step

Gather real timing/profiling data on current PFAM task runs before implementing
any of the above — wall time breakdown between DB I/O, `hmmscan`/`hmmsearch`
compute, and per-task startup overhead. Mirror the approach already used
elsewhere in this repo to ground performance decisions in data rather than
guessing:
- `analysis/funannotate_predict_stage_timing/` (FUNANNOTATE_PREDICT stage timing)
- `analysis/nextflow_memory_profile/` (memory/CPU profiling)

That data should drive whether scratch-copy, the `hmmsearch` switch, and/or
chunking are worth doing — don't implement all three speculatively.

## Relevant existing code (starting points, not yet audited for this)

- `nextflow/modules/BFD/RUN_PFAM/main.nf` (or equivalent) — the actual PFAM scan
  process; check whether it currently calls `hmmscan` or `hmmsearch`, and how
  the Pfam-A DB path is referenced (shared storage vs. staged).
- `nextflow/modules/BFD/MERGE_PFAM/main.nf` / `scripts/pfamtbl_to_long.py` — the
  merge step touched by T-014 (#27); downstream of the scan itself, not the
  performance question, but where the merged PFAM table is produced.

**Tags**: performance, pfam, hmmer, hmmscan, hmmsearch, scratch-partition, chunking, mpi, future-investigation
