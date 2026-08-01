# PFAM/HMMER scan performance — scratch-copy DB and MPI task scaling

**Status**: real A/B timing data collected 2026-08-01 (see §Results below) — both
open questions answered. Remaining work is implementation, not investigation:
fix `RUN_PFAM`'s MPI `--cpu-bind` bug and decide whether to enable
`pfam_tasks>1` in production given the modest (+15%) real-world gain.
**Added**: 2026-08-01, during T-014 (genome_stats/function storage reorg) work.
**Relationship to T-014**: separate concern — T-014 changed *where* PFAM's merged
output lands (bucketed storeDir, then Parquet staging via `MERGE_PFAM`), not how
the PFAM scan itself runs. This item is about the scan's own wall time, raised as
an adjacent open question while touching `MERGE_PFAM`/`pfamtbl_to_long.py` during
that work, not a blocker for or dependency of T-014.

## Confirmed starting facts (checked against real code/config, 2026-08-01)

- `nextflow/modules/BFD/PFAM/main.nf` (`RUN_PFAM`) already runs `hmmsearch`
  (query=Pfam-A HMM DB, target=protein set) — **not** `hmmscan`. This is the
  generally-recommended direction for scanning against a large HMM library, so
  there is no "switch to hmmscan" question here; a switch would likely be a
  performance *regression* (hmmscan processes one protein against the whole
  library at a time — worse cache locality for a big library — which is exactly
  why hmmsearch is preferred) on top of a real parsing/output-format migration
  cost. **Dropped from scope per explicit user decision (2026-08-01): stay with
  hmmsearch, focus on better use of MPI instead.**
- `$PFAM_DB` (from `module load db-pfam`) resolves to
  `/srv/projects/db/pfam/2026-01-27-Pfam38.2/` — shared/network storage, not
  node-local scratch. `Pfam-A.hmm` is 2.1 GB; the pressed index files
  (`.h3f`/`.h3i`/`.h3m`/`.h3p`) total another ~1.6 GB (already `hmmpress`-ed
  once, so no per-task pressing overhead — that's not part of what a
  scratch-copy would save).
- `params.pfam_nodes`/`params.pfam_tasks` both default to **1**
  (`conf/profile_BFD.config`). The MPI code path (`hmmer/3.4-mpi` module,
  `srun -N ${pfam_nodes} -n ${pfam_tasks} --mpi hmmsearch`) only activates when
  `pfam_tasks > 1` — which nothing in the current default config sets. **In
  normal production runs today, PFAM is not actually running in MPI mode at
  all** — every task is a single-process `hmmsearch --cpu 4` (the `pfam` label's
  `cpus = 4`), not a multi-task MPI job. "Better use of MPI" means: does turning
  on real `pfam_tasks`/`pfam_nodes` MPI parallelism beat the current single-task
  4-thread baseline, and by how much relative to the added node/task scheduling
  overhead?

## Open questions (narrowed scope)

1. **Scratch-copy the Pfam-A HMM database?** Would copying the ~3.7 GB Pfam-A
   HMM DB (pressed files included) to node-local scratch, instead of reading it
   from `/srv/projects/db/pfam/...` per-task, meaningfully reduce PFAM scan wall
   time at this cluster's scale? Don't assume — profile actual I/O wait vs.
   compute time in a real `hmmsearch` task first. Cheapest possible test: time
   one real `hmmsearch` run against `$PFAM_DB` in place vs. against a
   scratch-copied version, same genome, same node.
2. **MPI task/node scaling for `hmmsearch`?** Does enabling `pfam_tasks > 1`
   (real MPI `hmmsearch`, currently unused in production) reduce wall time
   meaningfully over the current single-task 4-thread baseline? If so, at what
   `pfam_tasks`/`pfam_nodes` combination does the benefit taper off relative to
   scheduling/startup overhead? This replaces the original "chunking" question —
   chunking the protein set for separate parallel tasks is one way to get more
   parallelism, but tuning `pfam_tasks`/`pfam_nodes` on the existing MPI code
   path is the more direct lever already built into `RUN_PFAM` and should be
   tried first.

## Results (real A/B timing, 2026-08-01)

Real `hmmsearch` runs on `highclock`, real Pfam-A DB, real protein set
(`Malassezia brasiliensis` CBS 14135, 3786 proteins — on the small side, see
caveat below). Full scripts/raw output in `analysis/pfam_hmmsearch_perf/`.

**Q1 — scratch-copy the DB: no benefit, don't implement.** Copying the ~3.7GB
DB to node-local scratch took 2.7s (fast NFS→NVMe), but wall time was
statistically identical reading from shared storage vs. scratch (5:00.35 vs
5:02.54 vs a 5:02.84 repeat — within ~1% noise). This task is compute-bound,
not I/O-bound. Also confirmed: `hmmsearch` runs fine against just the raw
`Pfam-A.hmm` file — it doesn't need the pressed `.h3f`/`.h3i`/`.h3m`/`.h3p`
index files at all (those are for `hmmscan`'s per-target random access, not
`hmmsearch`'s sequential HMM-as-query read). Not relevant now since scratch-copy
isn't happening, but worth remembering if this ever gets revisited alongside
`hmmscan`.

**Q2 — MPI task scaling: real gain, but small and non-monotonic; found a real
bug along the way.** First attempt at MPI timing — using `RUN_PFAM`'s *exact*
current invocation pattern (`srun -N ${pfam_nodes} -n ${pfam_tasks} --mpi
hmmsearch`) — failed immediately: `srun: error: CPU binding outside of job
step allocation ... Unable to satisfy cpu bind request`. **This means
`RUN_PFAM`'s existing MPI code path is currently broken as written** — it's
just never been hit because `pfam_tasks`/`pfam_nodes` both default to 1.
Root-caused (via a minimal compiled MPI hello-world, not guesswork) to a
`--cpu-bind` default vs. cgroup/allocation-mask mismatch for multi-task job
steps on this cluster, independent of `--mpi=<plugin>` choice. Fix: add
`--cpu-bind=none` to the `srun` call. See `.living/learnings.md` (2026-08-01)
and `~/.claude/skills/nextflow-hpcc/SKILL.md` for the full mechanism —
this is a general HPCC finding, not PFAM-specific.

With that fix applied, real scaling data vs. the current single-task 4-thread
baseline (5:00.35):

| `pfam_tasks` | real workers (rank 0 dispatches, doesn't search) | wall time | vs. baseline |
|---|---|---|---|
| 1 (current default) | — | 5:00.35 | baseline |
| 2 | 1 | 9:15.56 | **-85% (worse)** |
| 4 | 3 | 4:16.58 | **+15% faster (best)** |
| 8 | 7 | 4:49.62 | +4% faster |
| 16 | 15 | 5:11.42 | -3% (worse than baseline) |

Non-monotonic: peaks at `pfam_tasks=4`, degrades beyond that as coordination
overhead outgrows available parallel work. `pfam_tasks=2` is a specific trap —
1 real worker plus dispatch overhead is worse than no MPI at all. Domain
counts (7819) identical across every configuration — differences are pure
wall-time, not correctness.

**Caveat**: 3786 proteins is a smallish genome. The scaling curve (where it
peaks, how much headroom exists) should be re-checked against a large/typical
genome before picking a production default — a bigger protein set may sustain
more real workers before coordination overhead dominates.

## Recommendation

1. **Don't implement scratch-copy staging.** No measurable benefit, real cost
   (dev time, extra failure surface).
2. **If enabling `pfam_tasks>1` in production**: first fix `RUN_PFAM`'s
   `mpi_launch`/`srun` construction (`nextflow/modules/BFD/PFAM/main.nf`) to
   add `--cpu-bind=none` — without it, turning on `pfam_tasks>1` today would
   crash every PFAM task, not just underperform. Then default to
   `pfam_tasks≈4`, not a larger number — this data shows more tasks make it
   *worse*, not better, past that point. Re-validate against a large genome's
   protein set before committing to a specific default, since 3786 proteins is
   on the small end of what this pipeline processes.
3. **Given the modest realistic gain (+15% at best, on a small genome)**,
   weigh implementation/maintenance cost against the win before prioritizing
   this — it's a real but not dramatic speedup, not a step-change.

## Relevant existing code

- `nextflow/modules/BFD/PFAM/main.nf` (`RUN_PFAM`) — the scan process itself;
  `hmmsearch --cut_ga --noali --cpu ${task.cpus} $PFAM_DB/Pfam-A.hmm ${proteins}`,
  MPI-conditional on `params.pfam_tasks`.
- `nextflow/conf/profile_BFD.config` — `pfam_nodes`/`pfam_tasks` params (both
  default 1), `withLabel: 'pfam'` resource block (`cpus=4`, `queue=highclock`,
  `clusterOptions` computed from `pfam_tasks`/`pfam_nodes`).
- `nextflow/modules/BFD/MERGE_PFAM/main.nf` / `scripts/pfamtbl_to_long.py` — the
  merge step touched by T-014 (#27); downstream of the scan itself, not the
  performance question, but where the merged PFAM table is produced.

**Tags**: performance, pfam, hmmer, hmmsearch, mpi, scratch-partition, task-scaling, future-investigation
