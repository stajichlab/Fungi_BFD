# Reclaiming normalized RNA-seq reads after Trinity + training are done

Status: **tooling added, not yet applied** (2026-06-20). `rnaseq_reads/` holds ~787 GB of
normalized FASTQ. Once a species' transcript work is fully consumed, those reads are
disposable. This note documents the safe-deletion predicate, the stub mechanism, how it
interacts with `funannotate train`/`predict`, and how to recover.

See also: `doc/README_RNASeq_handling.md` (Trinity caching / staleness model),
`scripts/audit_rnaseq_reclaim.py`, `scripts/reclaim_rnaseq_reads.py`.

---

## TL;DR

- Normalized reads `rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz` are consumed
  **twice**: once by `RNASEQ_PREPARE` (→ per-species Trinity) and once per *strain* by
  `FUNANNOTATE_TRAIN` (→ PASA / BAM / stringtie). Both must finish before reads are
  disposable.
- **Safe-deletion predicate** (per `species_tag`):
  1. `rnaseq_data/<species_tag>.trinity-GG.fasta` exists and is **non-empty**, AND
  2. **every** assembly mapped to that `species_tag` in `samples.csv` has a non-empty
     `genome_annotation_training/<out>/training/funannotate_train.pasa.gff3`.
- We **stub** (truncate to zero bytes), not `rm`. `SRA_FETCH` uses `storeDir rnaseq_reads/`;
  deleting its outputs re-triggers an expensive re-download, while an empty file keeps
  `storeDir` satisfied and reads as "no data" everywhere (all branches gate on `.size() > 0`).
- The stub keeps the **original mtime**, so it can never look "newer than the GBK" and trip
  `staleRnaseq()`.
- A per-species provenance manifest `rnaseq_data/<species_tag>.reads_reclaimed.json` records
  sha256/bytes/mtime of each read, the Trinity sha256, the trained strains, and the SRA
  accessions used.

---

## Why stubs and not `rm`

`SRA_FETCH` / `SRA_FETCH_SE` declare `storeDir "${launchDir}/rnaseq_reads"` and are skipped
only when their **declared output files exist**. A bare `rm` makes the next run treat the
species as a `storeDir` miss and re-download + bbnorm + fastp from scratch. A zero-byte stub
satisfies `storeDir` and is interpreted as "no reads" by every downstream `.size() > 0`
branch — no churn, no re-download. (`reclaim_rnaseq_reads.py --mode delete` exists but warns.)

---

## How each process interprets a reclaimed (stubbed) species

`species_tag` = `SPECIES` with quotes stripped and whitespace runs → `_`
(`sed -E 's/[[:space:]]+/_/g'` in bash; `replaceAll(/\s+/, '_')` in Groovy).

| Process | Reads an input? | Behavior with empty stub reads (Trinity intact) |
|---|---|---|
| `SRA_FETCH` / `SRA_FETCH_SE` | produces them | `storeDir` hit (stub exists) → **skipped**, no re-download |
| `RNASEQ_PREPARE` | yes (repr only) | `storeDir` hit on the non-empty Trinity FASTA → **skipped** |
| `FUNANNOTATE_TRAIN` | yes (`path r1/r2/se`) | `pasa.gff3` exists + reads empty → RETRAIN guard `[ -s r1 ]` is false → **"Training already complete; skipping"** (`funannotate.nf:1141`) |
| `FUNANNOTATE_PREDICT` | **no** | does not stage reads at all; re-derives staleness from on-disk files, each guarded by `[ -s "$f" ]` — empty stubs are ignored; only the (untouched) Trinity could trip it. **No change.** |

### `funannotate predict` in detail — does it need changes? **No.**

`FUNANNOTATE_PREDICT` (`funannotate.nf:1293`) takes only metadata + `genome_fa` as input;
it never stages the normalized reads. It consumes the persistent `training/` directory (via a
symlink, `funannotate.nf:1366-1369`) — i.e. `pasa.gff3`, the BAM, and the stringtie GTF, all
of which our predicate guarantees already exist and which the reclaim never touches.

Its only contact with `rnaseq_reads/` is the staleness re-derivation at
`funannotate.nf:1337-1344`:

```bash
for f in .../${SPECIES_TAG}_norm_R1.fastq.gz \
         .../${SPECIES_TAG}_norm_SE.fastq.gz \
         .../${SPECIES_TAG}.trinity-GG.fasta; do
    if [ -s "$f" ] && [ "$f" -nt "$SKIP_GBK" ]; then STALE=1; fi
done
```

Every term is guarded by `[ -s "$f" ]` (**non-empty**). A zero-byte stub fails that guard and
is ignored, so a reclaimed species with a current GBK still short-circuits to
`predict.done`. The only file that can legitimately trip staleness is the Trinity FASTA —
which reclaim never modifies. Mtime preservation on the stub is belt-and-suspenders: even a
fresh mtime wouldn't matter because the size guard already excludes empty files.

**Conclusion: `FUNANNOTATE_PREDICT` needs no code changes.** The same holds for
`FUNANNOTATE_TRAIN`'s skip path for already-trained strains.

---

## The one real risk: a NEW strain added later

If a strain is added to `samples.csv` for an **already-reclaimed** species, its
`FUNANNOTATE_TRAIN` has no `pasa.gff3` yet, so it proceeds to assembly. With a non-empty
shared Trinity but empty stub reads it falls into the "PASA only, no reads" branch
(`funannotate.nf:1232`), which calls `funannotate train --trinity ... --left_norm <empty>
--right_norm <empty>` — likely to fail.

Mitigations:
- The reclaim predicate only fires once **all current strains** are trained, so this can only
  arise from a *future* `samples.csv` addition.
- The manifest's `strains_trained` list lets a pre-flight check detect a reclaimed species
  that has gained a new untrained strain.
- To genuinely refresh such a species, follow the recovery procedure below (restore reads).

---

## Workflow

```bash
# 1. Audit — read-only classification + reclaimable totals
python3 scripts/audit_rnaseq_reclaim.py
python3 scripts/audit_rnaseq_reclaim.py --list-safe       # itemize
python3 scripts/audit_rnaseq_reclaim.py --list-blocked    # Trinity-ready, strains pending

# 2. Validate the round-trip on ONE species (writes manifest + stubs)
python3 scripts/reclaim_rnaseq_reads.py --only Morchella_importuna --apply
cat rnaseq_data/Morchella_importuna.reads_reclaimed.json

# 3. Mass reclaim (dry-run is the default; --apply to act)
python3 scripts/reclaim_rnaseq_reads.py            # dry-run
python3 scripts/reclaim_rnaseq_reads.py --apply    # act
```

`reclaim_rnaseq_reads.py` re-verifies the full predicate at apply time (never trusts a stale
audit), writes the manifest **before** stubbing (so provenance survives an interrupt), and is
idempotent (skips species that already have a manifest unless `--force`).

As the training backlog drains, previously **blocked** species (Trinity built, some strains
untrained) automatically become eligible on the next run — the recurring savings exceed the
first-pass total.

### Snapshot (2026-06-20)

| Class | Species | Bytes |
|---|---|---|
| SAFE now | 311 | 163.6 GB |
| Blocked (Trinity ready, strains pending) | 544 | 292.6 GB |
| No Trinity yet | 5,089 | — |
| Empty Trinity (no RNA-seq) | 1,296 | ~0 |

---

## Recovery: forcing a true refresh of a reclaimed species

Reads are gone (stubbed), so a refresh means re-fetching. To pull new SRA data through:

1. Remove the Trinity cache and the manifest:
   `rm rnaseq_data/<species_tag>.trinity-GG.fasta rnaseq_data/<species_tag>.reads_reclaimed.json`
2. Remove the stub reads so `SRA_FETCH`'s `storeDir` re-downloads:
   `rm rnaseq_reads/<species_tag>_norm_*.fastq.gz`
3. Remove the training folders for the affected strains so PASA rebuilds:
   `rm -rf genome_annotation_training/<out>/training`
4. Re-run the pipeline. `SRA_FETCH` re-fetches → `RNASEQ_PREPARE` rebuilds Trinity →
   `FUNANNOTATE_TRAIN` rebuilds → `FUNANNOTATE_PREDICT` re-predicts (new Trinity mtime trips
   `staleRnaseq`).

The accession list in the manifest documents what the original Trinity was built from.
