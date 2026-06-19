# RNA-seq handling, Trinity caching, and staleness in `funannotate.nf`

Status: **current behavior, intentionally left as-is** (2026-06-18). This note documents
how RNA-seq reads, the Trinity transcript assembly, training, and prediction interact so
future decisions can be made deliberately. The primary concern driving this design is
**avoiding unnecessary churn**: once a Trinity transcript assembly exists for a species, we
do *not* want a new SRA accession appearing in the query to silently trigger an expensive
re-assembly.

See also: `doc/` plans, and the process definitions in `nextflow/funannotate.nf`.

---

## TL;DR

- The Trinity transcript FASTA for a species is built **once** and cached at
  `rnaseq_data/<species_tag>.trinity-GG.fasta`. It is **never regenerated automatically**,
  even when new SRA data becomes available. This is the churn guard you want.
- A "stale" signal (RNA-seq/Trinity newer than the prediction GBK) re-runs **train** and
  **predict**, but **not** the Trinity assembly. So a stale rerun refreshes the prediction
  **using the already-cached (possibly older) Trinity**.
- To force a genuine end-to-end refresh from new reads you must **manually delete**
  `rnaseq_data/<species_tag>.trinity-GG.fasta` (and usually the training folder). Nothing in
  the pipeline does this on its own — by design.

---

## The pieces and where their outputs live

| Process | Output location | Cache mechanism | Cost |
|---|---|---|---|
| `SRA_QUERY` / `SRA_QUERY_BATCH` | `rnaseq_reads/sra_query/<species_tag>.sra_query.csv` | `storeDir` / publishDir, reuses existing CSV | cheap (network only) |
| `SRA_FETCH` / `SRA_FETCH_SE` | `rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz` | `storeDir rnaseq_reads/` | expensive (download + bbnorm + fastp) |
| `RNASEQ_PREPARE` | `rnaseq_data/<species_tag>.trinity-GG.fasta` | `storeDir rnaseq_data/` | **very expensive** (genome-guided Trinity) |
| `FUNANNOTATE_TRAIN` | `genome_annotation_training/<out>/training/` | manual skip on `funannotate_train.pasa.gff3` | expensive (PASA; full Trinity only in fallback) |
| `FUNANNOTATE_PREDICT` | `genome_annotation/<out>/predict_results/` | direct-to-target + on-disk GBK guard | expensive (Augustus/GeneMark/EVM) |

`species_tag` = `SPECIES` with whitespace runs collapsed to `_` (`replaceAll(/\s+/, '_')`).
This transform is used identically in Groovy (`staleRnaseq`) and in the predict bash guard.

---

## Who writes the Trinity FASTA, and when

There is exactly **one** writer of `rnaseq_data/<species_tag>.trinity-GG.fasta`:
`RNASEQ_PREPARE`. Because it declares `storeDir "${launchDir}/rnaseq_data"`, **the entire
task is skipped whenever that file already exists.** That single fact is the churn guard.

When `RNASEQ_PREPARE` *does* run (i.e. the cached FASTA is absent), it has two internal
paths:

1. **Recover-from-training (cheap):** if
   `genome_annotation_training/<out>/training/funannotate_train.pasa.gff3` already exists, it
   copies the existing `training/trinity.fasta` out into `rnaseq_data/` — no re-assembly.
2. **Fresh assembly (expensive):** otherwise it runs
   `funannotate train --stop_after_trinity` in `$SCRATCH`, copies the resulting
   `trinity.fasta` to `rnaseq_data/`, and deletes the scratch dir.

`FUNANNOTATE_TRAIN` **never** writes to `rnaseq_data/`. It consumes the shared Trinity via
`--trinity <rnaseq_data FASTA>` (PASA-only path) and writes results into the training
folder. Its only cleanup is `rm -rf hisat2` and `rm -rf trinity_gg`. (There is a fallback
"full train, no shared Trinity" path that assembles Trinity inside the training folder, but
it still does not copy it back to `rnaseq_data/`.)

**Implication:** new SRA accessions discovered by a later `SRA_QUERY` do not propagate into
Trinity. Once `rnaseq_data/<species_tag>.trinity-GG.fasta` exists, the transcript assembly is
frozen until a human removes it.

---

## The staleness signal (`staleRnaseq`)

`staleRnaseq(out, species)` (`nextflow/funannotate.nf`) returns true when the prediction GBK
exists **and** any of the following is newer than the GBK:

- `rnaseq_reads/<species_tag>_norm_R1.fastq.gz`
- `rnaseq_reads/<species_tag>_norm_SE.fastq.gz`
- `rnaseq_data/<species_tag>.trinity-GG.fasta`

It gates two things:

- **`FUNANNOTATE_TRAIN`**: stale species are added to `train_todo` and re-run. NOTE: the
  `rm -rf training` inside the retrain branch is currently **commented out**, so a stale
  retrain resumes the existing training folder rather than rebuilding it from scratch.
- **`FUNANNOTATE_PREDICT`**: stale genomes are scheduled, and the process re-derives
  staleness from the same on-disk timestamps. If stale, it clears `predict_results/` and
  `predict_misc/` and runs a **fresh** prediction.

`staleRnaseq` does **not** gate `RNASEQ_PREPARE`, `SRA_FETCH`, or `SRA_QUERY`. Their
`storeDir` caches are invalidated only by deleting the cached files.

Because Trinity's own mtime is one of the staleness inputs, the only realistic way for it to
become "newer than the GBK" is for a human to replace/regenerate it — which then correctly
cascades into retrain+repredict.

---

## What actually happens on a stale rerun

Scenario: prediction already exists; normalized reads (or Trinity) are newer than the GBK.

| Step | Behavior | Effect on transcripts |
|---|---|---|
| `SRA_FETCH` | skipped (storeDir hit) | reads unchanged |
| `RNASEQ_PREPARE` | **skipped (storeDir hit)** | **Trinity unchanged — OLD assembly** |
| `FUNANNOTATE_TRAIN` | re-runs (stale) using `--trinity <old FASTA>`, resumes training dir | training refreshed against OLD Trinity |
| `FUNANNOTATE_PREDICT` | re-runs **fresh** (clears predict_results/predict_misc) | prediction refreshed against OLD Trinity |

So a stale rerun gives you a **fresh prediction built on the cached Trinity**. It does *not*
incorporate transcripts from any newly available SRA runs. This is the deliberate trade:
prediction can be cheaply refreshed (e.g. after a parameter change), without paying for a
Trinity re-assembly each time the SRA query turns up another accession.

---

## How to force a true end-to-end refresh (manual, intentional)

When you actually want new SRA data to flow all the way through for a species:

1. Remove the Trinity cache:
   `rm rnaseq_data/<species_tag>.trinity-GG.fasta`
2. Remove the training folder so PASA rebuilds against the new transcripts:
   `rm -rf genome_annotation_training/<out>/training`
   (optionally also clear the normalized reads in `rnaseq_reads/` if you want new accessions
   re-downloaded and re-normalized — otherwise the existing normalized reads are reused).
3. Re-run the pipeline. `RNASEQ_PREPARE` rebuilds Trinity → `FUNANNOTATE_TRAIN` rebuilds
   training → `FUNANNOTATE_PREDICT` re-predicts (the new Trinity mtime trips `staleRnaseq`).

There is intentionally no automatic trigger for this cascade.

---

## Open questions / decisions to revisit later

- **Selective refresh policy.** Is "refresh prediction on cached Trinity" the right default,
  or should certain events (e.g. a substantially larger/newer SRA accession) justify a
  Trinity rebuild? Today there is no heuristic; the cache is all-or-nothing per species.
- **Accession set drift.** `SRA_QUERY` may find new accessions over time, but the Trinity
  assembly reflects whatever accessions were used at first build. There is no record in
  `rnaseq_data/` of *which* accessions went into a given Trinity FASTA. If accession
  provenance matters for decisions, consider emitting a sidecar manifest
  (e.g. `rnaseq_data/<species_tag>.trinity-GG.accessions.txt`) at build time.
- **Train-side staleness.** The commented-out `rm -rf training` in `FUNANNOTATE_TRAIN` means
  a stale retrain resumes rather than rebuilds. Combined with the frozen Trinity, a stale
  rerun's training refresh is partial. Decide whether stale should truly rebuild training.
- **A churn-safe "consider new data" mode.** If desired, a future opt-in flag could compare
  the current `SRA_QUERY` accession set against a stored manifest and only then invalidate
  the Trinity cache — giving controlled refresh without churn on every run.
