# Design: handling representatives whose shared Trinity cannot train (RNA-seq)

Status: UNDER REVIEW — initial forensics COMPLETE (2026-08-26); verdict below
resolves the original Hypothesis A vs B question. Not yet implemented. The
document lays out the root-cause analysis of the
`Saccharomyces_cerevisiae_KCTC_13826BP` train failure, evaluates how to encode
"a representative too far from the transcriptome reference → either build its
own RNA-seq or skip RNA-seq-based training (ab-initio-only)" into the Nextflow
strategy, and defines the evidence package and question set for external review
by a bioinformatics expert and a software/workflow expert.

## Forensics verdict (2026-08-26): the "representative" is a misidentified genome

The initial framing ("divergent S. cerevisiae strain") is **superseded**. The
forensics (section below) show:

- `GCA_026225675.1_ASM2622567v1` ("KCTC 13826BP", submitted as
  "*Saccharomyces boulardii* (nom. inval.)") is **not S. cerevisiae at all**.
  NCBI's own ANI taxonomy check fails for it (`match_status: mismatch`,
  `is_atypical: true`, notes `unverified source organism`): best ANI match is
  **Nakaseomyces glabratus (Candida glabrata) at 99.93% ANI / 96.9% coverage**;
  match to S. cerevisiae is 87.2% ANI at 0.95 coverage (fail).
- Independent alignment forensics agree with NCBI:
  - Trinity (reads SRR3491422 = **BY4743, an S288C-lab strain**) aligns
    **99.8%** (4,259 high-quality ≥80%-cov) to S288C R64, but only **0.3%**
    (17 high-quality) to KCTC's genome.
  - Whole-genome KCTC vs S288C: only **5.1%** of KCTC's 12.9 Mb aligns, **zero**
    contigs ≥75% covered (strain-level divergence would be ≥95% covered at
    ≥97% identity; genuine S. cerevisiae sibling Y10 aligns 100.5% to S288C).
  - KCTC is absent from `tables/BUSCO.csv.gz`; all ~1,308 in-scope sibling
    S. cerevisiae rows have blank ANI to it (`reuse_eligible=False`) — the
    pipeline's own blank-ANI = divergence signal, here reflecting true
    taxonomic mismatch: the entire S. cerevisiae reuse group is keyed to a
    C. glabrata "representative".

Consequences for this document:

1. The gene-set problem is **not** "how do we train a divergent strain" — it is
   "how do we stop a misidentified genome from being promoted to representative
   and poisoning a species' shared Trinity/ab-initio plan." The generalizable
   guards in this document (alignment floor, storeDir invalidation,
   representative-species validation) still apply; the immediate remediation is
   **re-pick the S. cerevisiae representative** (and treat KCTC's genome as its
   own independent unit / reclassify), not "build KCTC its own RNA-seq".
2. The shared Trinity is **healthy and reusable** — it was built from BY4743
   reads and aligns to S288C-family genomes; it only fails because it was aimed
   at the misidentified genome. Rebuilding it against a correctly-picked
   S. cerevisiae representative will succeed.
3. Because every sibling already has `reuse_eligible=False`, no bad shared
   ab-initio parameters were propagated (the fail-closed design held) — the
   blast radius is limited to wasted compute + a wrong representative label.

Original Hypothesis A (stale `storeDir` artifact) and Hypothesis B (genuinely
divergent strain) are downgraded/replaced by this verdict; Hypothesis C
(misidentified/mislabeled genome promoted to representative) is the resolved
root cause. The workflow hardening proposed below (Options 1-3) is still what
makes this degrade gracefully and fail loud in general, and it is what an
expert should ratify alongside the data remediation.

Companion docs:

- `nextflow/docs/HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md` — how the shared
  `rnaseq_data/<species>.trinity-GG.fasta` is single-sourced from the
  representative and fanned out to every sibling's `FUNANNOTATE_TRAIN`.
- `nextflow/docs/GENEMARK_RUN_DESIGN.md` — the established pattern this project
  uses for design docs: evidence-grounded, code-referenced, with explicit
  "resolved/open" tracking.

## Executive summary

1. `Saccharomyces_cerevisiae_KCTC_13826BP` was picked (or retained) as the
   species representative for RNA-seq purposes, but the shared Trinity its
   training reused is effectively unusable against its own genome: PASA
   assigned **36 of 5,186** transcripts to loci, funannotate train aborted
   (exit 1).
2. **Forensics (2026-08-26) resolved the cause: the representative genome
   `GCA_026225675.1` is a misidentified Nakaseomyces glabratus (Candida
   glabrata)**, not S. cerevisiae — confirmed both by NCBI's taxonomy/ANI
   check (99.93% ANI to C. glabrata; submitted-as-S. cerevisiae fails) and by
   independent alignments (trinity/reads → S288C 99.8%, → KCTC 0.3%; KCTC
   whole genome → S288C only 5.1% covered). Read provenance:
   `SRR3491422` = S. cerevisiae **BY4743** (lab strain, Exeter study). The
   shared Trinity is healthy — it was merely aimed at the wrong genome.
3. The end-state logic we want to encode (and keep as a safety net): when a
   representative's shared Trinity alignment collapses, the strain should not
   crash — it should **degrade to ab-initio-only** and be flagged for
   maintainer review (re-pick representative / reclassify / build own RNA-seq).
4. Recommended implementation package (small, fail-closed, generalizable):
   a representative-path **alignment floor** (`train_min_pasa_loci`), a
   **storeDir invalidation sentinel** keyed on representative identity, and
   **representative-species validation** (data-driven exemption for blank-own-ANI
   representatives + species/ANI cross-check against NCBI-style taxonomy). The
   "build its own RNA-seq" track (Option 4) is **not** the answer for this
   specific case — the genome is the wrong species, so the RNA-seq would be
   built for a C. glabrata isolate, not S. cerevisiae.
5. The original Hypothesis A (stale artifact) vs Hypothesis B (divergent
   strain) framing is superseded by Hypothesis C (misidentified genome); the
   expert review is now asked to ratify the data remediation (re-pick
   representative) and the workflow hardening (Options 1-3), rather than to
   adjudicate whether to spend compute on an independent KCTC RNA-seq build.

## Background and observed failure

### The failing run

- Strain: `Saccharomyces_cerevisiae_KCTC_13826BP`
  (`GCA_026225675.1_ASM2622567v1`, BIOPROJECT PRJNA891108, NCBI taxid 4932).
- Sample row in `samples.csv`:
  `GCA_026225675.1_ASM2622567v1,Saccharomyces cerevisiae,KCTC 13826BP,PRJNA891108,4932,dikarya,Ascomycota,,Saccharomycetes,,Saccharomycetales,Saccharomycetaceae,Saccharomyces,Saccharomyces cerevisiae,1,F4ABBCDF`.
- Failing task: `work/funannotate/d0/eb840120d9c9d5f9109fb06220d2c0/`
  (`FUNANNOTATE_TRAIN`, representative branch, `pasa_tier=stringent`).
- `.command.out` (2026-08-26 13:51 / 01:56 PM):
  `PASA assembled 5186 transcripts` → `PASA assigned 36 transcripts to 36 loci (genes)`
  → `[ERROR] funannotate train failed for Saccharomyces_cerevisiae_KCTC_13826BP (exit 1)`.
  `pasa_asmbls_to_training_set.dbi` never gets created → `funannotate_train.pasa.gff3`
  never produced → train aborts.

### Precedent: this failure class has already occurred once (Ascochyta_rabiei)

`utils.nf:178-189` (`loadRnaseqRepresentativeOverride` header comment) documents
a directly analogous prior case, and it is the reason
`rnaseq_representative_override_csv` exists at all: `Ascochyta_rabiei`'s
ANI+BUSCO-picked representative is a perfectly good assembly (a pks1-deletion
construct genome, BUSCO 98.1%) but the species' real RNA-seq (SRR330019xx)
barely aligns to it — HISAT2 → 9 Trinity-GG clusters → **7 transcripts total**.
`GCF_004011695.2` (Me14, BUSCO 98.4%, the RefSeq reference for the species) is
not the ANI pick but works fine as the RNA-seq anchor. The override CSV was
built specifically because **ANI+BUSCO representative selection optimizes for
assembly quality, not for which genome the species' actual RNA-seq reads align
to**, and the two occasionally disagree badly.

This is stronger, more directly on-point grounding for Hypothesis B-shaped
failures (good assembly, wrong RNA-seq anchor) than the Hansenula rescue case
cited elsewhere in this doc — Hansenula was a rescue via relaxed PASA
thresholds on a shared Trinity that still aligned at some usable rate;
Ascochyta_rabiei is a case where the ANI-picked representative's genome should
never have been the RNA-seq anchor at all, and the existing override machinery
already solves exactly that. It also means Option 3 (data-driven representative
exemption) and Option 4 (independent per-strain track) are not purely
speculative designs — they extend a pattern the codebase has already needed
once before. It is an *intra-species* sibling of the KCTC failure; for KCTC
itself the alignment forensics (below) have since resolved the question — the
anchor mismatch was a **taxonomic misidentification** (Hypothesis C), a
stronger form of "this genome should never have been the species' anchor".

### The assignment context

- In `genome_annotation/_reuse_assignments/repr_assignments.tsv`:
  `Saccharomyces_cerevisiae_KCTC_13826BP` is `is_representative=True`,
  `representative_out=Saccharomyces_cerevisiae_KCTC_13826BP`,
  `ani_to_representative=100.0`, `reuse_eligible=False`.
- All ~1,364 sibling S. cerevisiae instances in the ASCO scope have **blank**
  `ani_to_representative` and `reuse_eligible=False`.
- Per the code's own convention (`nextflow/modules/funannotate/utils.nf`,
  `loadAbinitioReuseMap`, ~lines 158-163): *a blank ANI is itself a signal of
  exceptional divergence, not "not yet computed"*. So under the pipeline's
  semantics, the **entire S. cerevisiae sibling cohort already reads as
  "divergent from representative"** — and the representative itself has no ANI
  evidence against any peer (only against itself).
- Earlier speculation ("KCTC must have been picked for lack of ANI coverage /
  default-retained") is **superseded by the rep-pick audit (2026-08-26)**: the
  picker's actual data source is `tables/busco_genome.parquet` +
  `tables/asm_stats.parquet` (ASMID join), and `GCA_026225675.1` **is** present
  there — `complete_pct=75.80` (tied for the group maximum among the 1,308
  strains) and **N50=1,117,196 bp (the group maximum)**. Faithful replication
  of `pick_representative_strain.py` against the real ANI TSV
  (`all_pairs_merged.tsv`, Aug 17) + parquet returns
  **KCTC = rank 1 of 1,308**. It also had ANI coverage — its *only*
  within-species ANI pair is **100.00% to `GCA_051107375.1` ("MRD-KRBAY")**,
  i.e. it clusters 1:1 with a second co-mislabeled genome (Cluster 2, N=2 in
  the skani ANI report). **No algorithm bug**: the picker ranked exactly as
  designed (BUSCO → N50 among ANI-covered candidates); the failure is that two
  mislabeled *C. glabrata* genomes are members of the S. cerevisiae group and
  one of them happens to hold the top assembly-quality metrics.

### The RNA-seq artifact context (cache forensics)

- Shared Trinity used:
  `rnaseq_data/Saccharomyces_cerevisiae.trinity-GG.fasta` — 5,186 transcripts,
  **mtime 2026-05-31 20:59** (~33 MB).
- There is **no** `rnaseq_data/Saccharomyces_cerevisiae.funannotate-trinity.log`
  — the expected `RNASEQ_PREPARE` sidecar. By contrast,
  `rnaseq_data/Saccharomyces_x_bayanus.funannotate-trinity.log` **does** exist.
  Absence of the log is consistent with `RNASEQ_PREPARE` having been **skipped
  via Nextflow `storeDir`** caching (artifact already present → process not
  re-run), i.e. the Trinity predates this run.
- Reads (`rnaseq_reads/Saccharomyces_cerevisiae_norm_R1/R2.fastq.gz`, ~147 MB
  each) dated **2026-05-26**, symlinked back to
  `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/rnaseq_reads/`.
- `rnaseq_reads/sra_query/Saccharomyces_cerevisiae.sra_query.csv` is **empty**
  (header only: `species_tag,taxonid,sra_accession,spots,platform,layout`) —
  no accessions were resolved/fetched for S. cerevisiae in the current round.
- So we cannot currently attribute the reads: whose strain's RNA-seq is this
  Trinity built out of, against which genome? (Open question → Hypothesis A/B.)
- Known other thin/mis-anchored artifacts under `rnaseq_data/` by species are
  already reaped into `misc/poor_trinity/` (e.g.
  `Saccharomyces_eubayanus_x_Saccharomyces_uvarum.trinity-GG.fasta` 54 KB,
  `Saccharomyces_x_bayanus.trinity-GG.fasta` 252 B) — this S. cerevisiae case is
  the same failure class, but with a **count-compliant** (5,186 ≥ 2,000)
  Trinity, so the existing count-based raps (`scripts/pick_rnaseq_representative_override.py`,
  `--threshold 2000`) do **not** catch it.

### Forensics results (2026-08-26) — the verdict

Run on HPCC (minimap2 2.30 module, samtools 1.22.1) in
`/scratch/jstajich/27800197/opencode/rnaseq_forensics/`. Genomes:
KCTC = `input_clean_genomes/GCA_026225675.1_ASM2622567v1.fa.gz`; S288C =
`input_clean_genomes/GCF_000146045.2_R64.fa.gz`; sibling control Y10 =
`input_clean_genomes/GCA_000192375.1_Saccharomyces_cerevisiae_Y10-2.0.fa.gz`;
also S. kudriavzevii / S. eubayanus references.

| Evidence | Result | Read |
|---|---|---|
| Trinity (5,186 tx) → S288C R64 | **5,175 / 5,186 (99.8%)** aligned, mapQ≥1; **4,259 ≥ 80% cov** | healthy S288C-background transcriptome |
| Trinity → KCTC genome | only 1,273 aligned; **17 (0.3%) ≥ 80% cov** | Trinity is not from/for KCTC |
| KCTC genome → S288C (`asm5`) | only **5.1%** of KCTC's 12,908,789 bp aligns; **0 contigs ≥ 75% covered** | strain-divergence would be ≥95% covered at ≥97% id; this is foreign |
| KCTC → S. kudriavzevii / S. eubayanus | 1.59% / 4.04% | not a sensu-stricto censit |
| Y10 (genuine S. cerevisiae sibling) → S288C | **100.5%** | siblings are real yeast; KCTC is the sole outlier |
| Genome stats | KCTC 12.9 Mb, GC 38.53% | — |
| NCBI taxonomy/ANI auto-check on `GCA_026225675.1` | `is_atypical: true`, `match_status: "mismatch"`, notes `unverified source organism`; **submitted as "*Saccharomyces boulardii* (nom. inval.)"** (Kangwon National Univ., bioethanol/food isolate, India 2012); best ANI = **Nakaseomyces glabratus (C. glabrata) 99.93% ANI, 96.9% coverage**; vs S. cerevisiae 87.2% ANI / 0.95 cov (fail) | independent confirmation |
| BUSCO table footprint | 0 rows in `tables/BUSCO.csv.gz` (annotation-level, irrelevant to the pick) — but **present in `tables/busco_genome.parquet`** with `complete_pct=75.80` (tied **group max**) and **N50=1,117,196 bp (group max)** → this is what made it rank-1 representative | rep DID have BUSCO evidence; it topped it |
| ANI coverage footprint | KCTC's **only** within-species ANI pair = **100.00% to `GCA_051107375.1` ("MRD-KRBAY")** (co-mislabeled *C. glabrata*; skani Cluster 2, N=2, at <90% to the other 1,305) | its "ANI coverage" was a second mislabeled genome, not its species |
| Reuse footprint | **1,308** sibling rows in `abinitio_reuse_assignments.csv` point at KCTC rep; all blank ANI, `reuse_eligible=False` | whole S. cerevisiae reuse group keyed to a foreign genome |

**Interpretation.** The Trinity is healthy and S288C-anchored; the genome it
failed to train against (KCTC) is a misidentified/C. glabrata genome promoted to
species representative. This is neither a stale-artifact bug (A) nor a genuinely
divergent-but-real yeast (B) — it is a **taxonomic misidentification** (C). The
pipeline's fail-closed reuse semantics (blank ANI ⇒ not eligible) held, so no
wrong ab-initio params propagated; the blast radius is one mislabeled
representative + wasted representative-branch compute + a shared Trinity that
must be rebuilt against a correctly-picked representative.

## Why the pipeline currently cannot handle this (code-grounded)

### Gate inventory

1. **`train_min_trinity_transcripts`** (`nextflow/conf/profile_funannotate.config`
   ~line 191, default 2000; set to 0 to disable) — enforced inside
   `FUNANNOTATE_TRAIN/main.nf`. Counts FASTA records in the Trinity input.
   Guards against a too-thin assembly, but **cannot distinguish** a healthy
   5k transcript assembly that aligns to the genome from one built out of
   different reads/against a different genome. This gate passed on 5,186.
2. **`pasa_tier`** (`FUNANNOTATE_RNASEQ.nf`, `pasaTierFor(...)`, ~lines 272-279 —
   **not** `utils.nf`; corrected from an earlier draft of this doc):
   - `ani_to_representative >= 97.0` → `'stringent'` (PASA defaults; includes
     the representative itself, whose ANI is hardcoded to 100.0 in
     `loadAbinitioReuseMap`).
   - `90.0 <= ani_to_representative < 97.0` → `'relaxed'`
     (`--pasa_min_avg_per_id 85 --pasa_min_pct_aligned 70 --pasa_num_bp_splice 1`,
     `profile_funannotate.config` ~lines 161-163) — the cross-strain rescue
     tier that saved e.g. a Hansenula case with 19/14,319 transcripts.
   - `ani_to_representative < 90.0` **or blank/null** → `'skip'` (no shared
     Trinity; ab-initio-only), a deliberate fail-closed degradation via the
     `predict_no_rnaseq` branch (`FUNANNOTATE_RNASEQ.nf:249`).
   - Note: `ani_reuse_threshold` (default 99.0, `profile_funannotate.config:215`)
     is a **different, unrelated** parameter — it only feeds `reuse_eligible`
     in `PICK_REPRESENTATIVE_STRAIN` and `BACKFILL_ABINITIO_PARAMS` (ab-initio
     parameter reuse), not the PASA tier cutoffs above. An earlier draft of
     this doc conflated the two; any expert-facing brief should use the
     97.0/90.0 cutoffs above, not `ani_reuse_threshold`.
   **Gap:** the representative itself never lands in the `'skip'` branch —
   it always goes `'stringent'` with whatever shared Trinity exists. There is no
   representative-path equivalent of "this evidence is unusable → degrade".
3. **Empty-reads sentinel** (`RNASEQ_PREPARE/main.nf`): no reads → empty
   Trinity → empty `shared_fa` → FUNANNOTATE_TRAIN skips cleanly (**exit 0**,
   ab-initio-only). This is the only representative-path graceful degradation
   that exists, and it only triggers when there are *no* reads at all — not when
   reads/Trinity exist but don't align.
4. **`storeDir "${launchDir}/rnaseq_data"`** (`RNASEQ_PREPARE/main.nf:8`):
   the shared Trinity is cached by **file presence/path**, keyed on `species_tag`
   only. If the representative strain (or the reference used to build the
   Trinity) changes but the species tag doesn't, `storeDir` serves the stale
   artifact and `RNASEQ_PREPARE` never re-runs. **This is the mechanism that
   can serve a mis-anchored Trinity and there is no invalidation.**
5. **`rnaseq_representative_override_csv`** (`profile_funannotate.config`
   ~line 200, default `"${launchDir}/rnaseq_representative_override.csv"`) +
   `loadRnaseqRepresentativeOverride`/`saveRnaseqRepresentativeOverride`
   (`utils.nf`) + `scripts/pick_rnaseq_representative_override.py`: lets a
   maintainer force a *different strain* to build the species' shared Trinity.
   Only fires on the low-transcript-count signal. Current file contains one
   unrelated hybrid-species row, no S. cerevisiae override.
6. **ANI aggregate `storeDir` — FIXED 2026-08-26** (see Option 2 extension
   below for the implementation). Originally: (`nextflow/modules/ani/report/{MERGE_ANI,REPORT_ANI}/main.nf`)
   — the same failure class as gate 4, one layer upstream: both modules use
   `storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}"`,
   keyed on the ANI group's *name* (species or genus tag), not on the group's
   *membership* (which genomes are actually in it). Confirmed on this dataset:
   `results/ANI/skani/SPECIES/Saccharomyces_cerevisiae/`,
   `results/ANI/skani/GENUS/Saccharomyces/`,
   `results/ANI/skani/SPECIES/Nakaseomyces_glabratus/`, and
   `results/ANI/skani/GENUS/Nakaseomyces/` all already exist on disk. Once
   `samples.csv` is corrected (KCTC + MRD-KRBAY move from `Saccharomyces` to
   `Nakaseomyces`), re-running `compare_ani` will **not** regenerate any of
   these four directories — `storeDir` only checks output-path existence, so
   the S. cerevisiae report would keep reflecting the old (foreign-genome-
   inflated) membership and the Nakaseomyces report/genus dir would never gain
   the two corrected genomes, silently. This is the same bug as gate 4
   (`RNASEQ_PREPARE`'s `storeDir`) one step upstream in the pipeline, and it
   would have fed `PICK_REPRESENTATIVE_STRAIN` stale ANI tables even after
   `samples.csv` and `RNASEQ_PREPARE`'s cache were both correctly invalidated.
   The pairwise `*_COMPARE` modules (`SKANI_COMPARE`, `MASH_COMPARE`,
   `SOURMASH_COMPARE`, `FASTANI_COMPARE`) already avoid this — their code
   comments explicitly note they use `publishDir` + normal input-hash caching
   *because* `storeDir` "only checks output-path existence... regardless of
   `-resume`"; only the group-level `MERGE_ANI`/`REPORT_ANI` aggregation step
   was left on `storeDir`. Needs the same invalidation-sentinel treatment as
   Option 2 below (see Option 2 extension).

## The logic we want to encode

User-stated policy (paraphrase): *if the representative strain is too far from
the transcriptome reference to train on, then either (a) build that strain its
own RNA-seq, or (b) skip RNA-seq-based training for it and its siblings
(ab-initio-only).*

In pipeline terms that means three behaviors, in order of preference:

1. **Detect** the unusable-evidence state — a count-compliant Trinity whose
   alignment to the representative genome collapses below a floor (or a blank
   own-ANI representative) — and **record** it (report TSV, not a crash).
2. **Degrade gracefully**: the affected strain(s) flow to
   `predict_no_rnaseq` (ab-initio-only), the same path an empty-reads genome
   already takes today. Representative + siblings consistent (no mixed
   half-trained family).
3. **Offer the upgrade path**: an explicit opt-in per-strain "build its own
   RNA-seq" track (SRA fetch + Trinity + train keyed to that strain alone) for
   when a maintainer decides the strain matters enough to spend the compute.

## Strategy options evaluated

Ordered roughly by scope/risk. Recommendation at the end.

### Option 0 — one-off manual unblock (EXECUTED 2026-08-26)

- **DONE** — stale S. cerevisiae Trinity invalidated: moved
  `rnaseq_data/Saccharomyces_cerevisiae.trinity-GG.fasta` →
  `misc/quarantine/Saccharomyces_cerevisiae.trinity-GG.fasta_mislabeled_rep_2026-08-26`
  (preserved as evidence). Rebuilt on next `RNASEQ_PREPARE` once the
  representative is valid.
- **DONE** — read provenance: SRR3491422 = BY4743 (S288C background).
- **DONE** — forensics (verdict: Hypothesis C / misidentified genome; requires
  representative **re-pick**, not an independent KCTC RNA-seq build).
- **PENDING** — actual re-pick of S. cerevisiae representative + rerun of the
  S. cerevisiae annotation group (data change to reuse-assignment files;
  awaiting user decision).
- Scope/risk: data fix, no code. Does not prevent recurrence (that is Options
  1-3).

### Option 1 — representative-path alignment floor (recommended core)

- New param `train_min_pasa_loci` (e.g. default 100) in `FUNANNOTATE_TRAIN`.
- After PASA runs, parse the `PASA assigned N transcripts to N loci (genes)`
  value from the PASA/train output (need a reliable source — see
  Software/Workflow question Q-SW-2). If `N < floor`:
  - **exit 0** (graceful, matching the existing empty-reads skip path), and
  - write a record to e.g. `misc/rnaseq_training_failed_low_alignment/<out>.tsv`
    (or a single append-style TSV) so the strain's downstream predict routes
    via `predict_no_rnaseq` (ab-initio-only), and the report is visible.
- Also consider running the same floor check inside `RNASEQ_PREPARE` at build
  time (align fresh Trinity-GG back to the representative genome; if below
  floor, don't archive to `rnaseq_data/`, drop into `misc/poor_trinity/`
  pattern, and record "species needs rnaseq representative review").
- Calibration: use the Hansenula precedent (relaxed `85/70/1` rescued
  19/14,319) as the upper bound of "worth one rescue attempt"; anything below
  that is ab-initio-regime. Expert-validated floor asked in Bio question Q-BIO-3.
- Scope/risk: small module change; the crux is a *reliable* way to read PASA's
  assigned-loci count without fragile string-parsing.
- Verdict: **recommended**.

### Option 2 — storeDir invalidation keyed on group membership/identity (recommended core)

- `RNASEQ_PREPARE` records `rnaseq_data/.built_by/<species_tag>.tsv` holding
  (species_tag, representative `out`, genome fingerprint/ASMID, date).
- `RNASEQ_PREPARE` (or a preflight in `FUNANNOTATE_RNASEQ.nf`) compares the
  current representative's identity to the recorded one; on mismatch, move the
  stale artifact aside (preserve for forensics) and force a rebuild, bypassing
  `storeDir` cache for that species.
- Scope/risk: small; must respect Nextflow `-resume`/cache semantics (see
  Q-SW-3). Prevents the "rep changed, cache key unchanged" class of bug that
  produced the May 31 → Aug 26 stale serve.
- **Extension, `MERGE_ANI`/`REPORT_ANI` — IMPLEMENTED 2026-08-26**: same
  failure class as `RNASEQ_PREPARE`'s `storeDir`, one layer upstream (gate
  inventory item 6). The originally-sketched fix here was a membership-
  fingerprint file (`${group_name}.built_from.tsv` listing the ASMID set that
  produced the stored report, compared on each run) — but the pairwise
  `*_COMPARE` modules (`SKANI_COMPARE`, `MASH_COMPARE`, `SOURMASH_COMPARE`,
  `FASTANI_COMPARE`) had already solved this exact problem more simply:
  **drop `storeDir` entirely, keep only `publishDir`**, and let Nextflow's own
  input-hash caching do the invalidation (it already hashes `part_tsvs`/
  `ani_tsv`+`names_tsv`, which change whenever group membership changes — no
  custom fingerprint file needed). Applied the identical fix to `MERGE_ANI`
  (`nextflow/modules/ani/report/MERGE_ANI/main.nf` — was `storeDir`-only, now
  `publishDir`) and `REPORT_ANI` (`.../REPORT_ANI/main.nf` — had **both**
  `storeDir` and `publishDir` pointed at the same path; `storeDir`'s bare
  existence check silently won regardless, so the `publishDir` was
  vestigial). Verified with a stub-run before/after: pre-fix, a second
  invocation logged `[skipping] Stored process > COMPARE_ANI:REPORT_ANI`
  (blind skip); post-fix, the same invocation shows
  `[PROCESS ...] COMPARE_ANI:REPORT_ANI (Testus_fungus)` — it actually
  re-executes and picks up membership changes. **Trade-off, stated plainly**:
  `publishDir` + normal caching only skips recomputation within the *same*
  `work/` directory (via `-resume`); unlike `storeDir`, it provides no
  skip-if-unchanged behavior across a `work/` wipe between runs, so a fully
  fresh run now always recomputes every group's merge/report rather than
  reusing untouched ones. This is the same trade-off already accepted for the
  pairwise `*_COMPARE` modules; consistency with that established precedent
  was preferred over reintroducing a second stale-cache mechanism.
  `RNASEQ_PREPARE`'s own representative-identity sentinel (the `storeDir`
  call this Option was originally about) is **still design-only, not yet
  implemented** — it needs the fingerprint-file approach (or an equivalent
  `-resume`-safe check) since the Trinity-build output can't as easily switch
  to input-hash-only caching without recomputing large per-species
  transcriptome assemblies on every run.
- Verdict: **ANI side (`MERGE_ANI`/`REPORT_ANI`) done; `RNASEQ_PREPARE` side
  still recommended, not yet implemented.**

### Option 3 — data-driven representative exemption (recommended core)

- Extend `pasaTierFor`/add `rnaseqPolicyFor(out)` in `utils.nf` so a
  **representative** whose own condition is:
  - blank ANI row (the codebase's own divergence signal), OR
  - prior low-alignment record (from Option 1's report),
  is routed to the `'skip'`/`predict_no_rnaseq` branch — the same fail-closed
  degradation already granted to blank-ANI siblings today.

  Concretely the representative would (a) not be allowed to *build/share* a
  Trinity in this round (skip `RNASEQ_PREPARE` staging or produce the
  empty-reads sentinel), and (b) predict ab-initio-only. This directly encodes
  "skip rnaseq-based training for it" as an automatic, data-driven rule rather
  than a manual override.
- Fail-closed caveats: blank-ANI-on-everyone can be a skani/sourmash/coverage
  data gap, not biology — so this rule needs the same escape hatch
  `rnaseq_representative_override_csv` already provides (a maintainer can always
  force a different anchor strain back in). See Q-BIO-2 and Q-SW-5.
- Scope/risk: medium; central logic in `utils.nf`, must stay fail-closed and
  override-able.
- Verdict: **recommended** (as the automatic engine behind the policy; Option 1
  supplies the alignment evidence it consumes).

### Option 4 — opt-in independent RNA-seq track (gated, higher cost)

- New per-strain mode, e.g. `--rnaseq_independent_override_csv` (or extend the
  existing override CSV with a column/flag) → for listed `out`s, run a
  dedicated SRA query + fetch + Trinity + `FUNANNOTATE_TRAIN` with the strain's
  **own** reads against its **own** genome, bypassing the shared-assembly path.
- Default OFF; used only when the bioinformatics expert deems the strain
  biologically important enough and public reads exist (Q-BIO-4).
- Must reuse the existing per-species fetch machinery (`SRA_QUERY_BATCH`,
  `SRA_FETCH`) with a per-strain key instead of the species key.
- Scope/risk: moderate; new branch in `FUNANNOTATE_RNASEQ.nf` + param; cost
  gated by usage.
- Verdict: **optional, expert-gated**.

### Option 5 — decouple rnaseq representative from ab-initio representative (follow-up)

- The ab-initio/shared-param representative (assembly-quality pick: BUSCO, N50)
  and the RNA-seq anchor (must align a usable transcriptome) are conceptually
  different roles. The codebase already gestures at this (
  `rnaseq_representative_override`, `pick_rnaseq_representative_override.py`)
  but the selection math in `pick_representative_strain.py` /
  `species_reuse_clusters.py` still picks one strain for both.
- Longer-term: a strain that cannot serve as rnaseq anchor should lose rnaseq
  anchor duties (demote to "ab-initio-only leader"), possibly while retaining
  ab-initio reuse leadership for its close siblings.
- Scope/risk: deeper refactor of the representative-pick; needs the alignment
  evidence from Option 1 to be meaningful. Defer.
- Verdict: **follow-up**.

### Option 7 — proactive dataset-wide misidentification sweep (new, cluster-vs-label cross-check)

Everything above (Options 1-3, 6) is *reactive*: KCTC only surfaced because it
happened to get picked as a representative and crashed training. A genome that
is mislabeled but never gets picked (e.g. a low-BUSCO/low-N50 member of its
false species group) would sit quietly forever, silently poisoning
`reuse_eligible`/blank-ANI assignments for its (real) siblings, with no
tripwire at all. Researched 2026-08-26 (background investigation, no code
changes) against the actual `nextflow/modules/ani/` codebase; findings below.

**Two relevant mechanisms already exist in this codebase, wired for a
different trigger than "verify every label":**

- `MASH_PREFILTER` → `MASH_COMPONENTS` (connected-components clustering at
  `--prefilter_ani`, default 80%) → confirmatory comparison only *within* each
  component — a real divide-and-conquer cascade in `ANI_COMPARE_METHOD.nf`,
  but gated behind `ani_method == 'fastani' && params.fastani_prefilter`.
  Production runs `ani_method = 'skani'` (whose native `triangle` mode is
  already fast enough at this scale that the cascade was never turned on), so
  this exists but is unused — and even if enabled, it only prunes redundant
  comparisons *inside* an already-claimed group; it never questions the
  group's own label.
- `query_ANI.nf` + `SKANI_DIST_QUERY` already do asymmetric query-vs-reference
  placement (cost O(queries×references), not O(n²)), with a `--fallback` that
  walks FAMILY→ORDER→CLASS on a miss — exactly the right cost shape for a
  screening pass. But its query set is defined as **rows with a blank**
  taxonomy rank (orphans). KCTC/MRD-KRBAY had a fully-populated, confidently
  *wrong* `Saccharomyces cerevisiae` label, not a blank one, so `query_ANI`
  as written would never select them.
- No code anywhere cross-checks a species-level ANI cluster's majority label
  against its members' claimed `SPECIES` value. `species_reuse_clusters.py`
  and `pick_representative_strain.py` only ever cluster *within* a
  pre-grouped species tag — the group membership itself is trusted, never
  audited. This is the actual gap, and it is the cheapest one to close: the
  KCTC/MRD-KRBAY case was already visible as "Cluster 2, N=2, <90% to the
  1,305-genome majority cluster" inside data the pipeline had already
  computed (`results/ANI/skani/SPECIES/Saccharomyces_cerevisiae/batches/
  Saccharomyces_cerevisiae.full.ani.tsv`) — nothing needed to *look* at it.

**(a) — cheap, do now (no new Nextflow compute):** a post-hoc, off-Nextflow
script that reads the already-`storeDir`'d per-group ANI report/pair tables
under `results/ANI/skani/{SPECIES,GENUS}/`, runs single-linkage/connected-
components at the existing `cluster` threshold (95%, `profile_ANI.config`),
and flags any cluster whose member ASMIDs' claimed `SPECIES` labels are not
unanimous — reporting the minority-label members as misidentification
candidates (this is precisely the KCTC-class signal). Depends on Option 2's
storeDir-invalidation fix being in place first (or being run right after a
fresh `compare_ani`), so it isn't scanning stale group data. Zero new compute;
this alone would have caught KCTC+MRD-KRBAY dataset-wide without needing a
training crash as the tripwire.

**(b) — later, if (a) proves useful:** extend `query_ANI.nf` with a
`--verify_all` mode (or a query-selection sentinel beyond "blank rank") that
treats every genome in a compare-group as a query against a small per-cluster
centroid/reference set, turning the existing asymmetric machinery into a
general label-verification pass at O(N×centroids) instead of O(N²) — useful
once genus-level all-vs-all becomes too expensive to run comprehensively at
full dataset scale.

Phylogenetic placement (EPA-ng/pplacer onto a reference tree) was also
considered — no existing module does this (`nf_phyling` builds a BUSCO-marker
tree for a *fixed* input set, it doesn't place new queries), and ANI already
discriminates genus-level swaps at near-zero marginal cost once sketches
exist, so placement isn't the lowest-cost option for this specific failure
mode.

- Scope/risk: (a) is a standalone Python script over existing output files,
  no Nextflow/module changes, no risk to the running pipeline — safe to
  implement immediately. (b) is a real module extension to `query_ANI.nf`,
  deferred.
- Verdict: **(a) recommended now; (b) follow-up, gated on (a)'s findings.**

**(a) — IMPLEMENTED 2026-08-26**: `scripts/find_ani_label_mismatches.py`.
Reads `results/ANI/skani/SPECIES/all_pairs_merged.tsv` (already-computed,
whole-dataset merged pairwise ANI, 2,195,284 pairs / 1,770 species groups) +
`samples.csv`; for each claimed-species group, union-finds edges at
`--cluster-threshold` (default 95.0, matches `ani_cluster_threshold`) and
flags any group that splits into more than one connected component. Each
non-majority component is reported with the max ANI from any of its members
to the majority cluster, tiered `HIGH` confidence (max ANI < `--outlier-
threshold`, default 90.0, matching the already-documented
`ani_outlier_threshold` / README_ANI.md convention) vs `REVIEW` (91-95%,
plausibly just real divergent conspecific strains or assembly-quality noise,
not necessarily misidentification). Runs the whole 23k-genome dataset in
~10 seconds, zero new Nextflow compute.

First run (2026-08-26) found **170 HIGH-confidence flagged clusters across
~150 species groups** (410 total incl. REVIEW-tier, across 219 species).
Confirms the method: `Saccharomyces cerevisiae` → exactly
`GCA_026225675.1_ASM2622567v1,GCA_051107375.1_ASM5110737v1` at 0.00% ANI to
the 1,305-genome majority (the KCTC/MRD-KRBAY case this whole doc is about) —
reproduced automatically, no training crash needed as a tripwire. Also
surfaces new, previously-unreviewed HIGH-confidence candidates dataset-wide,
e.g. `Candidozyma auris` (2 separate 0.00%-ANI singleton outliers),
`Macrophomina phaseolina` (a 2-genome 0.00%-ANI pair), and — notably —
`Nakaseomyces glabratus` itself already has its own internal 2-genome
0.00%-ANI outlier pair, independent of and before KCTC/MRD-KRBAY are added to
that group by the Option 0/3a remediation. None of these new candidates have
been forensically confirmed the way KCTC was (NCBI atypical-check +
independent minimap2 alignment) — this script produces a triage list, not a
verdict; `HIGH`-confidence flags (nowhere near the cluster threshold, often
exactly 0.00% = no measurable ANI at all) warrant the same NCBI-taxonomy-
check + alignment forensics workflow used for KCTC, while `REVIEW`-tier
flags (91-95%) are more likely genuine divergent strains or fragmented/
low-quality assemblies than mislabeling. Output not yet triaged/actioned
beyond the KCTC/MRD-KRBAY case already in hand — follow-up task, not blocking
this doc's core remediation.

### Option 8 — skani divide-and-conquer prefilter (new, IMPLEMENTED 2026-08-26)

Fork research into Option 7 (see above) surfaced that a real
mash-prefilter → connected-components → confirmatory-comparison-within-
component cascade already existed in `ANI_COMPARE_METHOD.nf`, but was gated
behind `method == 'fastani'` — unreachable in production, which runs
`ani_method = 'skani'`. At `SPECIES` scope (~1,300 genomes/group max seen
here) a single `skani triangle` job over the whole group is already fast
enough that this never mattered; at `GENUS` scope, where a group can be
substantially larger, running the full pairwise comparison over genomes that
are obviously unrelated (e.g. via a cheap mash prefilter) is wasted compute —
and it's the same divide-and-conquer shape Option 7 already needs at genus
scale for the misidentification sweep to stay cheap.

**Implemented:**
- `nextflow/conf/profile_ANI.config` — new `skani_prefilter = false` param,
  independent of `fastani_prefilter` (so enabling one doesn't change the
  other's default), reusing the existing `prefilter_ani`/`min_group_size`
  params.
- `nextflow/modules/ani/compare/SKANI_COMPARE/main.nf` — added a `batch_tag`
  to its input/output signature (mirroring `FASTANI_COMPARE`), so multiple
  per-component invocations for the same `group_name` don't collide on the
  previously-hardcoded `${group_name}.full.ani.tsv` filename. Whole-group
  calls (prefilter off) pass `"full"`, preserving today's filenames/behavior
  exactly.
- `nextflow/subworkflows/local/ANI_COMPARE_METHOD.nf` — extracted the
  fastani branch's mash-prefilter/connected-components glue into a new
  `PREFILTER_COMPONENTS` subworkflow (a **subworkflow**, not a Groovy
  closure — Nextflow does not allow process calls inside closures, confirmed
  by a real compile error during implementation: `Processes cannot be called
  from within a closure`). Both the `fastani_prefilter` and new
  `skani_prefilter` branches now call this same subworkflow instead of
  duplicating the channel-wrangling. When `skani_prefilter=true`, per-component
  `SKANI_COMPARE` outputs are recombined via the existing `MERGE_ANI`
  (same shape fastani's prefilter/batch outputs already use).
- Validated via `-stub-run` (no real skani/mash binaries needed): confirmed
  (a) the default path (`skani_prefilter=false`) still produces
  `SKANI_COMPARE (group [full] n=N)` unchanged — no regression; (b) the new
  path (`skani_prefilter=true`) correctly fires
  `PREFILTER_COMPONENTS:MASH_SKETCH` ×N → `MASH_PREFILTER` → `MASH_COMPONENTS`
  with no channel/type errors. Full behavioral validation (real skani/mash
  binaries producing real components) not done in this environment — stub-run
  only confirms DSL2 wiring correctness, not numerical output.
- Scope/risk: opt-in (`skani_prefilter` defaults to `false`), so no change to
  current production behavior unless explicitly enabled.

### Recommendation summary

Implement **Option 1 + 2 + 3** as the core package: automatic detection with
graceful ab-initio-only degradation + report (1), cache-invalidation so stale
artifacts can't silently serve (2), and the data-driven representative exemption
that makes "too far from reference → skip rnaseq training" automatic and
fail-closed (3). Run **Option 0** now to unblock and to seed the evidence
package. Keep **Option 4** as an opt-in flag, gated on expert ruling. Track
**Option 5** as a follow-up refactor. **Option 7a** (cluster-vs-label sweep) is
a standalone, no-pipeline-risk script — implement it now, in parallel with the
Option 0 remediation, since it's the only proactive (dataset-wide) guard among
all the options and costs no new compute. **Option 7b** is a follow-up gated
on 7a's findings.

All of the above is contingent on the expert adjudication in the next section;
with the forensics verdict (Hypothesis C: KCTC is a misidentified C. glabrata
genome), the weighting changes from the original drafting:

- **Data remediation (the actual unblock, not optional):** re-pick the
  S. cerevisiae RNA-seq + ab-initio representative to a genuine S. cerevisiae
  (the ~1,308 siblings are now all ANI-comparable, e.g. S288C/Y10-class
  strains); the existing healthy Trinity (BY4743 reads) will rebuild cleanly
  against a real representative. KCTC's genome should be annotated as its own
  independent unit (and flagged for possible reclassification as N. glabratus),
  never as the species' shared reference. This makes Option 4 (independent
  KCTC RNA-seq) moot — it would only be appropriate if the misidentified genome
  were itself the object of study.
- Option 2 (storeDir invalidation) + Option 0 have already proven their point:
  a stale 5,186-transcript, count-compliant Trinity served by an unchanged
  species-key cache would have hidden this silently if `FUNANNOTATE_TRAIN` had
  not crashed. Option 1 (alignment floor) is the generalizable safety net that
  turns a future crash into a graceful ab-initio-only degrade + report.
- Newly-motivated generalizable guard beyond the original three: a
  **representative-species sanity check** so a misidentified genome cannot be
  promoted to representative and key an entire species' reuse group. The
  rep-pick audit pins down the concretest, cheapest rule: a representative must
  belong to its species' **dominant ANI cluster** — e.g. an ANI ≥ `cluster`
  threshold (95%) to some minimum fraction of its species peers, or membership
  in the largest ANI cluster / same-cluster as the majority. This specific
  check would have excluded KCTC outright: its only pair was a 100% self-cluster
  with MRD-KRBAY at <90% to the other 1,305, so it fails "representative must
  sit in the species' majority cluster". Fold into Option 3 as
  `rnaseqPolicyFor` + representative-pick hardening (as "Option 6"). A
  blank-own-ANI representative and the alignment floor (Option 1) remain the
  runtime safety nets for cases that slip past pick time.

## Expert evaluation plan

### Goal

**Superseded by the forensics verdict above** (Hypothesis C: misidentified
genome, resolved 2026-08-26) — the original A-vs-B adjudication is no longer
the open question. The evidence package and Q-BIO/Q-SW items below are kept
for reference (they document the diagnostic process and most of Q-SW still
applies to the general-purpose guards), but the expert session's actual scope
is now:
1. Ratify the data remediation: re-pick the S. cerevisiae representative
   (excluding `GCA_026225675.1`), reclassify/flag it as its own unit (likely
   *Candida glabrata* / *Nakaseomyces glabratus*, not *S. cerevisiae*), and
   confirm the rebuilt shared Trinity (BY4743-derived, already shown to align
   99.8% to S288C) trains cleanly against the new representative.
2. Ratify the generalizable workflow hardening (Options 1-3), now reframed as
   **representative-species validation** guards (catch a misidentified/
   mismatched-species representative before it poisons a shared Trinity) rather
   than purely an ANI-divergence-threshold problem.
3. Decide whether Option 4 (independent per-strain RNA-seq) is still worth
   keeping as a generic opt-in for *other*, genuinely-divergent future cases —
   it is **not** applicable to KCTC itself (building RNA-seq for a
   mislabeled-species genome doesn't fix anything).

### Evidence package (assemble before the session)

1. **Failure artifact**: `work/funannotate/d0/eb840120d9c9d5f9109fb06220d2c0/`
   `.command.out`/`.command.err` with the 36/5,186 line, `pasa_tier=stringent`.
2. **Assignment slice**: S. cerevisiae rows of
   `genome_annotation/_reuse_assignments/repr_assignments.tsv` (one rep, all
   siblings blank ANI / not eligible).
3. **Cache forensics**: the mtime asymmetry (Trinity May 31, reads May 26,
   run Aug 26), the missing `Saccharomyces_cerevisiae.funannotate-trinity.log`
   vs present `Saccharomyces_x_bayanus.funannotate-trinity.log`, empty
   `sra_query` CSV.
4. **Alignment forensics**: see "Forensics results" table above. Outcome was the
   `high↔canonical, low↔KCTC` row (Trinity is S288C-background; KCTC genome is
   foreign at 5.1% whole-genome coverage, 0.3% transcript coverage) — which the
   original matrix mapped to Hypothesis A, but the NCBI taxonomy/ANI check turns
   it into **Hypothesis C (misidentification)** as documented.
5. **Read attribution**: `Saccharomyces_cerevisiae_norm_R1/R2.fastq.gz` head
   read `@SRR3491422.7/1` → **SRR3491422 = S. cerevisiae BY4743 (S288C
   background)**, glycerol-evolved (dup2) lab strain, Exeter 2016 study. The
   reads and therefore the Trinity are genuine S. cerevisiae — they were simply
   aimed at the wrong (misidentified C. glabrata) genome. R1/R2 both ~147 MB.
   The SE pair is 0 bytes; `sra_query` CSV empty (headers only).
 6. **Rep-pick audit — RESOLVED (2026-08-26)**: KCTC became representative by
    the picker's *intended* ranking — highest `complete_pct` (75.80, tied group
    max) then highest N50 (1,117,196 bp, group max) among all 1,308 ANI-covered
    S. cerevisiae strains (verified by faithful replication of
    `pick_representative_strain.py` against `all_pairs_merged.tsv` +
    `busco_genome/asm_stats.parquet`: KCTC = rank 1/1,308). Its ANI coverage
    was a **single 100% pair to co-mislabeled `GCA_051107375.1` (MRD-KRBAY)**.
    Absence from `tables/BUSCO.csv.gz` is irrelevant (that's the
    annotation-level table; the picker reads `busco_genome.parquet`). No code
    bug — the mislabeled-genome data is the failure.

**Remediation representative candidate — VERIFIED (2026-08-26).** Re-running
the picker's ranking with `GCA_026225675.1` + `GCA_051107375.1` removed (1,306
remain, all ANI-covered) yields a genuine S. cerevisiae as new rep:
- **Picker-native winner: `GCA_003277715.1_ASM327771v1` (SX6, PRJNA396809,
  NCBI 4932)** — complete_pct 75.80 (group max), N50 994,905 bp (top); pairs
  all 1,305 in-group peers, ≥95% ANI to 1,302 (in the majority cluster). Caveat:
  at the `ani_reuse_threshold` 99.0 edge — only **15/1,305** peers ≥99%, so most
  siblings fail-closed to independent prediction (correct, but little reuse).
- Reuse-center alternative: `GCA_947344615.1_NCYC478` (cp 74.40, N50 911,387;
  ≥99% to ~1,270 peers) or `GCA_003273825.1_ASM327382v1` (cp 75.50, N50
  901,330; ≥99% to ~1,269) — maximal reuse-eligibility at near-top assembly.
- `GCA_024196135.1` (the 96.3-complete_pct outlier) is a *Pleurotus*-S.
  cerevisiae fusant, tagged as its own species — correctly excluded from-group.
- Either candidate is S288C-family and aligns the existing BY4743-derived
  Trinity cleanly, so `RNASEQ_PREPARE` + `FUNANNOTATE_TRAIN` are expected to
  succeed once the taxonomy is corrected. Choosing SX6 (picker-native) vs a
  reuse-center strain is a remediation-policy decision; see Next steps item 3.
7. **Existing rescue history**: the Hansenula case (19/14,319 → relaxed
   `85/70/1` rescued) for threshold calibration; list of prior
   `rnaseq_representative_override.csv`/`misc/poor_trinity/` uses.

### Bioinformatics expert (Q-BIO-1/2 resolved by forensics; ratifies remediation + thresholds)

Q-BIO-1 and the "divergent strain" framing of Q-BIO-4/5 are answered by the
forensics verdict: KCTC is not a divergent S. cerevisiae strain, it is a
different species. Kept below for the parts still open (threshold calibration,
all-blank-ANI as a systemic signal, general-case guidance).

- Q-BIO-1: Is KCTC_13826BP genuinely divergent from mainstream S. cerevisiae
  (ANI to S288C/R64 and to the ~1,364 in-scope siblings)? Is a bioethanol
  strain's transcriptome expected to diverge substantially in the genes that
  matter for PASA's 3% assignment?
- Q-BIO-2: Is **all-blank sibling ANI** in this species a skani/sourmash
  failure or a real divergence signal? (Consequences for Option 3's
  blank-ANI triggers.)
- Q-BIO-3: Is 36/5,186 salvageable at all with relaxed PASA (cf. Hansenula
  at 19/14,319)? What is a defensible `train_min_pasa_loci` floor, biologically
  speaking? (Sanity: a 12 Mb yeast genome has ~6k genes; what fraction of loci
  must be hit before transcript training beats ab-initio-only?)
- Q-BIO-4: For Option 4 — is there public RNA-seq for KCTC-class strains (or
  the specific isolate) worth fetching? If not, is ab-initio-only an acceptable
  quality outcome for this species' annotation (GeneMark ES + Augustus/SNAP +
  EVM, no hints)?
- Q-BIO-5: Is the missing BUSCO row for this isolate an assembly/QC red flag
  (contamination? incomplete?) that should itself disqualify it as any kind of
  representative?

### Software/workflow expert (ratifies placement + behavior)

- Q-SW-1: Where exactly to enforce the alignment floor so failure *degrades*
  (exit 0, `predict_no_rnaseq`, report) rather than crashing — inside
  `FUNANNOTATE_TRAIN` after PASA, or earlier (a minimap2-based pre-alignment in
  `RNASEQ_PREPARE`)? Interaction with `-resume`/cache so a fixed strain doesn't
  re-crash.
- Q-SW-2: The most reliable, non-fragile way to capture PASA's "assigned N
  transcripts to N loci" (parse PASA log file vs stdout vs the produced
  `*.transcripts.gff3` line count). Is `funannotate_train.pasa.gff3` only
  produced on success (so we must act on PASA's own log in the failure case)?
- Q-SW-2b: `pasa_asmbls_to_training_set.dbi` is created only after
  alignment; confirm the earliest reliable failure signal.
- Q-SW-3: `storeDir` invalidation design: `rnaseq_data/.built_by/<species>.tsv`
  sentinel + forcing rebuild on representative mismatch — does this survive
  `-resume` and multi-species batches? Does it need a `.nextflow` cache nudge?
- Q-SW-4: Where does the representative-exemption (`rnaseqPolicyFor`) live —
  alongside `pasaTierFor` in `FUNANNOTATE_RNASEQ.nf` (~272-279), not `utils.nf`
  (corrected; see gate-inventory item 2 above) — so it stays fail-closed and
  ordered correctly vs `BACKFILL_ABINITIO_PARAMS`
  (GeneMark ES→ET / ab-initio store must not break when a rep goes
  ab-initio-only)? Confirm no downstream consumer assumes
  `transcript.alignments.bam` exists for a rep (precedent:
  `trainingTranscriptBamFor` returns `''` on no-reads — safe).
- Q-SW-5: Escape-hatch consistency: `rnaseq_representative_override_csv` must
  be able to force-rescind an automatic exemption (both grounding the data and
  forcing independent training for a rep). Confirm the existing override
  machinery handles that round-trip or what new columns are needed.
- Q-SW-6: Report-surface design: where the `misc/rnaseq_training_failed_low_alignment/`
  + `predict_blocked_awaiting_representative.tsv`-style artifacts land, and how
  the `--allow_independent_fallback`-family flags interact.

### Session format and outputs

- 45-60 min per expert (or one joint review). Walk the evidence package first,
  then the options table, then the recommended package (Options 1+2+3, Option 0
  done, 4 opt-in, 5 follow-up).
- Decision matrix to complete:
  - Per-case (KCTC only, via override) vs generalizable rule (Option 3).
  - Cost: independent RNA-seq build (Option 4) vs annotation-quality cost of
    ab-initio-only.
- Deliverables (aligned with repo precedent in `GENEMARK_RUN_DESIGN.md`):
  1. An ADR-style design doc revision of this plan (status → IMPLEMENTED with
     the decisions recorded).
  2. Final param list: `train_min_pasa_loci`, and any new flags
     (`--rnaseq_independent_override_csv` etc.), added to
     `nextflow_schema.json` + `profile_funannotate.config` with defaults.
  3. Changeset list per module
     (`FUNANNOTATE_TRAIN`, `RNASEQ_PREPARE`, `FUNANNOTATE_RNASEQ.nf`,
     `FUNANNOTATE_PREDICTION.nf`, `utils.nf`, `pick_representative_strain.py` /
     `species_reuse_clusters.py`, `scripts/pick_rnaseq_representative_override.py`).
  4. Ticket list for implementation + validation (real end-to-end run on
     S. cerevisiae + a species with healthy Trinity to confirm no regression).
- Review outcome recorded here (Status line) after the session.

## Next steps

1. **(Done, 2026-08-26)** Alignment forensics (evidence item 4) and read
   attribution (item 5) — see "Forensics verdict" above.
2. **(Done, 2026-08-26)** Rep-pick audit (evidence item 6) — see "Rep-pick
   audit" above: KCTC won the intended BUSCO→N50 ranking (rank 1/1,308,
   replicated faithfully) on account of its mislabeled-C.-glabrata genome
   holding the group's max N50/tied-top BUSCO and pairing 100% to co-mislabeled
   `GCA_051107375.1`. **No code bug — the data (two mislabeled genomes inside
   the S. cerevisiae ANI group) is the failure.**
3. **(Immediate remediation)** Correct the taxonomy, then re-pick:
   a. Flag `GCA_026225675.1` (KCTC 13826BP) and `GCA_051107375.1` (MRD-KRBAY)
      in `samples.csv` as mislabeled *C. glabrata* / *Nakaseomyces glabratus* —
      remove them from the S. cerevisiae group (own species tag so they get
      annotated as their own unit, independently), NOT the S. cerevisiae
      representative.
   b. **No longer requires manual clearing** (superseded 2026-08-26 by the
      Option 2 extension fix): `MERGE_ANI`/`REPORT_ANI` now use `publishDir` +
      normal input-hash caching instead of `storeDir`, so once
      `samples.csv` changes KCTC/MRD-KRBAY's group membership, re-running
      `compare_ani` will automatically detect the changed genome list for
      `Saccharomyces_cerevisiae`/`Saccharomyces` (losing 2 members) and
      `Nakaseomyces_glabratus`/`Nakaseomyces` (gaining 2) and regenerate both
      — no manual `results/ANI/skani/{SPECIES,GENUS}/<group>/` directory
      move-aside step needed before rerunning. (Historical note: before this
      fix landed, the four directories —
      `results/ANI/skani/SPECIES/Saccharomyces_cerevisiae/`,
      `results/ANI/skani/GENUS/Saccharomyces/`,
      `results/ANI/skani/SPECIES/Nakaseomyces_glabratus/`,
      `results/ANI/skani/GENUS/Nakaseomyces/` — would have needed manual
      clearing first; this is now handled automatically.) Re-run `compare_ani`
      + `PICK_REPRESENTATIVE_STRAIN` for S. cerevisiae; the picker then
      selects a genuine S. cerevisiae (S288C-class, cluster-1; SX6 or
      NCYC478-class per the "Remediation representative candidate" analysis
      above).
   c. Re-run `RNASEQ_PREPARE` (rebuilds the shared Trinity from the existing
      BY4743 reads — already shown to align 99.8% to S288C) +
      `FUNANNOTATE_TRAIN` for the species.
4. Run the expert review session per the (now-narrowed) plan above; capture
   rulings on the workflow-hardening options and on whether a
   representative-species sanity check (new, see "Recommendation summary")
   should be added as Option 6.
5. Implement Options 1+2+3(+6) (+ 4 if gated in for future genuinely-divergent
   cases), add params/changesets per deliverables, validate end-to-end, update
   this doc's Status.

## Open questions (preliminary list)

- **Resolved by forensics (2026-08-26):** whose reads built the May 31
  Trinity and against which genome (BY4743 reads, correctly S288C-family;
  the *representative genome*, not the reads, was the mismatch); Q-BIO-1
  (KCTC's divergence — it's not intraspecies divergence, it's a different
  species).
- **Still open:** Q-BIO-2..5 (thresholds/general-case), Q-SW-1..6 (above).
- ~~Was KCTC picked by ANI-coverage-default rather than BUSCO rank?~~ —
  **Resolved (2026-08-26):** it was picked by BUSCO→N50 rank exactly as
  designed (rank 1/1,308; max N50, tied-top complete_pct in
  `busco_genome.parquet`); the absence from `tables/BUSCO.csv.gz` is
  irrelevant to the picker.
- Does the relaxed `pasa_shared_*` tier apply to the *representative* anywhere
  in the current code (believed no — reps always `stringent`)? Given the
  forensics verdict, this is now moot for KCTC specifically but still relevant
  as a general-case question.

## References (evidence / code locations)

- Failure: `work/funannotate/d0/eb840120d9c9d5f9109fb06220d2c0/`
- Assignments: `genome_annotation/_reuse_assignments/repr_assignments.tsv`,
  `abinitio_reuse_assignments.csv`
- Shared Trinity (quarantined): `rnaseq_data/Saccharomyces_cerevisiae.trinity-GG.fasta`
  (mtime 2026-05-31 20:59) — moved to
  `misc/quarantine/Saccharomyces_cerevisiae.trinity-GG.fasta_mislabeled_rep_2026-08-26`;
  missing log `rnaseq_data/Saccharomyces_cerevisiae.funannotate-trinity.log`
- Reads: `rnaseq_reads/Saccharomyces_cerevisiae_norm_R1/R2.fastq.gz` (May 26;
  SRR3491422 = BY4743 / S288C background)
- Forensics (evidence): `/scratch/jstajich/27800197/opencode/rnaseq_forensics/`
  (KCTC/S288C/Y10/kudr/eub references, PAFs); NCBI Datasets API report for
  `GCA_026225675.1` (atypical / ANI mismatch / best = Nakaseomyces glabratus)
- Empty query: `rnaseq_reads/sra_query/Saccharomyces_cerevisiae.sra_query.csv`
- Samples: `samples.csv` rows `GCA_026225675.1_ASM2622567v1` (KCTC) +
  `GCA_051107375.1_ASM5110737v1` (MRD-KRBAY), both labeled `Saccharomyces cerevisiae`
- ANI data (rep-pick audit): `results/ANI/skani/SPECIES/all_pairs_merged.tsv`
  (KCTC's only pair = 100.00 vs MRD-KRBAY);
  `results/ANI/skani/SPECIES/Saccharomyces_cerevisiae/{Saccharomyces_cerevisiae_ANI_report.txt,
  genomic_names.tsv, batches/Saccharomyces_cerevisiae.full.ani.tsv}` (KCTC+MRD-KRBAY =
  Cluster 2, N=2, <90% to the 1,305-strain majority cluster)
- BUSCO: picker sources `tables/busco_genome.parquet` + `tables/asm_stats.parquet`
  (KCTC complete_pct 75.80, N50 1,117,196); `tables/BUSCO.csv.gz` has no row (irrelevant)
- `nextflow/modules/funannotate/rnaseq/RNASEQ_PREPARE/main.nf:8` (`storeDir`)
- `nextflow/modules/ani/report/MERGE_ANI/main.nf`,
  `nextflow/modules/ani/report/REPORT_ANI/main.nf` — **fixed 2026-08-26**:
  group-level `storeDir` (same invalidation gap as `RNASEQ_PREPARE`, gate
  inventory item 6) replaced with `publishDir` + normal input-hash caching,
  matching `nextflow/modules/ani/compare/{SKANI,MASH,SOURMASH,FASTANI}_COMPARE/main.nf`,
  which already used this pattern for the identical reason (see their
  in-file comments)
- `nextflow/conf/profile_ANI.config` (`skani_prefilter`, new 2026-08-26,
  mirrors `fastani_prefilter`); `nextflow/modules/ani/compare/SKANI_COMPARE/main.nf`
  (added `batch_tag` to its signature so multiple per-component invocations
  don't collide on filename); `nextflow/subworkflows/local/ANI_COMPARE_METHOD.nf`
  (new `PREFILTER_COMPONENTS` subworkflow, factored out of the fastani
  branch's mash-prefilter cascade so `skani_prefilter` can reuse it) — Option 7's
  divide-and-conquer prefilter now available for skani (the production
  `ani_method`), not just fastani; opt-in (`skani_prefilter=false` default),
  useful at GENUS-scale groups where a single skani-triangle-over-everything
  job would be wasteful
- `scripts/find_ani_label_mismatches.py` (Option 7a, implemented 2026-08-26) —
  reads `results/ANI/skani/SPECIES/all_pairs_merged.tsv` + `samples.csv`,
  flags claimed-species groups that split into multiple ANI-connected
  components; first run output not yet triaged, see Option 7a discussion above
- `nextflow/modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf`
  (`train_min_trinity_transcripts`)
- `nextflow/modules/funannotate/utils.nf`
  (`loadAbinitioReuseMap` ~138-176; `loadRnaseqRepresentativeOverride` ~178-211,
  Ascochyta_rabiei precedent in its header comment)
- `nextflow/subworkflows/local/FUNANNOTATE_RNASEQ.nf` (header 1-9;
  `repr_ch` 217-233; `pasaTierFor` 272-279; `.combine` 284; train gate 295-299;
  `predict_no_rnaseq` 300)
- `nextflow/subworkflows/local/FUNANNOTATE_PREDICTION.nf` (`predict_scope` 45;
  `forceIndependentSet`; `predict_blocked_awaiting_representative.tsv`)
- `nextflow/conf/profile_funannotate.config`
  (`pasa_shared_*` ~161-163; `train_min_trinity_transcripts` ~191;
  `rnaseq_representative_override_csv` ~200; `ani_reuse_threshold` 99.0 ~215)
- `params_predict_representatives.yaml`
- `scripts/pick_rnaseq_representative_override.py`; `scripts/check_rnaseq_training.py`;
  `scripts/fix_low_trinity.py`
- Rescue precedent: relaxed `pasa_shared_min_avg_per_id 85 / min_pct_aligned 70 /
  num_bp_splice 1` saving Hansenula (19/14,319 transcripts)
