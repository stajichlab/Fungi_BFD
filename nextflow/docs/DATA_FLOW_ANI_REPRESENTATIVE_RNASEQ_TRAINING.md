# Data flow: ANI clustering → representative pick → RNA-seq training → predict

Reference doc, written 2026-08-28. Code-grounded (file:line references throughout)
so it can be re-verified against the repo, not treated as a standing description.
Companion docs, more narrowly scoped:

- `HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md` — why `FUNANNOTATE_TRAIN` runs on every
  strain even under `--predict_scope representative_only`.
- `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` — the KCTC misidentified-representative
  incident; proposes an alignment-floor gate and representative-species validation
  (Options 1-3), not yet implemented. Read that doc's "gate inventory" section
  before touching any of the logic documented here.
- `HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md` — proposal building on this doc, for
  `Genus_x_Genus` hybrid-cross species specifically.

## Why this exists

Two independent problems get solved by mostly-shared machinery, and the code
comments across `FUNANNOTATE_RNASEQ.nf`/`FUNANNOTATE_PREDICTION.nf`/`utils.nf`
assume the reader already holds the whole picture in their head. This doc lays
out the full pipe, stage by stage, with the actual code path for each decision
point.

**Problem 1 — ab-initio parameter reuse** ("does gene-model training need to
run independently for every strain, or can closely related strains share
AUGUSTUS/SNAP/GeneMark parameters?"). Driven by ANI clustering +
`PICK_REPRESENTATIVE_STRAIN`.

**Problem 2 — RNA-seq/PASA evidence sharing** ("does every strain need its own
RNA-seq fetched and assembled, or can strains of the same species share one
Trinity assembly for PASA training?"). Driven by the *same* representative
pick, but is a functionally separate reuse axis with its own tiering and its
own override.

## Stage 0 — samples.csv → per-genome metadata

`subworkflows/local/ANI_SAMPLES.nf` (used by both `compare_ani` and
`query_ani` entrypoints) reads `samples.csv`, filters by taxon rank/suppress
list, resolves each row's `ASMID` to a genome file in `params.genome_dir`, and
emits `tuple(group_key, meta)` where `group_key` is `sanitizeTag(row[compare_rank])`
(usually `SPECIES`, sometimes `GENUS`/`FAMILY` for higher-rank comparisons) and
`meta` carries `species`, `strain`, `asmid`, `genome`, taxonomy columns, and
`is_query`.

`funannotate.nf` separately builds `predict_genome_ch` (the 10-tuple
`out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid`)
that everything downstream in this doc actually consumes. `species_tag` (used
throughout `FUNANNOTATE_RNASEQ.nf`) is always `species.replaceAll(/\s+/, '_')` —
computed inline wherever needed, not carried as a channel field.

## Stage 1 — ANI comparison

`ANI_COMPARE_METHOD.nf` runs one of `skani` (production default) / `mash` /
`sourmash` / `fastani` over each `compare_rank` group (all-vs-all within
`SPECIES`, or within `GENUS` for cross-species placement). Production sizing:
`SPECIES`-scope groups top out around 1,300 genomes (S. cerevisiae); a single
`skani triangle` call handles that directly. `skani_prefilter`/`fastani_prefilter`
(mash-sketch → connected-components → confirm-within-component) exist for
`GENUS`-scope groups where a naive all-vs-all would be wasteful — off by
default, opt-in per `ani_method`.

Output: per-group pairwise ANI tables, merged by `MERGE_ANI`/`REPORT_ANI` into
`results/ANI/skani/{SPECIES,GENUS}/**` and ultimately
`results/ANI/skani/SPECIES/all_pairs_merged.tsv` (whole-dataset). As of
2026-08-26 `MERGE_ANI`/`REPORT_ANI` use `publishDir` + normal Nextflow
input-hash caching (not `storeDir`) specifically so a group's report
regenerates when its membership changes — see
`DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` gate 6 for why that mattered.

Relevant thresholds (`conf/profile_ANI.config`):
- `ani_cluster_threshold = 95.0` — used by `scripts/find_ani_label_mismatches.py`
  for post-hoc cluster-vs-label auditing (not in the live Nextflow path).
- `ani_outlier_threshold = 90.0` — same script, "HIGH-confidence misidentification"
  cutoff.

## Stage 2 — representative pick (`PICK_REPRESENTATIVE_STRAIN`)

`modules/ani/report/PICK_REPRESENTATIVE_STRAIN/main.nf` runs
`bin/pick_representative_strain.py` over the merged ANI TSV + `tables/busco_genome.parquet`
+ `tables/asm_stats.parquet`. Per species (only species with ≥2 strains are
candidates):

1. Restrict to strains that have **at least one ANI pairing** in the merged
   TSV (a strain with zero ANI coverage cannot be picked — see the KCTC
   postmortem for why "has ANI coverage" is not itself sufficient).
2. Rank by BUSCO `complete_pct` (descending), tiebreak by assembly N50
   (descending), final tiebreak alphabetical by `out`.
3. The top-ranked strain becomes `is_representative=True` for the species.
4. Every other strain's `ani_to_representative` is looked up from the ANI TSV;
   `reuse_eligible = ani_to_representative >= ani_reuse_threshold` (default
   `99.0`, `conf/profile_funannotate.config:225`) — note this is a **different**
   threshold from the `stringent`/`relaxed`/`skip` PASA tiers in Stage 4.
   **Fail-closed**: no ANI row for a strain (sparse ANI methods omit
   low-identity pairs) → `ani_to_representative` blank → not eligible. A blank
   value is a divergence signal, not "not yet computed" (see
   `utils.nf:158-163`).

Output: `abinitio_reuse_assignments.csv` (`out, species, reuse_eligible,
is_representative, ani_to_representative, ...`) and `repr_assignments.tsv`
(human-readable audit copy). Loaded once per pipeline run via
`loadAbinitioReuseMap()` (`modules/funannotate/utils.nf:138`) into
`out -> [species, reuse_eligible, is_representative, ani_to_representative]`
(representative's own `ani_to_representative` hardcoded to `100.0`), passed as
a plain Groovy map (not a channel) into both `FUNANNOTATE_RNASEQ` and
`FUNANNOTATE_PREDICTION`.

Whole feature gated by `params.share_abinitio_params` (off → empty map → every
strain trains independently, unchanged from pre-reuse behavior).

`rnaseq_representative_override.csv` (`species_tag -> out`, loaded by
`loadRnaseqRepresentativeOverride()`, `utils.nf:190`) can force a **different**
strain to be the RNA-seq anchor specifically, without touching the ab-initio
representative pick — see Stage 4. Populated by
`scripts/pick_rnaseq_representative_override.py`, presently used for one
Ascochyta_rabiei-class case (ANI+BUSCO picks an assembly whose real RNA-seq
barely aligns).

## Stage 3 — ab-initio parameter reuse at predict time (independent of RNA-seq)

`subworkflows/local/FUNANNOTATE_PREDICTION.nf` classifies every incoming
strain via `abinitioReuseMap[out]` into `representative` / `eligible_sibling`
/ `independent` (lines 70-82). This branching is **orthogonal to RNA-seq** —
it governs whether `FUNANNOTATE_PREDICT` reuses the representative's trained
AUGUSTUS/SNAP/GeneMark `parameters.json` (`BACKFILL_ABINITIO_PARAMS`) instead
of training its own ab-initio models from scratch. `--predict_scope`
(`all` default, or `representative_only`) further restricts *which* strains
actually run `FUNANNOTATE_PREDICT` in a given invocation — see
`HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md` for the important gotcha that
`predict_scope` does **not** gate `FUNANNOTATE_TRAIN` (Stage 4 always runs for
every strain with RNA-seq, regardless of predict_scope).

`forceIndependentSet`/`forceIndependentGenemarkSet` are maintainer-supplied
species/strain lists that opt specific groups out of ab-initio reuse even when
`reuse_eligible=True` — an escape hatch parallel to the RNA-seq override.

## Stage 4 — RNA-seq acquisition and shared-Trinity training (`FUNANNOTATE_RNASEQ.nf`)

Gated by `params.run_sra_fetch`; off → every strain passes through untrained
(9-tuple, no RNA-seq, straight to Stage 3/predict as ab-initio-only).

When on:

1. **Per-species SRA query** (`SRA_QUERY_BATCH`): one query per `species_tag`
   (all strains of a species share the same taxon-level SRA search — reads are
   never queried per-strain). Results cached to
   `rnaseq_reads/sra_query/<species_tag>.sra_query.csv`;
   `params.skip_sra_query` reuses cached CSVs without submitting new jobs.
2. **Classification and fetch**: each species' CSV is routed to `SRA_FETCH`
   (has ≥1 usable PAIRED accession), `SRA_FETCH_SE` (SE-only, or PE forced to
   SE via `rnaseq_blacklist.csv`'s `SE_trinity` action, or genuinely
   SINGLE-layout when `enable_single_end=true`), or `WRITE_EMPTY_READS` (no
   usable accessions at all) — `FUNANNOTATE_RNASEQ.nf:152-166`.
3. **Representative pick for RNA-seq purposes**
   (`FUNANNOTATE_RNASEQ.nf:217-233`): all strains of a species are
   `groupTuple`'d into one row; the RNA-seq anchor is
   `rnaseqRepOverride[species_tag]` if present, else the same
   `abinitioReuseMap[out].is_representative` pick as Stage 2, else index 0 of
   the group (only when neither exists). **Exactly one** strain per species is
   chosen — this is the strain whose genome+reads get assembled.
4. **`RNASEQ_PREPARE`** (`modules/funannotate/rnaseq/RNASEQ_PREPARE/main.nf`)
   runs `funannotate train --stop_after_trinity --no_trimmomatic` on that one
   representative, producing `rnaseq_data/<species_tag>.trinity-GG.fasta` —
   `storeDir`-cached keyed on `species_tag` alone (no representative-identity
   invalidation yet; see `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` Option 2).
   A species whose representative has zero reads gets an empty FASTA written
   directly by the Nextflow driver (no SLURM job).
5. **Fan-out to every strain** (`FUNANNOTATE_RNASEQ.nf:283-320`): the single
   shared Trinity-GG is joined back to every strain of the species by
   `species_tag`. **Every** strain — not just the representative, and
   regardless of `--predict_scope` — runs `FUNANNOTATE_TRAIN --trinity
   <shared assembly>` to PASA-align that shared transcriptome against its own
   genome. This is deliberate pre-staging (see
   `HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md`).
6. **PASA tiering** (`pasaTierFor`, `FUNANNOTATE_RNASEQ.nf:272-279`) sets how
   strict PASA's alignment thresholds are, based on
   `abinitioReuseMap[out].ani_to_representative` (representative itself always
   `100.0`):
   - `>= 97.0` → `stringent` (PASA defaults).
   - `90.0–97.0` → `relaxed` (`pasa_shared_min_avg_per_id=85`,
     `pasa_shared_min_pct_aligned=70`, `pasa_shared_num_bp_splice=1` —
     `conf/profile_funannotate.config:171-172`).
   - `< 90.0` or blank/no-ANI-signal → `skip`: the strain bypasses
     `FUNANNOTATE_TRAIN` entirely and flows straight to `predict_no_rnaseq`
     (ab-initio-only), same path as a genuinely RNA-seq-less strain.
   - When `abinitioReuseMap` is empty (ANI-reuse feature off), every
     shared-Trinity strain falls back to `relaxed` unconditionally — there's
     no per-strain ANI signal to tier on.
   - **Known gap**: the representative's own row always tiers `stringent`
     (its `ani_to_representative` is hardcoded `100.0`) — there is no
     "this representative's own evidence is unusable" tier. See
     `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` for the incident this produced
     and the proposed alignment-floor fix (not yet implemented).
7. Additional gates inside `FUNANNOTATE_TRAIN`
   (`modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf`):
   `train_min_trinity_transcripts` (default 2000, set 0 to disable) skips
   training when the shared Trinity has too few transcripts — a signal of
   "assembled against the wrong reference strain" per that module's header
   comment. Per-strain staleness re-triggers training when the genome or
   RNA-seq reads are newer than the existing `funannotate_train.pasa.gff3`.

## Stage 5 — predict

`FUNANNOTATE_PREDICT` runs per strain, reusing shared ab-initio params when
Stage 3 marked it eligible, and reading whatever PASA training Stage 4
produced (real evidence, or nothing if the strain routed to
`predict_no_rnaseq`/`skip`). `BACKFILL_ABINITIO_PARAMS` writes the shared
`parameters.json` store from a representative's own completed prediction so
later runs (or later strains in the same run) can reuse it.

## Two reuse axes, side by side

| | Ab-initio params (Stage 2/3) | RNA-seq/PASA (Stage 4) |
|---|---|---|
| Unit shared | AUGUSTUS/SNAP/GeneMark `parameters.json` | Trinity-GG transcriptome FASTA |
| Picked by | ANI+BUSCO representative (`PICK_REPRESENTATIVE_STRAIN`) | Same pick, unless `rnaseq_representative_override.csv` overrides it |
| Eligibility threshold | `ani_reuse_threshold` (99.0) → binary `reuse_eligible` | `pasaTierFor` (97.0/90.0) → 3-way `stringent`/`relaxed`/`skip` |
| Escape hatch | `forceIndependentSet` / `forceIndependentGenemarkSet` | `rnaseq_representative_override.csv` |
| Gated by | `params.share_abinitio_params` | `params.run_sra_fetch` |
| Failure mode when off | every strain trains its own ab-initio models | every strain trains without shared Trinity (still runs its own `FUNANNOTATE_TRAIN` if it has its own reads, otherwise ab-initio-only) |

## Known open gaps (as of 2026-08-28)

1. Representative-path alignment floor — `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md`
   Options 1-3, not implemented.
2. `RNASEQ_PREPARE`'s `storeDir` has no representative-identity invalidation
   sentinel (same doc, Option 2, ANI-side already fixed, RNA-seq side not).
3. Hybrid-cross species (`Genus_a species_a x Genus_b species_b[...]`) are
   treated as an ordinary species for both reuse axes — no special handling.
   See `HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md`.
4. `predict_scope`-aware filtering of the Stage 4 fan-out (train siblings
   lazily instead of pre-staging all of them) — noted as a possible future
   option in `HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md`, not requested/implemented.
