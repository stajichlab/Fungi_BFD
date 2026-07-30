# CHANGES — `refactor_modules` → `main`

Summary of the `refactor_modules` branch (70 commits ahead of `main`, ~213 files
changed) being merged into `main`.

## Highlights

### 1. Monolithic Nextflow workflows decomposed into `modules/` + `subworkflows/`

The previous single-file pipelines (`nextflow/BFD.nf`, `nextflow/funannotate.nf`,
`nextflow/compare_ANI.nf`, `nextflow/query_ANI.nf`, `nextflow/comparative_genomics.nf`,
`nextflow/earlgrey_mask.nf`, `nextflow/phyling.nf` — thousands of lines each) were
retired and replaced with:

- `nextflow/main.nf` — single entry point
- `nextflow/workflows/*.nf` — thin top-level workflow definitions (`BFD.nf`,
  `funannotate.nf`, `compare_ANI.nf`, `query_ANI.nf`, `comparative_genomics.nf`,
  `earlgrey_mask.nf`, `backfill_abinitio.nf`)
- `nextflow/subworkflows/local/*.nf` — reusable stages (`FUNANNOTATE_GENOME_PREP`,
  `FUNANNOTATE_PREDICTION`, `FUNANNOTATE_RNASEQ`, `FUNANNOTATE_ANNOTATION`,
  `ANI_SAMPLES`, `ANI_COMPARE_METHOD`, `ANI_REPRESENTATIVE_SELECT`,
  `BFD_GENOME_STATS`, `BFD_FUNCTIONAL`, `BFD_MERGE`, `PHYLING_ALIGN`,
  `PREPARE_COMPARATIVE`, `CLUSTER_MCL`/`CLUSTER_MMSEQS2`/`CLUSTER_ORTHOFINDER`,
  `INPUT_SETUP`)
- `nextflow/modules/{BFD,ani,common,comparative,earlgrey,funannotate,phyling}/**/main.nf`
  — one process per module (~90 new module files)

`nextflow/lib/SampleUtils.groovy` was dropped in favor of `nextflow/modules/common/utils.nf`
and `nextflow/modules/funannotate/utils.nf`.

Uses nf-schema for parameter and `samples.csv` validation (`nextflow/assets/schema_input.json`).

### 2. ANI (Average Nucleotide Identity) pipeline rework

- `compare_ANI` / `query_ANI` migrated to modules, sharing common comparison
  backends: MASH, skani (updated for the v0.3.x API), sourmash, FastANI.
- Representative-strain selection (`PICK_REPRESENTATIVE_STRAIN`,
  `ANI_REPRESENTATIVE_SELECT`) added to pick one strain per species/cluster for
  downstream annotation-parameter reuse, backed by DuckDB for fast load/query
  instead of SQLite.
- `asmid_filter` plumbed through `ANI_SAMPLES`/`compare_ANI`/`query_ANI` so
  `--asmid` filtering is now respected consistently everywhere.
- **cf93c83** — `ANI_SAMPLES` replaced up to two sequential S3
  `file().exists()` checks per genome with a single directory listing, resolved
  in-memory. Verified live on a phylum-scale run (Basidiomycota, 3173 genomes):
  eliminates 20-30+ minutes of serial S3 round-trip latency before the first
  compare task can start.
- **dc086bc** — Same fix applied to `BFD.nf`: its 8 per-genome input channels
  (pep/cds/gff/genome, including `resolveGenomeFile()`'s asm-stats/BUSCO-genome
  callers) now resolve against a cached `dirIndex(genomeDir)` map instead of
  doing up to three `.exists()` calls per genome. Not S3-latency-bound today
  (`profile_BFD.config` points these dirs at local paths), but `genome_dir`
  already shares `input_clean_genomes`'s S3 convention with the ANI k8s
  profile, so the redundant-stat cost is latent, not hypothetical.

### 3. Ab-initio parameter reuse across strains

- `species_reuse_clusters.py` / `backfill_abinitio_params.py` (new
  `nextflow/bin/` scripts) identify ANI-qualified representative strains per
  species and let related strains reuse trained AUGUSTUS/SNAP/GeneMark-ES
  parameters instead of retraining from scratch (`--share_abinitio_params`,
  `--abinitio_reuse_csv`, `--gene_prediction_shared_abinitio`).
- `pick_representative_strain.py` writes per-species
  `abinitio_reuse_assignments.{species}.csv` with an atomic-write + merge
  strategy to avoid clobbering concurrent runs.
- `FUNANNOTATE_RNASEQ` picks its RNA-seq representative from the same
  `abinitioReuseMap`.

### 4. Nautilus / Kubernetes (k8s) offload pilot for `compare_ANI`

New `nextflow/k8/` directory: RBAC/PVC manifests, head-pod config, staging
scripts (`stage_repo.sh`, `stage_data_ani.sh`), and a split
compute/gather-phase design (`run_ani_compute.nf`, `run_ani_gather.nf`) for
distributing ANI comparison work across the Nautilus (NRP) cluster, plus a
`QUICKSTART.md` cheat sheet and `README_compare_ani.md` documenting
operational gotchas (e.g. the namespace's hard 6h `activeDeadlineSeconds`
pod deadline, S3-vs-PVC-vs-pod-local path handling).

### 5. Test/data cleanup

- Removed ~111k lines of vendored AUGUSTUS species-model binaries
  (`lib/augustus/3.5/config/species/**`) that don't belong in version control.
- Reworked `nextflow/tests/data/**` fixtures (funannotate + ANI stub-run test
  data) and extended `run_test.sh` stub-run coverage to `compare_ani` and
  `funannotate`, including `--asmid` filtering assertions.

### 6. Misc fixes folded in

- Deadlocks and shell-escaping bugs in representative-strain picking fixed by
  moving string construction into Groovy instead of interpolating inside
  `script:`/`shell:` blocks.
- `BFD.nf` genome-stats channels now resolve by ASMID when the
  post-annotation `meta.id` naming convention isn't present (raw/clean
  assemblies), instead of silently dropping genomes.
- `CONCAT_ANI_TSVS` no longer loses the ANI manifest to `collectFile()`
  swallowing it.
- `MERGE_AA_FREQ` / `MERGE_CODON_FREQ` nondeterminism fixed.
- Sample identity threaded through `BFD.nf` as a `meta` map rather than ad hoc
  strings.

## Fix applied while preparing this PR

`nextflow/subworkflows/local/FUNANNOTATE_GENOME_PREP.nf` includes three
modules — `GENOME_CLEAN`, `GENOME_CLEAN_BATCH`, `MASKREPEAT_TANTAN_RUN` — whose
`main.nf` files existed on disk but had never been `git add`ed, so `main.nf`
failed to parse from a fresh checkout of this branch. Added them to this PR;
`nextflow main.nf --help` now parses and runs cleanly.

## Full commit list

See `git log main..refactor_modules` for the complete 70-commit history.
