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

## Fixes applied while preparing this PR

- `nextflow/subworkflows/local/FUNANNOTATE_GENOME_PREP.nf` includes three
  modules — `GENOME_CLEAN`, `GENOME_CLEAN_BATCH`, `MASKREPEAT_TANTAN_RUN` —
  whose `main.nf` files existed on disk but had never been `git add`ed, so
  `main.nf` failed to parse from a fresh checkout of this branch. Added them
  to this PR; `nextflow main.nf --help` now parses and runs cleanly.
- `BFD_MERGE.nf`'s run-mode `MERGE_INTERGENIC`/`MERGE_GENE_STATS`/
  `MERGE_ASM_STATS` calls passed `.collect()`-ed channels into `toManifest()`,
  which expects one file per emission; `.collect()` emits a single list
  instead, so `toManifest()` crashed trying to `file()` the whole list
  (`No signature of method: nextflow.util.ArrayBag.getFileSystem()`). Fixed
  by dropping `.collect()`, matching the working `MERGE_AA_FREQ`/
  `MERGE_CODON_FREQ` calls beside them.

## Test status

- `nextflow/run_test.sh` stub-run (BFD test profile) — **passing** as of this
  PR, including the `toManifest()` fix above.
- Real HPCC smoke test of a full `BFD.nf`/`funannotate.nf` run — still
  pending before merge.

## Full commit list

See `git log main..refactor_modules` for the complete 70-commit history.

---

# `main` — 2026-08-26: taxonomy misidentification incident + ANI/representative-picking hardening

Separate, later changeset — not part of the `refactor_modules` PR above.
Full narrative, evidence, and open items: `nextflow/docs/DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md`.

## What happened

`FUNANNOTATE_TRAIN` crashed for `Saccharomyces_cerevisiae_KCTC_13826BP`
(PASA assigned 36 of 5,186 transcripts to loci) because that genome had been
picked as the species' shared RNA-seq/ab-initio representative. Root-caused
via NCBI's own ANI-taxonomy check + independent minimap2 alignment: the
genome (`GCA_026225675.1`, "KCTC 13826BP") and a second genome
(`GCA_051107375.1`, "MRD-KRBAY") are **not** *Saccharomyces cerevisiae* —
they are misidentified *Nakaseomyces glabratus* (*Candida glabrata*, NCBI
taxid 5478; 99.93% ANI to *N. glabratus* vs. 87.2%/failing to *S. cerevisiae*).
The representative-picker itself worked exactly as designed (ranked KCTC
1st of 1,308 by BUSCO completeness then N50); the two mislabeled genomes
were the actual defect, not the pick logic. No downstream ab-initio
parameters were corrupted — the pipeline's fail-closed blank-ANI semantics
already prevented that — but the shared RNA-seq training resource was wrong
for the whole species group, and the underlying gap (nothing validates a
representative's *species*, and ANI-group storeDir caches don't invalidate
on membership change) is generalizable, not KCTC-specific.

## Data/curation fix

- `data/curation/overrides.csv`: added `SPECIES_IN`/`GENUS`/`SPECIES`/
  `NCBI_TAXONID` overrides for both `GCA_026225675.1_ASM2622567v1` and
  `GCA_051107375.1_ASM5110737v1`, reclassifying them from *Saccharomyces
  cerevisiae* to *Nakaseomyces glabratus* (taxid 5478, confirmed against the
  local NCBI taxdump via `taxonkit`). `samples.csv` regeneration
  (`scripts/create_samples_file.py --overrides data/curation/overrides.csv`)
  and copying the result into `Fungi_BFD_runs/samples.csv` (a separate,
  non-symlinked copy that Nextflow actually reads at runtime) are still
  **pending** — not yet executed as of this writing.

## Code changes

1. **ANI aggregate `storeDir` → `publishDir`** (the actual bug preventing
   automatic recovery from the above): `nextflow/modules/ani/report/MERGE_ANI/main.nf`
   and `nextflow/modules/ani/report/REPORT_ANI/main.nf` used `storeDir`,
   which only checks output-path *existence* — a species/genus group whose
   membership changed (this incident's reclassification, or any future one)
   would have kept silently serving its pre-change merged ANI table/report
   forever, requiring a manual directory move-aside before every affected
   rerun. Switched both to `publishDir` + Nextflow's normal input-hash
   caching, matching the pattern the pairwise `*_COMPARE` modules
   (`SKANI_COMPARE`/`MASH_COMPARE`/`SOURMASH_COMPARE`/`FASTANI_COMPARE`)
   already used for this exact reason. `REPORT_ANI` previously had **both**
   `storeDir` and `publishDir` on the same path — `storeDir`'s existence
   check silently won regardless, making the `publishDir` vestigial; now
   `storeDir` is removed and `publishDir` is the only cache mechanism.
   Verified via stub-run: pre-fix a second invocation logged
   `[skipping] Stored process > COMPARE_ANI:REPORT_ANI`; post-fix the
   identical invocation shows the process actually executing.
2. **skani divide-and-conquer prefilter** (`nextflow/conf/profile_ANI.config`,
   `nextflow/modules/ani/compare/SKANI_COMPARE/main.nf`,
   `nextflow/subworkflows/local/ANI_COMPARE_METHOD.nf`): a
   mash-prefilter → connected-components → confirmatory-comparison cascade
   already existed for `fastani` (gated behind `fastani_prefilter`) but was
   unreachable for `skani`, the actual production `ani_method`. Added a new
   `skani_prefilter` param (default `false`, independent toggle) and
   extracted the shared prefilter/component-explosion logic into a new
   `PREFILTER_COMPONENTS` subworkflow (not a Groovy closure — Nextflow
   disallows process calls inside closures) that both the `fastani` and
   `skani` branches now call, removing the previous duplication.
   `SKANI_COMPARE`'s signature gained a `batch_tag` (mirroring
   `FASTANI_COMPARE`) so per-component invocations don't collide on
   filename; whole-group calls pass `"full"`, preserving today's behavior
   byte-for-byte when the new flag is off (the default). Validated via
   `-stub-run`: default path unchanged, new path fires the full
   sketch→prefilter→components cascade with no channel/type errors. Opt-in;
   intended for `GENUS`-scale comparisons where a single skani-triangle
   job over an entire (large) genus would be wasted compute, and as
   infrastructure a future genus-scale version of the proactive
   misidentification sweep (below) can reuse.
3. **Proactive dataset-wide misidentification sweep** — new script
   `scripts/find_ani_label_mismatches.py`. Reuses the already-computed
   whole-dataset merged ANI table (`results/ANI/skani/SPECIES/all_pairs_merged.tsv`,
   2,195,284 pairs / 1,770 species groups) to flag any claimed-species group
   whose ANI graph splits into more than one connected component at the
   95% cluster threshold — exactly the KCTC/MRD-KRBAY signature (a 2-genome
   cluster at 0.00% ANI to the 1,305-genome majority), reproduced
   automatically with no training crash needed as a tripwire. Tiers findings
   `HIGH`/`REVIEW` confidence using the existing `ani_outlier_threshold`
   (90%) convention. Zero new Nextflow compute (~10s over the full 23k-genome
   dataset). First run: 170 HIGH-confidence + 240 REVIEW-tier candidates
   across 219 species, including 2 dataset-wide candidates beyond
   KCTC/MRD-KRBAY not yet forensically confirmed (`Candidozyma auris`,
   `Macrophomina phaseolina`) — triage tracked as `todo/TODO_REGISTRY.md` T-030.
4. **Documentation fixes**: corrected two factual errors in
   `nextflow/docs/DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` inherited from an
   earlier draft — `pasaTierFor()` actually lives in
   `FUNANNOTATE_RNASEQ.nf` (not `utils.nf`), with real cutoffs `ANI>=97 ->
   stringent, 90<=ANI<97 -> relaxed, <90/null -> skip` (not the
   `ani_reuse_threshold=99.0` value, which is an unrelated parameter feeding
   `reuse_eligible` in `PICK_REPRESENTATIVE_STRAIN`/`BACKFILL_ABINITIO_PARAMS`).

## Not yet done

- `samples.csv` regeneration + sync to `Fungi_BFD_runs/samples.csv` (data fix
  above).
- Re-running `compare_ani` → `PICK_REPRESENTATIVE_STRAIN` →
  `RNASEQ_PREPARE` → `FUNANNOTATE_TRAIN` for *S. cerevisiae* against a
  correctly-picked representative (candidates identified: `GCA_003277715.1`
  "SX6", the picker-native winner, vs. `GCA_947344615.1`/`GCA_003273825.1`
  for far better `reuse_eligible` coverage — see the plan doc's
  "Remediation representative candidate" analysis for the trade-off).
- `RNASEQ_PREPARE`'s own representative-identity `storeDir` sentinel (Option
  2's original scope) — still design-only.
- Option 1 (representative-path alignment floor) and Option 3 (data-driven
  representative exemption) from the plan doc — still design-only, pending
  expert review per the doc's (now-narrowed) evaluation plan.
- Triage of the 410 candidates from `find_ani_label_mismatches.py`'s first
  run (T-030).

---

# `main` — 2026-08-28: paralogoscope pipeline (wgd whole-genome duplication dating)

New standalone Nextflow pipeline that dates whole-genome duplication (WGD)
events per species using [wgd](https://github.com/arzwa/wgd) (Tiley et al.
2021) inside the `container_wgd2_complete` SIF.

## What was added

- **Pipeline**: `nextflow/workflows/paralogoscope.nf` +
  `nextflow/subworkflows/local/PARALOGOSCOPE.nf` — per-species
  `wgd dmd` (DIAMOND all-vs-all + MCL paralog families) → `wgd ksd`
  (per-family Ka/Ks via mafft + codeml, Ks distribution plots) → optional
  `wgd syn` (i-ADHoRe synteny anchors, gated behind `--run_wgd_syn true`).
- **Modules**: `nextflow/modules/comparative/wgd/{WGD_DMD,WGD_KSD,WGD_SYN}/main.nf`.
- **Final merge** (`MERGE_WGD_KSD`, label 'merge'): after the last `wgd ksd`
  task, globs *all* ks.tsv published so far off disk and streams them into
  `tables/wgd.ks.parquet` (`bin/merge_wgd_ks.py`; 16-column schema — genome +
  species_prefix + pair/family/g1/g2/gene1/gene2/node + N/S/dN/dN/dS/dS/
  alignmentlength/t — zstd; ~0.8–1 GB for the full dataset vs ~45 GB raw) +
  per-genome `tables/wgd.ks.summary.parquet`, following the BFD
  MERGE_*/tablesDir() convention. Disk-globbing (and a completion-only gate)
  makes it wave- and `-resume`-safe: a fully-cached rerun still rebuilds the
  merged table. (`node` = wgd tree-node id the pair was dated at; added so
  the Ks-peak summariser can use wgd's own node-averaged mixture models.)
- **BFD duckdb exposure**: `nextflow/bin/build_BFD_duckDB.sh` gained guarded
  optional tables — `wgd_ks` (indexed on species_prefix/genome) and
  `wgd_ks_summary` (unique on genome) are loaded into `db/BFD.duckdb` only
  when `tables/wgd.ks*.parquet` exist (paralogoscope-run pipelines only);
  plain BFD rebuilds print a SKIP line instead of erroring.
- **Per-genome Ks-peak summary framework** (the "number of duplicates + mean
  Ks of the peaks" ask): `analysis/WGD_PERFORMANCE_ANALYSIS/scripts/{
  wgd_ks_mix_fit.py,build_wgd_ksd_summary.py}` derive
  `tables/wgd_ksd_summary.parquet` + `tables/wgd_ksd_density.parquet` from
  the merged `wgd.ks*.parquet` using **wgd's own mixture models** (the CLI's
  `mix`/`peak` never export component means, so the worker calls
  `wgd.mix.fit_gmm` inside the SIF and captures them). wgd.ks schema gained
  `node` (16 cols) so the node-averaging filter reproduces exact `wgd mix`
  semantics. Validated on the 2-genome synb pair (fits match a direct
  in-container fit; arxii 3 peaks 0.080/1.559/3.500, pasadenensis 4 peaks
  0.053/0.853/2.449/3.766). Full-4,365 execution tracked under T-032.
- **Profile/launcher**: `nextflow/conf/profile_paralogoscope.config` +
  `nextflow/run_paralogoscope.sh` (SLURM submit, **preempt** partition —
  `-A preempt -p preempt`, 8 cpus / 32 GB / 24 h per wgd task).
- **Inputs** follow the BFD structure (`input/cds/{Species_Strain}.cds-transcripts.fa`,
  `input/gff3/` for syn); selection mirrors the other pipelines (`--taxon`,
  `--group`, `--ignore`, `--n_test`). Synthetic test fixtures in
  `nextflow/tests/data/input_wgd/{cds,gff3}`.

## Container/runtime fixes (real-data validation)

- **Workdir symlink staging**: Nextflow stages inputs as symlinks into the
  task workdir pointing at `/bigdata/...`; the SIF originally only bound
  `${PWD},${projectDir},$TMPDIR`, so real (non-`tests/`) inputs failed with
  `Path ... does not exist`. All three modules now bind
  `/bigdata:/bigdata,/rhome:/rhome`.
- **Host user-site shadowing**: the container's Python resolved
  `/rhome/<user>/.local/lib/python3.9/site-packages` before `/opt/conda`, so
  a newer host Biopython shadowed the container's Bio 1.76 (`Bio.Alphabet`
  import error). Fixed with `--env PYTHONNOUSERSITE=1`.
- **OpenMPI env hygiene** (syn): `wgd syn` needs
  `LD_LIBRARY_PATH=/opt/conda/lib` and `SLURM_*/PMI_*` env vars stripped
  inside the container, plus `-f mRNA -a ID` for funannotate-style GFF3s.
  "No anchors found" is a legitimate rc=0 (no `anchors.csv` written); the
  module guard tolerates a non-empty `iadhore-out/anchorpoints.txt` instead.

## Output layout (BFD hash sub-folder convention)

`wgd dmd`/`ksd`/`syn` outputs publish via `storeDir` (BFD-style SHA-1 hash
bucket keyed on the stable LOCUSTAG, 256 buckets), so no per-genome dir
accumulates too many files and reruns reuse stored results:

```
paralogoscope/wgd_dmd/{bucket}/{id}.cds-transcripts.fa.tsv
paralogoscope/wgd_ksd/{bucket}/{id}.cds-transcripts.fa.tsv.ks.tsv
paralogoscope/wgd_ksd/{bucket}/{id}.cds-transcripts.fa.tsv.ksd.pdf
paralogoscope/wgd_syn/{bucket}/anchors.csv
```

Natural wgd filenames (which embed the source cds id) are preserved so a
re-annotated sample can never be served a stale `storeDir` hit — the
staleness hazard BFD's `clearIfStale` machinery exists for is avoided by
construction.

## Test status

- Stub-run 6/6 (dmd×2, ksd×2, syn×2) clean.
- Real-data smoke test (2 Aaosphaeria genomes): `wgd dmd` 2/2 exit 0,
  published TSVs confirmed; `wgd ksd` validated in-container (slow on real
  genomes, still running) — full run proceeds via SLURM.

