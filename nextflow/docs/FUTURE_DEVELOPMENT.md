# Future development

A running list of things to consider later. Add new entries at the top of the list. Keep each entry short. State what was observed, and say what was not tested.

Entry format: date, area, what to consider, why, and status.

## 2026-10-21 — Container for gene-stats and tRNA steps (bedtools, pybedtools, pyfaidx, biopython)

- **Area:** `CALC_GENE_STATS`, `TRNASCAN` (modules/BFD).
- **Now:** both load the `biopython` and `bedtools/2.30.0` HPC modules. The decision is to keep the modules for now.
- **Consider:** one container with Python, biopython, pybedtools, pyfaidx and bedtools together. This removes the dependence on the HPC module setup.
- **Why:** the old `biopython/1.87` modulefile ran `export PATH=...` with a fixed list. It replaced `PATH` and removed `bedtools` from it. The module is fixed now and uses `prepend-path`. A container would avoid this class of problem.
- **Observed:** the cached image `depot.galaxyproject.org-singularity-bedtools-2.31.1--hf5e1c6e_0.img` has `bedtools` v2.31.1 only. It has no `python3`, so it has no `pybedtools`, `pyfaidx` or `Bio`.
- **Not done:** the cache has several `mulled-v2-*` images. I did not check whether one has all of the packages. Building an image from a small spec is the alternative.
- **Status:** open.

## 2026-10-21 — Bundle the short functional-annotation tools

- **Area:** `conf/profile_BFD.config`, BFD_FUNCTIONAL.
- **Consider:** run several genomes per job for `cazy`, `merops`, `wolfpsort`, `swissprot` and `predgpi`.
- **Why:** in the 2026-10-21 run (99 genomes) the median runtimes were 4 min (CAZy), 7 min (MEROPS), 4.5 min (WoLF PSORT), 1.7 min (SwissProt) and 2.5 min (predGPI). At about 23,500 genomes that is about 23,500 jobs per tool. The global guidance is about 1 to 1.5 h of work per job.
- **Check first:** none of these tools use a GPU (GPU labels: `signalp`, `idp`, and `deeptmhmm` only when `deeptmhmm_gpu` is true).
- **Status:** open. No design yet.

## 2026-10-21 — Queue settings need review after a full-size run

- **Area:** `conf/profile_BFD.config`.
- **Now:** `tmhmm`, `targetp`, `busco_pep` and `busco_genome` (attempt 1) use `preempt` with `--account=preempt`. `swissprot` uses `stajichlab`.
- **Consider:** the `preempt` account has a cpu=384 and mem=1T cap (see the 2026-10-08 rebalance note). About 23,500 genomes per tool may hit it. Check the pending reason after the first full run.
- **Consider:** the `stajichlab` partition has 3 abu_dhabi nodes (c01-c03) that lack BMI2, and 2 milan nodes. I did not test `diamond` (SwissProt) on the abu_dhabi nodes. If it crashes with SIGILL, add `--constraint=milan` to `swissprot`.
- **Consider:** one `TARGETP` task took 179 min against a 4 h limit. Preemption adds restarts, so raise the limit or split the work if tasks approach it.
- **Status:** open. Needs real timings from a full run.

## 2026-10-21 — Gzip audit for inputs read from `input/` (fixed; open items below)

- **Context:** `predict_results` files are compressed (`.proteins.fa.gz`, `.gff3.gz`, `.cds-transcripts.fa.gz`). `.scaffolds.fa` is never compressed. All 23,448 entries in `input/pep` and `input/gff3` are `.gz` links.
- **Fixed 2026-10-21:**
  - `bin/bfd_common.py` has `gff_stem()` and `is_gff_name()`, which accept `X.gff3` and `X.gff3.gz` and return `X`. `build_genestats_table.py` and `calculate_intergenic.py` use them (a plain `replace(gffext, "")` gave `X..gz`). Tested on a real genome: the `.gz` and plain GFF give identical intergenic output.
  - `workflows/paralogoscope.nf` uses `lookupFile()` for cds and gff3. `WGD_DMD`, `WGD_KSD` and `WGD_SYN` inflate a `.gz` input to its plain name first, so output names are unchanged (`merge_wgd_ks.py` parses `.cds-transcripts.fa`). I checked the rendered `WGD_DMD` script. I did not run `wgd` itself.
  - `subworkflows/local/PHYLING_ALIGN.nf` falls back to `<name>.gz`.
  - `STAGE_FILES` links a plain file and inflates a `.gz` to the plain `.faa` / `.cds.fa` name. Tested with plain, `.gz` and missing inputs.
  - `SETUP_SYMLINKS` links the tRNA GFF3 as plain or `.gz`, and removes stale links. Tested on a fake tree with plain, `.gz`, no-tRNA and stale-link cases.
- **Open:** `resolveGenomeFile()` in `modules/common/utils.nf` has no `.gz` fallback for `<id>.scaffolds.fa`. That is fine while scaffolds stay uncompressed.
- **Consider:** other scripts that strip file-name suffixes with `replace()` or `removesuffix()` may have the same bug. I only checked the ones that read GFF files.
- **Status:** mostly done.

## 2026-11-01 — Move `BUSCO_GENOME` from the `busco` module to a pinned container

- **Context:** `BUSCO_GENOME` failed for two genomes (`GCA_051941305.1_KEW_Xawe201`, strain `LIQ201XMSP`; `GCA_051941665.1_KEW_Poca`, strain `LIQ167XMSP`). BBTools 39.8x+ (`javasetup.sh`, `parseJavaArgs`) treats any `stats.sh` argument that contains `xms` or `xmx` (case-insensitive, anywhere in the string) as a Java heap flag. BUSCO passes the absolute input path to `stats.sh`, so the strain ID in the file name broke the JVM launch (`Invalid maximum heap size: -Xmxin=/bigdata/`).
- **Fixed 2026-11-01:** `modules/BFD/BUSCO_GENOME/main.nf` symlinks the staged genome to `genome.fa` (or `genome.fa.gz` when the input is compressed) and passes that name to `busco`. The work-dir hash is hex, so the rest of the path cannot contain `xms`/`xmx`. Verified by hand on the Gallowayella assembly: the bbtools step passed and miniprot started. Keep this rename after any container switch; every BUSCO 6.1.0 build tested has the bug.
- **Containers tested (same bug in all):** module `busco/6.1.0` (bbmap 39.81), `ezlabgva-busco-v6.1.0_cv3` (bbmap 40.02), biocontainers `busco:6.1.0--pyhdfd78af_1` (bbmap 39.81). Modules `busco/6.0.0` (bbmap 39.37) and `busco/5.8.0` (bbmap 39.10) use the old parser and are not affected. Report upstream to BBTools and BUSCO.
- **Plan (next wave boundary, never under a live run):**
  - Pin an image by digest and pre-pull it into `params.singularity_cache` as a `.sif`, the same way `aaftf_sif` is handled. Do not let Nextflow convert OCI to SIF on demand; that race corrupted the shared cache before (see `nextflow.config`). Candidate: `docker://ezlabgva/busco:v6.1.0_cv3` or the biocontainers build; both are already in the cache as `.img`.
  - Add a `busco_sif` param and switch `BUSCO_GENOME` (and `BUSCO_PEP`) to `apptainer exec` with the `SING_BINDS` pattern used in `modules/comparative/wgd/*`.
  - **Bind `/srv/projects/db/BUSCO` into the container.** The current bind list covers `/bigdata`, `/rhome`, `$PWD` and `$TMPDIR` only. Check that `/srv/projects/db` is mounted on every node class the `preempt` and `stajichlab` queues can land on.
  - Before adopting: run module vs container on about 20 genomes from `samples_test20.csv`. Compare wall time and the short-summary scores, not just exit codes. The official image carries a newer bbmap, so `stats.sh` output may differ.
- **Why:** reproducibility (pinned digest, citable), portability (NRP and other clusters), and no exposure to module-env changes mid-wave (the 6.0.0 to 6.1.0 bbmap jump is what introduced this bug).
- **Status:** open. No timing or score comparison has been run.
