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
