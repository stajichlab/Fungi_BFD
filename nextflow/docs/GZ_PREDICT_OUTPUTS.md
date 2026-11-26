# Gzipped predict/update/annotate results

Status 2026-10-07: branch `gz-predict-outputs`, not committed. The PREDICT part and the UPDATE/ANNOTATE
part are tested to different degrees. Read "Tested" and "NOT tested" before relying on it.

## What changed

- `FUNANNOTATE_PREDICT` gzips these files in `predict_results/` at the end (`pigz -f`):
  `.gbk`, `.gff3`, `.cds-transcripts.fa`, `.proteins.fa`, `.tbl`, `.txt`, `.mrna-transcripts.fa`.
  `.scaffolds.fa` stays plain.
- `SETUP_SYMLINKS` links a compressed result as `<name>.gz` when the plain file is missing.
- `SIGNALP_RUN` and `INTERPROSCAN_RUN` inflate `proteins.fa.gz` into `$SCRATCH` when only the `.gz` exists.
- `FUNANNOTATE_UPDATE` and `FUNANNOTATE_ANNOTATE` use the shell helpers in `modules/funannotate/utils.nf`
  (`gzShellHelpers()`):
  - Before the tool runs, `gz_inflate` writes plain copies next to the `.gz` (mtime copied from the `.gz`).
  - The `.gz` stays as the master copy. After the tool exits, `gz_cleanup` removes only the copies it made.
  - After a good run, `gz_compress` gzips the new `update_results/` or `annotate_results/` files.
  - UPDATE inflates `.gbk .gff3 .scaffolds.fa` in `predict_results/`.
  - ANNOTATE inflates `.gbk .gff3 .tbl` in `update_results/` if that dir exists, else `predict_results/`
    (this follows `annotate.py`).
- The empty-`update_results` cleanup in UPDATE now treats `.gbk.gz` as a valid result.

Why: funannotate `update.py` and `annotate.py` choose input files with `file.endswith(".gbk")` and similar,
so they do not see `.gz` files. See `~/projects/funannotate/funannotate-live/TODO_GZ_SUPPORT.md`.

## Tested

- Wave 1 and 2 do not run UPDATE or ANNOTATE (`run_update: false`, `run_annotate: false`).
- Stub run: a `.gbk.gz` result is treated as done; with it removed, PREDICT is scheduled.
- The `find ... pigz -f` command compresses a real result set, round-trips the gbk, and overwrites a stale `.gz`.
- The shell helpers, under `set -ue`: inflate round-trips, mtime is kept, cleanup removes only the inflated copies,
  an empty list does not fail.
- `nextflow lint` reports no errors for the four changed funannotate modules.

## NOT tested (needs testing before UPDATE or ANNOTATE is enabled)

- A real `funannotate update` and `funannotate annotate` run on a gzipped genome folder.
- A process killed mid-run (time limit, preemption): the inflated plain copies then stay next to the `.gz`.
  They hold the same content, and the plain file is preferred by `gbkResult()`. Delete them if disk matters.
- Two jobs on the same genome folder at once (in-place inflate would race). One job per folder is assumed.
- The stale-result checks (`-nt` against the GBK) after an inflate. The mtime is copied from the `.gz` for this reason.
- `FUNANNOTATE_ANNOTATION.nf` and `ANTISMASH_RUN` with a gzipped `update_results` or `annotate_results`.
  `ANTISMASH_RUN` uses `zcat` per the earlier check; not re-run.
- `.scaffolds.fa` is not compressed; `update` needs it plain.

## Still reading plain files

`workflows/BFD.nf` (`cdsIndex`), `STAGE_FILES`, `PHYLING_ALIGN`, `paralogoscope.nf`, `bin/merge_wgd_ks.py`,
and the `k8` and `development` InterProScan scripts read plain `input/cds`, `input/gff`, `input/pep`.
Fix them (read `.gz` or inflate into `$SCRATCH`) before running the comparative steps on gzipped output.

## genemark_store GTF compression (2026-11-24)

`genemark_store/<out>/<out>.genemark.gtf` may be stored as `<out>.genemark.gtf.gz`.

- Reader: `storedGenemarkFor()` (utils.nf) accepts either form. If both exist, the newer file wins.
- Predict: `FUNANNOTATE_PREDICT` inflates a `.gz` into `$TMPDIR` and passes that copy to `--genemark_gtf`. It deletes the copy after predict. The store keeps only the `.gz`.
- Writers (`GENEMARK_RUN`, `GENEMARK_TRAIN`, `GENEMARK_MERGE`) compress before publish. An `afterScript` makes `<out>.genemark.gtf.gz` at the end of the task. `publishDir` `saveAs` drops the plain GTF, so the store gets only the `.gz`. Steps in the same run still receive the plain GTF. A new output `gtf_gz` is declared (optional; no consumers). An empty GTF also gets a `.gz`.
- Existing plain GTFs: run `scripts/one-off/compress_genemark_store.py --store <genemark_store> --apply`. It is dry-run by default, verifies each `.gz` by sha256, and keeps the original mtime. Add `--remove-superseded` to delete an older plain GTF that sits next to a newer verified `.gz` (a pipeline rerun leaves these).
- Sizes: `scripts/one-off/genemark_store_census.py` (read-only; output in `analysis/genemark_store_census/`).
- `.mod` and `.other.gff3` stay plain.
- Tested 2026-11-24/25: reader probe (gz only, plain only, both, stale genome); converter on 3 real folders (10.7x); inflate block in bash; a real `funannotate predict` on Saccharomyces_x_bayanus_NBRC1948 from a `.gz` store entry (GeneMark contributed 4,889 models; GBK gene count 5,144, same as the earlier plain-GTF run); stub runs of GENEMARK_RUN, GENEMARK_TRAIN and GENEMARK_MERGE (store gets the `.gz` only). NOT tested: a real (non-stub) GeneMark run through the new `afterScript`, Prodigal/other.gff3 genomes, RNA-seq genomes, the chunked route end to end.
