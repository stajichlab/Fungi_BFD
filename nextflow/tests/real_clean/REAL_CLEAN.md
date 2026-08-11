# REAL CLEAN — AAFTF-SIF genome-cleaning validation

End-to-end validation of the AAFTF v0.7.0-beta.2 **SIF** integration in the
genome-cleaning step (`GENOME_CLEAN` / `GENOME_CLEAN_BATCH` / `CALC_ASM_STATS`),
replacing the old `module load AAFTF` / `module load taxonkit`.

**Status**: PASSED (2026-08-10). The SIF wiring is the only active path in all
three module files; the module-load form is retired.

## Objective

Prove the in-script `singularity exec` form works on real data end-to-end:
FCS-GX DB staged into host `/dev/shm`, `AAFTF fcs_gx_purge` + taxonkit phylum
resolution running inside the SIF against the cached taxdump, and real cleaned
outputs produced — before declaring the module replacement permanent.

## Test inputs (`samples.real_clean.csv`, `source/`)

3 small real assemblies symlinked from the 1KFG tree (Ascomycota: S. cerevisiae
`GCF_000146045.2_R64`, S. pombe `GCF_000002945.2_ASM294v3`; Basidiomycota:
C. neoformans `GCA_056632885.1_ASM5663288v1`), chosen to cover both divisions
and to keep the full clean run short.

## Configuration (`conf/test_clean_real.config`)

- samples/source overridden to `tests/real_clean/samples.real_clean.csv` + `source/`
- `clean_batch_size = 3` (one batch)
- `only_clean = true` (run only SETUP_TAXONDB + GENOME_CLEAN_BATCH; RNASEQ/PREDICT
  stay `[-]` — verified emptied downstream channel)
- `run_ani_reuse = false` — REQUIRED: `workflows/funannotate.nf` hard-errors at
  build time (`WorkflowScriptErrorException`, exit 1) if `run_ani_reuse=true`
  while `params.abinitio_reuse_csv` is unset / not relevant for `only_clean`.

## Launch

Interactive (proven pattern — do NOT launch real runs via sbatch, see below):

```bash
cd nextflow
screen -dmS real_clean bash -lc 'export NXF_OPTS="-Xms512m -Xmx4g"; \
  nextflow run main.nf -c nextflow.config -c conf/test_clean_real.config \
  -profile funannotate --pipeline funannotate -resume \
  > logs/nextflow/real_clean_driver.log 2>&1'
tail -25 logs/nextflow/real_clean_driver.log
```

`run_real_clean.sh` (sbatch wrapper) is kept for reference but is known-broken:
nextflow launched from an sbatch job fails deterministically on
`.nextflow/history.lock` (NoSuchFileException) and top-level `mkdir` (Permission
denied) even though interactive/`screen` shells work from the same directory.
See `.living/learnings.md` (2026-08-10).

## Results (job 27332111, h04, 18:05 runtime, pipeline `Succeeded: 2`)

| Check | Result |
|-------|--------|
| SETUP_TAXONDB | storeDir short-circuit; real taxdump at `lib/taxdump` (533 MB: `names.dmp` 303 MB, `nodes.dmp` 224 MB) |
| FCS-GX staging | 498.6 GB → host `/dev/shm` in 914 s (~520 MB/s), logged `logs/nextflow/fcs_gx_shm_timing.tsv`; freed after job (`FreeMem` back to ~1003 GB) |
| Cleaned genomes | `input_clean_genomes/<asmid>.fa.gz` ×3 + `clean/<asmid>.purge.fasta.gz` + `*.purge.fcs_gx-taxonomy.tsv.gz` ×3 |
| Taxonomy | S288C `fung:ascomycetes` + `budding yeasts`, agg-cov 0.998; C. neoformans `fung:basidiomycetes`, agg-cov 0.979 (FCS DB 2023-01-24 build, 3.0M seqs / 709 Gbp) |
| Suppression | `TO_ADD_TO_SUPRESS.csv` 0 B (nothing below min-assembly length; expected for these yeasts) |
| Downstream | RNASEQ/PREDICT correctly skipped (`only_clean`) |

Notes: `clean_batch_1.manifest.tsv` is written in the task work dir and scrubbed
on success because `profile_funannotate` sets `cleanup = true` — expect it
absent after a clean run (no downstream consumer in `only_clean`).

## Caveats observed

- Batch demands 16 cpus + **500 GB** pinned to `-w h04,h05,h06`; with those
  nodes already ~70–80% memory-allocated by other jobs the task sat PENDING
  (`Resources`) ~6 h before h04 freed enough memory. Expect an unbounded queue
  wait on a busy cluster.
- `lib/taxdump` must be free of 0-byte stub `.dmp` files before a real run
  (stub runs create empty dump files that a real run would treat as cached →
  silent taxonkit failure). Cleared pre-run.

## Reproduce

`screen` launch above from `nextflow/` (launchDir must stay `nextflow/` for the
output paths and `params.taxondb` → `lib/taxdump` resolution to match
production). `-resume` reuses SETUP_TAXONDB; re-trigger the batch after the
current `input_clean_genomes/*.fa.gz` exist can be skipped via the existing
storeDir-style skip, so delete `input_clean_genomes/` first to force a true
fresh clean.
