# Funannotate Pipeline (`funannotate.nf`)

Nextflow DSL2 pipeline for fungal genome prediction and annotation on SLURM HPC.
Covers the full funannotate workflow: genome cleaning → repeat masking → RNA-seq
discovery → download → train → predict → (optional) annotate.

Source synced from `../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf`.

---

## Quick start

```bash
# Default run: clean + mask + SRA query + fetch + train + predict (no annotate)
sbatch nextflow/run_funannotate.sh

# Survey RNA-seq availability only (produces samples.rnaseq_sra.csv, no downloads)
sbatch nextflow/run_funannotate.sh --stop_after_sra_query true

# Stop after SRA download/normalisation (skip train/predict)
sbatch nextflow/run_funannotate.sh --stop_after_sra_fetch true

# Predict + annotate
sbatch nextflow/run_funannotate.sh --run_annotate true

# Genome cleaning only
sbatch nextflow/run_funannotate.sh --only_clean true

# Restrict to a taxonomic group
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Ascomycota
sbatch nextflow/run_funannotate.sh --taxon CLASS:Sordariomycetes --n_test 5

# Resume a stopped run
sbatch nextflow/run_funannotate.sh -resume
```

---

## Pipeline stages

```
samples.csv
    │
    ▼
SETUP_TAXONDB         Download NCBI taxdump (storeDir; once per deployment)
    │
    ▼
GENOME_CLEAN          Decompress genome, run NCBI FCS-GX contamination purge,
    │                 filter short contigs → input_clean_genomes/<asmid>.fa.gz
    │                 (gzip-compressed to save space; storeDir; skipped if a clean
    │                  .fa.gz — or legacy uncompressed .fa — already exists)
    ▼
MASKREPEAT_TANTAN_RUN Soft-mask repeats with tantan via funannotate mask
    │                 → input_clean_genomes/<asmid>.masked.fasta.gz (gzip-compressed)
    │                 (storeDir; skipped if masked file exists; skippable with
    │                  --run_repeatmasker false, falls back to clean .fa[.gz])
    ▼
SRA_QUERY_BATCH       Batched NCBI SRA survey — up to sra_query_batch_size species
    │                 per SLURM job (default 100), max 10 concurrent jobs.
    │                 Runs esearch/efetch runinfo query per species; no downloading.
    │                 Species already cached in rnaseq_reads/sra_query/ are reused
    │                 (file copy, no NCBI call); only uncached species are queried.
    │                 Primary query: PAIRED[Layout]; records up to 5 candidates.
    │                 SE fallback: if no PE hits found AND --enable_single_end true,
    │                 a second SINGLE[Layout] query runs and records up to
    │                 --max_rnaseq_se_runs candidates.
    │                 Output: rnaseq_reads/sra_query/<species_tag>.sra_query.csv
    │                 Columns: species_tag, taxonid, sra_accession, spots, platform, layout
    │                 (publishDir copy, overwrite: false; delete CSV to re-query)
    │                 (only runs when --run_sra_fetch true)
    ▼
COLLECT_SRA_QUERY     Merges all per-species query CSVs into a single named manifest
    │                 written alongside the input samples file:
    │                   <stem>.rnaseq_sra.csv   (e.g. samples.rnaseq_sra.csv)
    │                 Columns: species_tag, taxonid, sra_accession, spots, platform, layout
    │                 ── stop here with --stop_after_sra_query true ──
    ▼
    ├── (species with PE accessions not overridden to SE)
    │       ▼
    │   SRA_FETCH         Download up to --max_rnaseq_runs PE accessions (read from
    │       │             the per-species query CSV — no NCBI call at this stage).
    │       │             Per-accession: parallel-fastq-dump → enforce read-pair
    │       │             length (EBI FTP fallback on mismatch) → bbnorm (target=30,
    │       │             ecc=t) → fastp QC trim.
    │       │             Writes rnaseq_reads/<species_tag>_norm_{R1,R2}.fastq.gz
    │       │             and a zero-byte _norm_SE.fastq.gz stub.
    │       │             (storeDir; skipped if all three files already exist)
    │       │             Accessions with blacklist action SE_trinity are skipped
    │       │             here — they route to SRA_FETCH_SE instead.
    │
    ├── (species with SE_trinity blacklist overrides OR SINGLE-layout hits
    │    when --enable_single_end true; PE always takes priority)
    │       ▼
    │   SRA_FETCH_SE      Download up to --max_rnaseq_se_runs SE accessions.
    │       │             SE_trinity: pfd --split-files, take _1 only (real SE data).
    │       │             SINGLE layout: pfd gives ACC.fastq.gz or ACC_1.fastq.gz.
    │       │             EBI FTP fallback on pfd failure.
    │       │             Per-accession: fix_fastq_header_trinity → bbnorm (target=30,
    │       │             ecc=f) → fastp QC trim.
    │       │             Writes zero-byte _norm_{R1,R2}.fastq.gz PE stubs and
    │       │             rnaseq_reads/<species_tag>_norm_SE.fastq.gz.
    │       │             (storeDir; skipped if all three files already exist)
    │
    └── (species with no usable SRA hits)
            ▼
        WRITE_EMPTY_READS  Writes zero-byte placeholder files without a SLURM
                          download job:
                          rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz
                          (storeDir; skipped if placeholder files already exist)
    │
    │  ── stop here with --stop_after_sra_fetch true ──
    ▼
RNASEQ_PREPARE        Run funannotate train --stop_after_trinity on the
    │                 representative (first) assembly per species.
    │                 PE mode:  --left_norm / --right_norm + --jaccard_clip
    │                 SE mode:  --single_norm (no --jaccard_clip)
    │                 Archives Trinity-GG FASTA to rnaseq_data/<species_tag>.trinity-GG.fasta
    │                 (storeDir). All other strains of the same species reuse this.
    ▼
FUNANNOTATE_TRAIN     Run funannotate train on every assembly (4 branches):
    │                 - Shared Trinity + PE reads: PASA only with --left_norm/--right_norm
    │                 - Shared Trinity + SE reads: PASA only with --single_norm
    │                 - No Trinity + PE reads:    full train with --left_norm/--right_norm
    │                 - No Trinity + SE reads:    full train with --single_norm
    │                 Removes large intermediates (hisat2/, trinity_gg/) after completion.
    │                 Skipped at channel level if funannotate_train.pasa.gff3 already exists.
    ▼
FUNANNOTATE_PREDICT   Run funannotate predict (Augustus + EvidenceModeler).
    │                 Uses training data linked from genome_annotation/<out>/training/.
    │                 Skipped if predict_results/<out>.gbk (or .gbk.gz) already exists.
    │                 Output: genome_annotation/<out>/predict_results/
    ▼
  (optional post-predict steps, run in parallel)
  ├── ANTISMASH_RUN       antiSMASH BGC prediction
  ├── INTERPROSCAN_RUN    InterProScan 5 XML
  └── SIGNALP_RUN         SignalP 6 (GPU)
    │
    ▼
FUNANNOTATE_ANNOTATE  funannotate annotate — functional annotation, GO terms,
                      BUSCO evidence; optionally incorporates antiSMASH and
                      InterProScan results.
                      Output: genome_annotation/<out>/annotate_results/
```

---

## RNA-seq discovery and download: two-phase design

The SRA workflow is split into a **lightweight query phase** (SRA_QUERY_BATCH) and a
**heavy download phase** (SRA_FETCH / SRA_FETCH_SE). This avoids allocating large SLURM
jobs for species that have no RNA-seq data, and makes the discovery results inspectable
before committing to downloads.

### Phase 1 — SRA_QUERY_BATCH (survey)

- Groups up to `--sra_query_batch_size` unique species (default 100) per SLURM job;
  at most 10 batch jobs run concurrently (`maxForks 10`).
- Per species within each batch: sends a single `esearch | efetch -format runinfo`
  call to NCBI SRA. Species already cached in `rnaseq_reads/sra_query/` are reused
  (file copy only — no NCBI call). Only uncached species hit the network.
- **Primary (PE) query filter criteria:**
  - `txid<taxonid>[Organism:noexp]` — exact taxon, no descendants
  - `RNA-Seq[Strategy]`, `PAIRED[Layout]`
  - Read length 75–300 bp
  - Illumina or BGI platform
  - ≥250,000 read pairs
  - Records **up to 5** candidates sorted by spot count descending.
- **SE fallback** (only when `--enable_single_end true` **and** no PE hits found):
  Re-runs the query with `SINGLE[Layout]` (Illumina only) and records **up to
  `--max_rnaseq_se_runs`** (default 3) candidates with `layout=SINGLE`.
- Output CSV columns: `species_tag, taxonid, sra_accession, spots, platform, layout`
- Results published to `rnaseq_reads/sra_query/<species_tag>.sra_query.csv`
  (`publishDir overwrite: false`). Delete the per-species CSV to force re-query.
- Resources: 1 CPU, 4 GB RAM, 4h per batch job (short queue); up to 2 retries.
- **Backward compatibility:** old cached CSVs with 5 columns (no `layout`) are handled
  gracefully — the pipeline treats the missing column as `PAIRED`.

### Phase 2 — COLLECT_SRA_QUERY (merge manifest)

- Concatenates all per-species CSVs into a single project-level manifest.
- Written alongside the input samples file as `<stem>.rnaseq_sra.csv`
  (e.g., `samples.rnaseq_sra.csv` when `--samples samples.csv`).
- Re-generated each run from the cached per-species CSVs — always reflects the
  current set of queried species.
- **Use `--stop_after_sra_query true`** to halt here and inspect the manifest
  before committing to downloads. Useful for large initial runs or when checking
  RNA-seq coverage for a new clade.

### Phase 3 — SRA_FETCH / SRA_FETCH_SE / WRITE_EMPTY_READS (download)

The channel performs a **three-way branch** based on the per-species CSV and the
`rnaseq_blacklist.csv` override file:

- **PE accessions present (not overridden) → SRA_FETCH**: downloads up to
  `--max_rnaseq_runs` paired-end accessions. Per-accession:
  `parallel-fastq-dump → enforce_seqpair_readlen (EBI FTP fallback on mismatch)
  → bbnorm (target=30, ecc=t) → fastp`. Writes `_norm_{R1,R2}.fastq.gz` plus a
  zero-byte `_norm_SE.fastq.gz` stub. PE always wins — if a species has PE
  accessions, SE is never fetched.
- **SE_trinity override OR SINGLE layout (if enabled) → SRA_FETCH_SE**: downloads
  up to `--max_rnaseq_se_runs` single-end reads. Per-accession:
  `parallel-fastq-dump → fix_fastq_header_trinity → bbnorm (target=30, ecc=f)
  → fastp`. Writes `_norm_SE.fastq.gz` plus zero-byte `_norm_{R1,R2}.fastq.gz` PE stubs.
- **No usable accessions → WRITE_EMPTY_READS**: writes zero-byte placeholder files
  for all three slots (`_norm_R1`, `_norm_R2`, `_norm_SE`) in a trivial job.

All three branches write to `rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz`
and use `storeDir`, so the result is cached identically regardless of path taken.

---

## Per-accession overrides: `rnaseq_blacklist.csv`

The file `rnaseq_blacklist.csv` (project root) allows per-accession overrides that
change how SRA_FETCH or SRA_FETCH_SE handles a specific run.

**Format:** CSV with a header row; column 1 is the SRA accession, column 4 is the action.

```
accession,taxonid,species_tag,action
SRR17583173,,Hemileia_vastatrix,SE_trinity
ERR123456,,Botrytis_cinerea,skip
SRR654321,,Neurospora_crassa,rename_headers
```

| Action | Effect |
|---|---|
| `skip` | Exclude the accession entirely from all processing |
| `rename_headers` | Download normally; replace FASTQ headers with sequential integers before pairing (fixes BGI/MGISEQ block-splitting desync; applied automatically for BGISEQ platform) |
| `SE_trinity` | Treat this accession as single-end despite SRA metadata saying PAIRED. SRA_FETCH skips it; SRA_FETCH_SE downloads `_1.fastq.gz` only and discards any `_2`. **Bypasses `--enable_single_end`** — the entry takes effect regardless of that flag. Useful for mislabeled runs where pfd produces mismatched R1/R2. |

**Routing priority:**
1. A species with at least one usable PE accession (not `skip` / not `SE_trinity`) → **SRA_FETCH**.
   Any `SE_trinity` entries for that species are silently skipped — PE wins.
2. A species where *all* PE accessions are `skip`/`SE_trinity` → **SRA_FETCH_SE** (always,
   regardless of `--enable_single_end`).
3. A species with only `SINGLE`-layout accessions and `--enable_single_end true` → **SRA_FETCH_SE**.
4. No usable accessions → **WRITE_EMPTY_READS**.

---

## SRA query cache

### How it is built

On the first run with `--run_sra_fetch true`, `SRA_QUERY_BATCH` runs for every
unique `species_tag` in the filtered sample set. Species are collated into batches of
`sra_query_batch_size` (default 100); each batch is a single SLURM job. Within each
batch, uncached species are queried sequentially and results written to:

```
rnaseq_reads/sra_query/<species_tag>.sra_query.csv
```

The header row is always written. If no qualifying accessions are found, the file
contains only the header (0 data rows). Failed queries (network error, timeout) are
retried up to 3 times with exponential backoff before an empty CSV is written.

### How it is used

On subsequent runs, `SRA_QUERY_BATCH` still runs for each batch, but species whose
per-species CSV already exists in `rnaseq_reads/sra_query/` are served from the cache
(file copy, no NCBI call). Only species without a cached CSV query NCBI.

The channel branch that drives `SRA_FETCH` vs `WRITE_EMPTY_READS` is evaluated from
the cached CSV at channel-evaluation time (`.readLines().size() > 1`). If `SRA_FETCH`
has also already completed for a species (its `_norm_R1/R2.fastq.gz` files exist in
`rnaseq_reads/`), it too is skipped via its own `storeDir` check.

### Cache invalidation

| Scope | Action |
|---|---|
| Single species | `rm rnaseq_reads/sra_query/<species_tag>.sra_query.csv` then re-run |
| All species | `rm -r rnaseq_reads/sra_query/` then re-run |
| Force re-download (keep query cache) | `rm rnaseq_reads/<species_tag>_norm_*.fastq.gz` then re-run |
| Full SRA reset | Remove both `rnaseq_reads/sra_query/` and `rnaseq_reads/<tag>_norm_*.fastq.gz` |
| Re-route PE species to SE (SE_trinity) | Add entry to `rnaseq_blacklist.csv`, then `rm rnaseq_reads/<tag>_norm_*.fastq.gz` and re-run. storeDir will block re-fetch unless all three `_norm_*` files are removed. |
| Enable SE for PE-absent species | Set `--enable_single_end true`, delete the per-species sra_query CSV (to trigger SE fallback query), then re-run. |

> **Note:** `rnaseq_reads/` is under `cleanup = true` scope in the Nextflow work
> directory, but the actual output files live outside the work directory (storeDir),
> so they are **not** removed by Nextflow cleanup. Delete them manually when needed.

---

## Output structure

```
genome_annotation/
  <Species_Strain>/
    predict_results/         ← primary output used by BFD.nf
      <name>.gbk             GBK (may be compressed to <name>.gbk.gz to save space;
                             both forms count as complete — see Skip/cache behavior)
      <name>.proteins.fa
      <name>.cds-transcripts.fa
      <name>.gff3
      <name>.scaffolds.fa
    predict_misc/
      trnascan.no-overlaps.gff3
      ab_initio_parameters/
    training/
      funannotate_train.pasa.gff3
      *.bam, *.bai
    annotate_results/        ← produced only if --run_annotate true
    update_results/          ← produced only if --run_update true
    antismash_local/         ← produced only if --run_antismash true

input_clean_genomes/
  <asmid>.fa                 FCS-GX cleaned genome
  <asmid>.masked.fasta       tantan soft-masked genome
  clean/                     FCS-GX intermediates (.purge.fasta.gz, .fcs_gx-taxonomy.tsv.gz)

rnaseq_reads/
  sra_query/
    <species_tag>.sra_query.csv   per-species SRA survey results (cached)
                                  columns: species_tag, taxonid, sra_accession, spots, platform, layout
  <species_tag>_norm_R1.fastq.gz  PE: bbnorm+fastp normalized R1 (or 0-byte stub if SE/no data)
  <species_tag>_norm_R2.fastq.gz  PE: bbnorm+fastp normalized R2 (or 0-byte stub if SE/no data)
  <species_tag>_norm_SE.fastq.gz  SE: bbnorm+fastp normalized SE reads (or 0-byte stub if PE/no data)

<samples_stem>.rnaseq_sra.csv    merged SRA survey manifest (alongside samples.csv)
                                  columns: species_tag, taxonid, sra_accession, spots, platform, layout

rnaseq_data/
  <species_tag>.trinity-GG.fasta  shared Trinity assembly per species

logs/nextflow/
  funannotate_trace.txt
  funannotate_report.html
  funannotate_timeline.html
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` | Master species/genome table |
| `--target` | `genome_annotation/` | Output root for funannotate folders |
| `--source` | `/bigdata/.../NCBI_ASM` | Directory of `<asmid>/<asmid>_genomic.fna.gz` inputs |
| `--taxon` | `""` (all) | Restrict to `RANK:VALUE`, e.g. `PHYLUM:Ascomycota` |
| `--n_test` | `0` (all) | Limit to first N samples after taxon filter |
| `--suppress` | `""` | File of ASMIDs to skip (one per line, first comma-delimited field) |
| `--only_clean` | `false` | Stop after GENOME_CLEAN |
| `--run_repeatmasker` | `true` | Run tantan soft-masking |
| `--run_sra_fetch` | `true` | Run SRA_QUERY_BATCH + SRA_FETCH/SRA_FETCH_SE (both phases) |
| `--stop_after_sra_query` | `false` | Halt after COLLECT_SRA_QUERY; produces `<stem>.rnaseq_sra.csv` only |
| `--stop_after_sra_fetch` | `false` | Halt after SRA_FETCH/SRA_FETCH_SE; skip train/predict |
| `--sra_query_batch_size` | `100` | Species per SRA_QUERY_BATCH job; reduce if batches time out |
| `--max_rnaseq_runs` | `2` | Max PE accessions downloaded per species (SRA_FETCH; SRA_QUERY_BATCH surveys up to 5) |
| `--enable_single_end` | `false` | Enable SE fallback: query SINGLE-layout SRA accessions when no PE found; required for SINGLE-layout species to be fetched. SE_trinity blacklist overrides bypass this flag. |
| `--max_rnaseq_se_runs` | `3` | Max SE accessions downloaded per species (SRA_FETCH_SE) |
| `--run_annotate` | `false` | Run funannotate annotate after predict |
| `--run_update` | `false` | Run funannotate update (requires `--run_sra_fetch true`) |
| `--run_antismash` | `false` | Run antiSMASH before annotate |
| `--run_interpro` | `false` | Run InterProScan 5 before annotate |
| `--run_signalp` | `false` | Run SignalP 6 GPU before annotate |
| `--pasa_mysql` | `true` | Use per-task MariaDB for PASA (recommended) |
| `--max_intronlen` | `3000` | Max intron length passed to funannotate |
| `--min_intronlen` | `10` | Min intron length passed to funannotate |
| `--min_contig_len` | `2000` | Contigs shorter than this are dropped after FCS-GX |
| `--debug` | `false` | Print extra channel and task diagnostics |

---

## SLURM resource allocation

| Process | Queue | CPUs | Memory | Notes |
|---|---|---|---|---|
| `SETUP_TAXONDB` | short | 1 | 1 GB | Once per deployment; storeDir-cached |
| `GENOME_CLEAN` | highmem | 16 | 500 GB | FCS-GX loads /dev/shm/gxdb on h04/h05/h06 |
| `MASKREPEAT_TANTAN_RUN` | short | 2 | 16 GB | |
| `SRA_QUERY_BATCH` | short | 1 | 4 GB | Batched esearch/efetch; up to 10 concurrent; up to 2 retries per batch |
| `COLLECT_SRA_QUERY` | short | 1 | 1 GB | Merges per-species CSVs; fast |
| `WRITE_EMPTY_READS` | short | 1 | 1 GB | Zero-byte placeholder; no-data species only |
| `SRA_FETCH` | short → epyc | 24 → 32 | 48 → 192 GB | PE species only; retry bumps memory and queue |
| `SRA_FETCH_SE` | short → epyc | 24 → 32 | 48 → 192 GB | SE species only (SE_trinity or SINGLE layout); same retry profile as SRA_FETCH |
| `RNASEQ_PREPARE` | epyc | 16 | 96 GB (+48 GB/retry) | Trinity-GG assembly; up to 3 retries |
| `FUNANNOTATE_TRAIN` | epyc | 8 → 24 | 96 → 192 GB | PASA alignment; up to 3 retries |
| `FUNANNOTATE_PREDICT` | epyc | 16 | 32 GB | |
| `SIGNALP_RUN` | exfab | 16 | 64 GB | Requires `--gres=gpu:1` |
| `INTERPROSCAN_RUN` | epyc | 8 | 32 GB | |
| `ANTISMASH_RUN` | epyc | 8 | 24 GB | |
| `FUNANNOTATE_ANNOTATE` | epyc | 16 | 64 GB | |

---

## Skip/cache behavior

All expensive one-time steps use `storeDir`, which means Nextflow skips the process
if all declared output files already exist on disk — even without `-resume`:

| Step | Cache location | Skip condition |
|---|---|---|
| `SETUP_TAXONDB` | `taxondb/` | `taxondb/names.dmp` exists |
| `GENOME_CLEAN` | `input_clean_genomes/` | `<asmid>.fa.gz` (or legacy `<asmid>.fa`) exists |
| `MASKREPEAT_TANTAN_RUN` | `input_clean_genomes/` | `<asmid>.masked.fasta.gz` (or legacy `<asmid>.masked.fasta`) exists |
| `SRA_QUERY_BATCH` | `rnaseq_reads/sra_query/` | per-species CSV exists (checked inside the batch script; `publishDir overwrite: false`) |
| `WRITE_EMPTY_READS` | `rnaseq_reads/` | all three `<tag>_norm_{R1,R2,SE}.fastq.gz` exist (0-byte) |
| `SRA_FETCH` | `rnaseq_reads/` | all three `<tag>_norm_{R1,R2,SE}.fastq.gz` exist |
| `SRA_FETCH_SE` | `rnaseq_reads/` | all three `<tag>_norm_{R1,R2,SE}.fastq.gz` exist |
| `RNASEQ_PREPARE` | `rnaseq_data/` | `<tag>.trinity-GG.fasta` exists |

`FUNANNOTATE_TRAIN` and `FUNANNOTATE_PREDICT` skip via channel-level file-existence
checks (`funannotate_train.pasa.gff3` and `predict_results/<out>.gbk`), so they also
skip gracefully when re-running over partially completed datasets.

### Compressed input genomes (`.fa.gz` / `.masked.fasta.gz`)

To save space, cleaned and repeat-masked genomes in `input_clean_genomes/` are stored
gzip-compressed: `GENOME_CLEAN`/`GENOME_CLEAN_BATCH` emit `<asmid>.fa.gz` and
`MASKREPEAT_TANTAN_RUN` emits `<asmid>.masked.fasta.gz`. The genome path threaded
downstream (`genome_fa`) may therefore point at a `.gz` file.

- **Resolution:** all existence/skip gates use the `genomeFile()` helper, which accepts
  either the compressed (`.gz`) or a legacy uncompressed form, preferring the compressed
  one. So a folder compressed after a prior run is **not** re-cleaned/re-masked, and
  genomes left uncompressed from older runs still work.
- **Consumers decompress on the fly:** `MASKREPEAT_TANTAN_RUN`, `RNASEQ_PREPARE`,
  `FUNANNOTATE_TRAIN`, and `FUNANNOTATE_PREDICT` each inflate a `.gz` genome to a local
  `genome_input.fa` in the task work dir before calling `funannotate ... -i` (funannotate
  cannot read a gzipped FASTA via `-i`). The work-dir copy is discarded with the task.
- **EarlGrey pipeline:** the separate `nextflow/earlgrey_mask.nf` curated repeat-masking
  pipeline gets the same treatment — it resolves clean genomes via its own `genomeFile()`,
  inflates `.fa.gz` before EarlGrey/RepeatMasker, and now delivers its soft-masked output
  **compressed** as `input_clean_genomes/<asmid>.masked.fasta.gz` (which the tantan
  `MASKREPEAT_TANTAN_RUN` storeDir and `FUNANNOTATE_PREDICT` both accept). The
  `select_repeat_representatives.py` existence check also accepts `.fa.gz`.

### Compressed GenBank outputs (`.gbk.gz`)

To save space, a step's GenBank output may be stored gzip-compressed (`.gbk.gz`)
instead of plain `.gbk`. All completion/skip gates accept **either** form (via the
`gbkResult()` helper), so compressing a finished folder does **not** force a re-run.
This applies to predict (`predict_results/<out>.gbk`), update
(`update_results/<out>.gbk`), and annotate (`annotate_results/<out>.gbk`), plus the
in-script skip guards in `FUNANNOTATE_TRAIN`/`FUNANNOTATE_PREDICT` and the RNA-seq
staleness test below. `ANTISMASH_RUN` transparently inflates a gzipped predict GBK
to a local copy before running.

**Caveats before compressing:**
- `funannotate annotate -i <dir>` and `funannotate update -i <dir>` read
  `predict_results/<out>.gbk` *internally* and cannot consume a `.gbk.gz`. Only
  compress `predict_results/<out>.gbk` for genomes that are fully done through any
  annotate/update you intend to run. With the defaults
  (`run_annotate`/`run_update`/`run_antismash` all `false`) predict is the terminal
  artifact, so this is safe. `update_results` and `annotate_results` GBKs are terminal
  and always safe to compress.
- Skip/staleness checks compare RNA-seq mtimes against the GBK mtime, so preserve the
  original mtime when compressing (`gzip` keeps it by default; use `gzip -n` to be safe)
  to avoid spuriously triggering a stale re-predict.
- Post-run assertions (verifying funannotate just wrote output) still expect plain
  `.gbk`, since funannotate always writes uncompressed — compression is a separate
  housekeeping step applied after a genome completes.

### RNA-seq staleness detection (train + predict)

The pipeline also detects when RNA-seq data has been refreshed *after* a genome was
already predicted, and automatically re-runs training and prediction for those genomes.

**What "stale" means:** A prediction is considered stale when any of the following is
newer than `genome_annotation/<out>/predict_results/<out>.gbk` (or its `.gbk.gz`):
- `rnaseq_reads/<species_tag>_norm_R1.fastq.gz` (PE reads), **or**
- `rnaseq_reads/<species_tag>_norm_SE.fastq.gz` (SE reads), **or**
- `rnaseq_data/<species_tag>.trinity-GG.fasta` (shared Trinity assembly).

In both cases the reads used to build the current prediction are no longer the most
recent ones on disk.

**How the pipeline responds:**

1. At the channel level, any assembly whose prediction is stale is moved from the
   "training complete" bucket into the "needs training" bucket, even though
   `funannotate_train.pasa.gff3` already exists.  This forces `FUNANNOTATE_TRAIN`
   to run for that assembly.

2. Inside the `FUNANNOTATE_TRAIN` script, before the normal "training already
   complete; skipping" guard, the script checks whether the R1 file is newer than
   the existing GBK (`-nt` test).  If it is, the old `training/` directory is
   deleted so funannotate runs a clean re-train rather than silently re-using stale
   PASA output.

3. At the predict channel filter, the same staleness test overrides the normal
   "GBK exists; skip predict" guard, so `FUNANNOTATE_PREDICT` re-runs after the
   fresh training finishes.

This means a normal pipeline re-run (`nextflow run funannotate.nf -resume`) will
automatically pick up any new RNA-seq reads without any manual intervention.

`SRA_QUERY_BATCH` uses `publishDir` (not `storeDir`) — the SLURM job always runs,
but the per-species cache check inside the script means only uncached species incur
an NCBI query. `COLLECT_SRA_QUERY` also always re-merges the per-species CSV cache
into a fresh `<stem>.rnaseq_sra.csv`, so the manifest always reflects the current
run's species set.

---

## Required modules

```
# Genome cleaning
miniconda3, AAFTF, taxonkit

# Repeat masking + prediction + annotation
funannotate

# SRA discovery (SRA_QUERY_BATCH only)
ncbi_edirect

# SRA download + normalisation (SRA_FETCH and SRA_FETCH_SE)
sratoolkit, parallel-fastq-dump, fastp, BBTools, workspace/scratch
# EBI FTP fallback (loaded on demand inside SRA_FETCH / SRA_FETCH_SE)
aria2

# Optional post-predict annotation
signalp/6-gpu         — GPU signal peptide prediction
interproscan          — domain annotation
antismash             — BGC prediction
singularity           — MariaDB container for PASA (if pasa_mysql=true)
```

---

## Relationship to BFD.nf

`funannotate.nf` produces `genome_annotation/<Species_Strain>/predict_results/`.
`BFD.nf` reads from that directory via its `SETUP_INPUT` step, which creates
symlinks in `input/pep/`, `input/cds/`, `input/gff3/`, `input/dna/`, and `input/trna/`.
Run `funannotate.nf` first, then `BFD.nf`.

---

## Cleanup script: checking for missing or outdated training

`scripts/check_rnaseq_training.py` is a standalone audit tool that scans the project
and reports every `genome_annotation/` folder that needs a (re-)training run.

### When to use it

Run it before submitting a pipeline job to get a plain-English summary of what is
missing or stale, without actually executing anything:

```bash
python scripts/check_rnaseq_training.py --tsv
```

### What it checks

For every entry in `samples.csv` that has a **non-zero**
`rnaseq_reads/<species_tag>_norm_R1.fastq.gz` file *and* a corresponding
`genome_annotation/<out>/` folder, the script evaluates two conditions:

| Status | Meaning |
|---|---|
| `MISSING_TRAINING` | `genome_annotation/<out>/training/` directory does not exist — training has never run for this genome |
| `NEW_READS` | The `training/` directory exists but the R1 file is newer than the directory — RNA-seq reads were refreshed after training last ran |
| `OK` | Training exists and is up to date (shown only with `--report-ok`) |

### How species and folder names are matched

The script derives both the `<out>` folder name and the `<species_tag>` file prefix
using the same rules as the Nextflow pipeline:

The `STRAIN` field from `samples.csv` is sanitised before the name is assembled:

1. Leading/trailing whitespace stripped; quotes (`'` `"`) removed.
2. Only the first `;`-delimited token is kept (some entries list synonyms after a semicolon).
3. **Colons (`:`) replaced with a space** — e.g. `CBS:123` → `CBS 123`.

Then:

- `<out>` = `{SPECIES}_{STRAIN}` with spaces → `_` and special characters
  (`[ ] * ? { }`) → `_` (mirrors the `make_out` Groovy closure in `funannotate.nf`)
- `<species_tag>` = `{SPECIES}` with spaces → `_`

Matching is driven by `samples.csv` rather than directory-name prefix scanning, so
species like *Aspergillus sp.* and *Aspergillus sp. 2663* — which share a common
prefix but are distinct organisms — are never confused.

### Options

```
--project-dir DIR    Project root (default: current directory)
--tsv                Print a TSV header line before output
--report-ok          Also print rows with status OK
```

Exit code equals the number of folders that need action (0 = nothing to do), making
the script safe to use in shell conditionals or CI checks.
