# Funannotate Pipeline (`funannotate.nf`)

Nextflow DSL2 pipeline for fungal genome prediction and annotation on SLURM HPC.
Covers the full funannotate workflow: genome cleaning → repeat masking → RNA-seq
download → train → predict → (optional) annotate.

Source synced from `../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf`.

---

## Quick start

```bash
# Default run: clean + mask + SRA fetch + train + predict (no annotate)
sbatch nextflow/run_funannotate.sh

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
    │                 filter short contigs → input_clean_genomes/<asmid>.fa
    │                 (storeDir; skipped if clean .fa already exists)
    ▼
MASKREPEAT_TANTAN_RUN Soft-mask repeats with tantan via funannotate mask
    │                 → input_clean_genomes/<asmid>.masked.fasta
    │                 (storeDir; skipped if masked file exists; skippable with
    │                  --run_repeatmasker false, falls back to clean .fa)
    ▼
SRA_FETCH             Search NCBI SRA for paired-end RNA-seq (≥250k read pairs,
    │                 75–300 bp, Illumina/BGI); download up to --max_rnaseq_runs sets.
    │                 Per-accession: parallel-fastq-dump → enforce read-pair length
    │                 → bbnorm (target=30) → fastp QC trim.
    │                 Writes rnaseq_reads/<species_tag>_norm_{R1,R2}.fastq.gz
    │                 (storeDir; empty files written when no SRA data found so the
    │                  cache is populated and downstream skips gracefully)
    │                 (skippable with --run_sra_fetch false)
    ▼
RNASEQ_PREPARE        Run funannotate train --stop_after_trinity on the
    │                 representative (first) assembly per species.
    │                 Archives Trinity-GG FASTA to rnaseq_data/<species_tag>.trinity-GG.fasta
    │                 (storeDir). All other strains of the same species reuse this.
    ▼
FUNANNOTATE_TRAIN     Run funannotate train on every assembly:
    │                 - With shared Trinity: only PASA runs (fast path)
    │                 - Without shared Trinity: full train (normalize + Trinity + PASA)
    │                 Removes large intermediates (hisat2/, trinity_gg/) after completion.
    │                 Skipped at channel level if funannotate_train.pasa.gff3 already exists.
    ▼
FUNANNOTATE_PREDICT   Run funannotate predict (Augustus + EvidenceModeler).
    │                 Uses training data linked from genome_annotation/<out>/training/.
    │                 Skipped if predict_results/<out>.gbk already exists.
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

## Output structure

```
genome_annotation/
  <Species_Strain>/
    predict_results/         ← primary output used by BFD.nf
      <name>.gbk             GBK (compressed .gz after predict)
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
  <species_tag>_norm_R1.fastq.gz    bbnorm + fastp normalized reads
  <species_tag>_norm_R2.fastq.gz
  rnaseq_manifest.tsv               provenance: species → SRA accessions

rnaseq_data/
  <species_tag>.trinity-GG.fasta    shared Trinity assembly per species

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
| `--run_sra_fetch` | `true` | Download RNA-seq from NCBI SRA |
| `--stop_after_sra_fetch` | `false` | Halt after SRA_FETCH (skip train/predict) |
| `--max_rnaseq_runs` | `4` | Max SRA accessions downloaded per species |
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
| `SRA_FETCH` | short → epyc | 24 | 48 → 192 GB | Retry bumps memory and queue |
| `RNASEQ_PREPARE` | epyc | 16 | 96 GB (+48 GB/retry) | Trinity-GG assembly; up to 3 retries |
| `FUNANNOTATE_TRAIN` | epyc | 8 → 24 | 96 → 192 GB | PASA alignment; up to 3 retries |
| `FUNANNOTATE_PREDICT` | epyc | 16 | 32 GB | |
| `SIGNALP_RUN` | exfab | 16 | 64 GB | Requires `--gres=gpu:1` |
| `INTERPROSCAN_RUN` | epyc | 8 | 32 GB | |
| `ANTISMASH_RUN` | epyc | 8 | 24 GB | |
| `FUNANNOTATE_ANNOTATE` | epyc | 16 | 64 GB | |

---

## Skip/cache behavior

All expensive one-time steps use `storeDir`, which means Nextflow skips the process if
all declared output files already exist on disk — even without `-resume`:

| Step | Skip condition |
|---|---|
| `SETUP_TAXONDB` | `taxondb/names.dmp` exists |
| `GENOME_CLEAN` | `input_clean_genomes/<asmid>.fa` exists |
| `MASKREPEAT_TANTAN_RUN` | `input_clean_genomes/<asmid>.masked.fasta` exists |
| `SRA_FETCH` | `rnaseq_reads/<tag>_norm_R1.fastq.gz` exists (even if 0-byte/no data) |
| `RNASEQ_PREPARE` | `rnaseq_data/<tag>.trinity-GG.fasta` exists |

`FUNANNOTATE_TRAIN` and `FUNANNOTATE_PREDICT` skip via channel-level file-existence
checks (`funannotate_train.pasa.gff3` and `predict_results/<out>.gbk`), so they also
skip gracefully when re-running over partially completed datasets.

---

## Required modules

```
miniconda3, AAFTF, taxonkit   — genome cleaning (FCS-GX)
funannotate                   — train, predict, annotate, mask
fastp, BBTools                — RNA-seq QC and normalization
sratoolkit, ncbi_edirect, parallel-fastq-dump  — SRA download
signalp/6-gpu                 — GPU signal peptide prediction (for annotate)
interproscan                  — domain annotation (for annotate)
antismash                     — BGC prediction (for annotate)
singularity                   — MariaDB container for PASA (if pasa_mysql=true)
```

---

## Relationship to BFD.nf

`funannotate.nf` produces `genome_annotation/<Species_Strain>/predict_results/`.
`BFD.nf` reads from that directory via its `SETUP_INPUT` step, which creates
symlinks in `input/pep/`, `input/cds/`, `input/gff3/`, `input/dna/`, and `input/trna/`.
Run `funannotate.nf` first, then `BFD.nf`.
