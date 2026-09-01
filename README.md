# Fungi 22k project

Comparative genomics analysis of **~22k re-annotated fungal genomes** from NCBI (dataset frozen 2025-05).

Nextflow workflows to support automated and updating genome dataset. Synchronizing completed structure to [gs://stajichlab-fungi-bfd](gs://stajichlab-fungi-bfd)

1. Downloads from [@1kfg/NCBI_Fungi](https://github.com/1KFG/NCBI_fungi) project.
   * An automatic LOCUS_TAG prefix generated based on last 8 characters of MD5 hash of name of assembly from NCBI (for unique prefix and not-reusing NCBI LOCUS tags to avoid clashes) 
2. Genome cleanup and contaminant removal with [@stajichlab/AAFTF](https://github.com/stajichlab/AAFTF) and [NCBI FCS-GX](https://github.com/ncbi/fcs-gx)
3. Genome masking (for large scale using tantan for speed, if higher accuracy or TE comparison and compute allows will deploy RepeatMasker/RepeatModeler or EDTA or EarlGrey pipelines)
4. SRA query for RNA-Seq datasets for a SPECIES is are downloaded and run through a cleanup. Downloading up to 5 RNASeq sets. Limiting each to 50M reads per to avoid too many reads unnecessary for trinity assembly
   * Some of the issues - forcing headers to be formatted for trinity, dealing with BGI seq numbering issues
   * Forcing length of reads to be same for FWD and REV of a pair (requried for bbnorm)
   * custom scripts were written for this [fix_fastq_header_trinity])(https://github.com/hyphaltip/fix_fastq_header_trinity) and [enforce_seqpair_readlen](https://github.com/hyphaltip/enforce_seqpair_readlen). these will be folded into another sequence utility package later I expect
   * Running bbnorm.sh on the read pairs instead of using trinity normalization for speed (much faster- a summary of profiling run-time and diff in trinity results should be evaluated)
5. run Trinity through funannotate train - added a stop-after-trinity flag to funannotate to enable this. This way a single trinity transcript assembly is available to the 1-many strain assemblies for a species
6. funannotate train will run with this trinity transcript set and normalized reads as requested producing a PASA mysqldb locally for speed and this can be reused if UPDATE is run to get UTRs
7. funnannotate predict will run with the trained PASA models or just BUSCO models if no RNASeq. Will also skip genemark runs on genome assemblie
   * Fungal-only swissprot proteins were extracted as the input target proteins for informing gene models to avoid running full swissprot DB which saves time when running across 20k datasets. These protein models help inform slighly the gene models when running EVM but it isn't clear how much they really elp.

Command line options
===
* in funannotate.nf and other workflows --taxon "GENUS:Saccharomyces" to specify running on a taxonomic group. --asmid ABC123 will restrict to a specific assembly.
* --supress gives list of assembly accessions to not run on (based on previously determined criteria, eg too short, contaminated, etc) allows skipping without having to edit samples.csv file, can skip microsporidia this way for example


Goals
=====
 * BFD workflow runs protein annotation, summary stats on assemblies, codon usage, 
 * [todo] Deploy RHIEPA as devleoped in [@alanmoses](https://github.com/alanmoses) lab
 * [todo] phyling runs on 
---

## Repository layout

```
samples.csv              ← master species table (primary input for all pipelines)
input/
  pep/                   ← protein FASTAs  ({Species_Strain}.proteins.fa)
  cds/                   ← CDS transcript FASTAs  ({Species_Strain}.cds-transcripts.fa)
  gff3/                  ← GFF3 annotation files  ({Species_Strain}.gff3)
  dna/                   ← genome FASTA files  ({Species_Strain}.scaffolds.fa)
  trna/                  ← tRNA GFF3 files  ({Species_Strain}.trna.gff3)
genome_annotation/       ← funannotate predict output  ({Species_Strain}/predict_results/)
tables/                  ← consolidated CSV.gz files loaded into DuckDB (Nextflow output)
bigquery/                ← consolidated CSV.gz files (legacy SLURM script output)
results/                 ← per-tool output files
nextflow/                ← Nextflow pipelines
scripts/                 ← Python/R helper scripts
sql/schema.sql           ← DuckDB table definitions
functionalDB/            ← function.duckdb  (full annotation database)
intronDB/                ← intron_db.duckdb (gene-structure database)
Phylogeny/               ← phylogenetic pipeline and outputs
```

The `input/` subdirectories are populated by the `SETUP_INPUT` step of `BFD.nf`,
which creates per-species symlinks into `genome_annotation/{Species_Strain}/predict_results/`.
Run `funannotate.nf` first to produce those predict results, then `BFD.nf`.

---

## Primary input file: `samples.csv`

Every pipeline derives its sample list from this file. Key columns:

| Column | Description |
|---|---|
| `ASMID` | NCBI assembly accession |
| `NCBI_TAXONID` | NCBI taxonomy ID |
| `BUSCO_LINEAGE` | BUSCO lineage for quality assessment |
| `PHYLUM` … `GENUS` | Taxonomic classification |
| `SPECIES` | Binomial species name |
| `STRAIN` | Strain identifier (`;`-delimited for multiple) |
| `LOCUSTAG` | **8-character hex string** — stable genome identifier used in all gene IDs and DB foreign keys |

---

## Nextflow pipelines

All pipelines are launched from the single entry point `nextflow/main.nf`:

```bash
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile <profile> \
    --pipeline <pipeline> \
    -resume
```

| `--pipeline` | `-profile` | Launcher script | Description |
|---|---|---|---|
| `funannotate` | `funannotate` | `run_funannotate.sh` | Gene prediction + functional annotation |
| `BFD` | `BFD` | `run_functional.sh` | Functional annotation + genome statistics |
| `compare_ani` | `ani` | `run_ANI.sh` | All-vs-all ANI clustering |
| `query_ani` | `ani_query` | `run_ANI.sh` | ANI query against existing sketches |
| `earlgrey_mask` | `earlgrey` | `run_earlgrey.sh` | EarlGrey repeat masking |
| `comparative` | `comparative` | `run_comparative.sh` | Comparative genomics clustering |
| `paralogoscope` | `paralogoscope` | `run_paralogoscope.sh` | Whole-genome duplication dating (wgd dmd → ksd [+syn]) |
| `phyling` | `phyling` | `run_phyling.sh` | PHYling phylogenomics |

All workflows accept `--taxon RANK:VALUE` (e.g. `--taxon PHYLUM:Ascomycota`) and
`--n_test N` to restrict to the first N samples.

### Funannotate (`--pipeline funannotate`)

Genome cleaning → repeat masking → RNA-seq discovery/download → funannotate train → funannotate predict → (optional) annotate.

```bash
# Full run
sbatch nextflow/run_funannotate.sh

# Restrict to a clade
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Ascomycota

# Pilot: first 2 samples only
sbatch nextflow/run_funannotate.sh --n_test 2
```

See [README_funannotate.md](README_funannotate.md) for details.

### BFD functional annotation (`--pipeline BFD`)

Runs 9 functional annotation tools across all species, then merges results into `tables/`.

```bash
# Full run (run after funannotate completes)
sbatch nextflow/run_functional.sh

# Restrict to a clade
sbatch nextflow/run_functional.sh --taxon PHYLUM:Ascomycota
```

See [README_BFD.md](README_BFD.md) for details.

### ANI comparison (`--pipeline compare_ani`)

All-vs-all average nucleotide identity within taxonomic groups.

```bash
# All genera, default thresholds
nextflow run nextflow/main.nf -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani -resume

# Only Ascomycota, at family level
nextflow run nextflow/main.nf -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    --taxon PHYLUM:Ascomycota --compare FAMILY -resume

# Via the SLURM wrapper (recommended for production)
sbatch nextflow/run_ANI.sh
```

See [README_ANI.md](README_ANI.md) for details.

### Paralog detection & WGD dating (`--pipeline paralogoscope`)

Per-species detection and dating of whole-genome duplication events using
[wgd](https://github.com/arzwa/wgd) (container `container_wgd2_complete`,
`params.wgd_sif`). Runs the wgd pipeline inside an apptainer container on the
quadruple-decoded CDS set:

1. `wgd dmd` — all-vs-all DIAMOND + MCL to delineate paralog families
2. `wgd ksd` — Ka/Ks (Ks) estimation per family via MSA + codeml, giving the
   Ks distribution used to date WGD events
3. `wgd syn` — synteny / i-ADHoRe anchor detection for collinear blocks
   (optional, behind `--run_wgd_syn true`; heavy)

Inputs follow the BFD structure (`input/cds/{Species_Strain}.cds-transcripts.fa`
plus `input/gff3/` for syn); selection mirrors the other pipelines
(`--taxon`, `--group`, `--ignore`, `--n_test`).

```bash
# Full run, restricted to a clade (submitted to SLURM, preempt queue)
sbatch nextflow/run_paralogoscope.sh --taxon CLASS:Sordariomycetes

# Pilot: first 2 samples only
sbatch nextflow/run_paralogoscope.sh --n_test 2

# Enable the heavy i-ADHoRe synteny step
sbatch nextflow/run_paralogoscope.sh --run_wgd_syn true
```

Outputs follow the BFD hash sub-folder convention (no per-genome directory
accumulates too many files): results are fanned out by SHA-1 hash bucket keyed
on the stable LOCUSTAG under `paralogoscope/wgd_{dmd,ksd,syn}/{bucket}/`

```text
paralogoscope/wgd_dmd/{bucket}/{Species_Strain}.cds-transcripts.fa.tsv            ← paralog families
paralogoscope/wgd_ksd/{bucket}/{Species_Strain}.cds-transcripts.fa.tsv.ks.tsv    ← per-family Ks
paralogoscope/wgd_ksd/{bucket}/{Species_Strain}.cds-transcripts.fa.tsv.ksd.pdf   ← Ks distribution plot
paralogoscope/wgd_syn/{bucket}/anchors.csv                                       ← i-ADHoRe anchors
```

Setting these outputs as `storeDir` lets reruns reuse already-stored files
(it doubles as an output cache), and the natural wgd filenames mean a
re-annotated sample produces a new filename — a stale storeDir hit is never
possible.

After the last `wgd ksd` task the workflow's terminal merge (`MERGE_WGD_KSD`)
globs all published ks.tsv and streams them into
`tables/wgd.ks.parquet` (16-col zstd, `genome` = sampletag, `species_prefix` =
LOCUSTAG; includes the wgd tree-node id `node`, needed by wgd's own
mix/peak node-averaging) plus `tables/wgd.ks.summary.parquet` (n_pairs /
n_pairs_with_ds / n_families per genome) — the BFD MERGE_*/tablesDir()
pattern, so a fully-cached `-resume` still rebuilds the merged table. Those
parquet files also load into `db/BFD.duckdb` (`wgd_ks`, `wgd_ks_summary`)
when present. Per-genome Ks-peak summaries (number of duplicates + mean Ks
of each peak, via wgd's own GMM) are derived with
`analysis/WGD_PERFORMANCE_ANALYSIS/scripts/build_wgd_ksd_summary.py`.

---

## Recommended run order

```bash
# 1. Syntax check (no SLURM, no tools)
bash nextflow/run_lint.sh

# 2. Stub-run: validate DAG and output structure without real tools
bash nextflow/run_test.sh

# 3. Genome annotation (funannotate) → produces genome_annotation/
sbatch nextflow/run_funannotate.sh --n_test 2   # pilot
sbatch nextflow/run_funannotate.sh               # full run

# 4. Functional annotation + input setup (BFD)
sbatch nextflow/run_functional.sh --n_test 2    # pilot
sbatch nextflow/run_functional.sh               # full run

# 5. Per-species sequence statistics
sbatch nextflow/run_seqstats.sh --n_test 2      # pilot
sbatch nextflow/run_seqstats.sh                  # full run

# 6. ANI clustering (all-vs-all within taxonomic groups)
sbatch nextflow/run_ANI.sh                        # full run

# 7. Load results into DuckDB
sbatch pipeline/db/02_build_functional.sh
```

All Nextflow commands support `-resume` to restart from the last successful checkpoint.

---

## Loading results into DuckDB

After Nextflow pipelines complete, `tables/` contains one `.csv.gz` per table.

```bash
sbatch pipeline/db/02_build_functional.sh

module load duckdb
duckdb functionalDB/function.duckdb
```

See `sql/schema.sql` for all table definitions and indexes.

---

## Other analyses

**MMseqs2 protein clustering** (SLURM scripts, run from project root):

```bash
sbatch pipeline/01_cluster_mmseqs.sh           # cluster at 30% ID / 70% coverage
sbatch pipeline/02_process_cluster_pairwise.sh  # per-cluster pairwise distances
sbatch pipeline/03_make_mmseq_orthogroups.sh    # assign orthogroup IDs
```

**Phylogenetics** (run from `Phylogeny/`):

```bash
sbatch pipeline/01_phyling.sh
sbatch pipeline/02_phyling_filter_msa.sh
sbatch pipeline/03_phyling_make_tree.sh
sbatch pipeline/04_make_concatpartition.sh
```

---

## Notes

Large parts of this were written with Claude and some kimi model. These were mainly focused on porting previous bash-script arrayjob workflow to nextflow.  Logical flow and ability to keep up with nextflow language changes and some migration of logic were most easily done with these AI language help.

## Contact

Jason Stajich — jasonst@ucr.edu — Stajich Lab
