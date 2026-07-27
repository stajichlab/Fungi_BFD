# compare_ANI — Average Nucleotide Identity Workflow

Compares all genomes within a taxonomic group (e.g. all Fusarium genomes) all-vs-all,
then reports pairwise clusters and outliers.

## Identity methods (`--ani_method`)

fastANI all-vs-all is O(N²) with per-pair fragment mapping, which is too slow at
clade level (CLASS / PHYLUM = thousands of genomes). The workflow now supports four
methods; **skani is the default**. All sketch-based methods sketch each genome
**once** into a shared cache, so the expensive step is linear in N.

| `--ani_method` | Tool | Speed | Notes |
|---|---|---|---|
| `skani` *(default)* | [skani](https://github.com/bluenote-1577/skani) | fast | Accurate ANI down to ~80%; sketch-once + sparse `triangle`. Best for clades. |
| `mash` | [mash](https://github.com/marbl/Mash) | fastest | MinHash distance → ANI = 100·(1−dist); rough below ~90%. |
| `sourmash` | [sourmash](https://sourmash.readthedocs.io) | fast | Containment-ANI matrix (`compare --ani`); dense N×N (memory-heavy on big clades). |
| `fastani` | [fastANI](https://github.com/ParBLiSS/FastANI) | slow | Classic mapping ANI. Add `--fastani_prefilter true` to mash-cluster first and run fastANI only *within* components. |

Every method emits a normalized `query⇥reference⇥ANI` table and runs through the
same report step, so clustering/outlier output is identical in format regardless of
method.

### Sketch cache

skani/mash/sourmash write one sketch per genome under `--sketch_cache`
(default `work/ANI/sketch_cache/<method>/<params>/`) using Nextflow `storeDir`.
A genome is sketched once and the sketch is reused across groups, across `-resume`
runs, and across different `--compare` levels (e.g. a GENUS run then a CLASS run).
Delete the cache to force re-sketching.

---

## Quick Start

All commands use the single entry point `nextflow/main.nf` with `--pipeline compare_ani`:

```bash
# All genera, default thresholds
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    -resume

# Only Ascomycota, at family level
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    --taxon PHYLUM:Ascomycota --compare FAMILY \
    -resume

# Clade level with the fastest method
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    --compare CLASS --ani_method mash \
    -resume

# fastANI, but mash-prefiltered so it only compares within components
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    --ani_method fastani --fastani_prefilter true \
    -resume

# Using a params file
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani --pipeline compare_ani \
    -params-file nextflow/params_ani.yaml \
    -resume

# Via the SLURM wrapper (recommended for production)
sbatch nextflow/run_ANI.sh
```

---

## Parameters

All parameters can be set on the command line (`--param value`) or in a YAML params
file (`-params-file nextflow/params_ani.yaml`).

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` | Master metadata CSV |
| `--genome_dir` | `input/dna` | Directory of `*.scaffolds.fa` genome files |
| `--compare` | `GENUS` | Taxonomic rank to group by. One of: `PHYLUM`, `SUBPHYLUM`, `CLASS`, `SUBCLASS`, `ORDER`, `FAMILY`, `GENUS` |
| `--taxon` | _(all)_ | Pre-filter rows before grouping. Format: `RANK:VALUE`, e.g. `PHYLUM:Ascomycota` |
| `--ani_method` | `skani` | Identity method: `skani`, `mash`, `sourmash`, `fastani` |
| `--sketch_cache` | `work/ANI/sketch_cache` | Per-genome sketch cache (skani/mash/sourmash) |
| `--ani_cluster_threshold` | `95.0` | ANI% cutoff for cluster formation in the report |
| `--ani_outlier_threshold` | `90.0` | ANI% below which a genome is flagged as an outlier |
| `--min_group_size` | `2` | Skip groups with fewer than this many genomes |
| `--outdir` | `results/ANI` | Root output directory |
| `--n_test` | `0` | Limit number of **groups** processed for testing (0 = all) |
| **skani** | | |
| `--skani_preset` | `medium` | `fast` / `medium` / `slow` (must match sketch + triangle) |
| `--skani_min_af` | `15` | Minimum aligned fraction (%) to report a pair |
| `--skani_compression` | `0` | `-c` override (0 = preset default) |
| `--skani_sketch_chunk` | `50` | Genomes sketched per skani job (batched to reduce job count) |
| **mash** | | |
| `--mash_kmer` | `21` | k-mer size |
| `--mash_sketch_size` | `10000` | sketch size (hashes per genome) |
| **sourmash** | | |
| `--sourmash_kmer` | `21` | k-mer size |
| `--sourmash_scaled` | `1000` | scaled factor (lower = more hashes = finer) |
| **fastANI** | | |
| `--fastani_fraglen` | `3000` | fastANI fragment length |
| `--fastani_kmer` | `16` | fastANI k-mer size (16 species-level, 13 divergent) |
| `--ani_batch_size` | `0` | Batch large groups into chunks (0 = single threaded job per group) |
| `--fastani_prefilter` | `false` | mash pre-cluster, then fastANI within components only |
| `--prefilter_ani` | `80.0` | mash-ANI% floor to join genomes in a prefilter component |

---

## params_ani.yaml Example

```yaml
samples:   "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
genome_dir: "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/input/dna"
compare:   "GENUS"
taxon:     "PHYLUM:Ascomycota"
ani_cluster_threshold:  95.0
ani_outlier_threshold:  90.0
min_group_size: 2
outdir:    "results/ANI"
```

---

## Output Structure

Output is segregated by method so switching `--ani_method` never reuses another
method's results:
```
results/ANI/
└── skani/                                   # value of --ani_method
    └── GENUS/                               # value of --compare
        ├── Fusarium/
        │   ├── batches/
        │   │   └── Fusarium.full.ani.tsv        # normalized q/r/ANI table
        │   ├── Fusarium_genome_names.tsv        # filename → species name mapping
        │   └── Fusarium_ANI_report.txt          # clustering + outlier report
        └── ...
```

For fastANI, when batching (`--ani_batch_size > 0`) or the prefilter cascade is
enabled, the `batches/` subfolder holds the per-batch / per-component TSVs and a
merged `<GROUP>.ani.tsv` is written at the group root.

---

## Report Format

```
=== ANI Report: Fusarium ===
Compare level  : GENUS
Genomes        : 1941
Cluster threshold  : 95.0%
Outlier threshold  : 90.0%
Generated      : 2026-06-05

--- Clusters (threshold >= 95.0%) ---

Cluster 1  N=1823  median_ANI=98.42%  min_ANI=95.12%
  Fusarium_oxysporum_f._sp._lycopersici_4287.scaffolds.fa   Fusarium oxysporum f. sp. lycopersici 4287
  Fusarium_oxysporum_race_1.scaffolds.fa                    Fusarium oxysporum race 1
  ...

Cluster 2  N=52  median_ANI=96.10%  min_ANI=95.02%
  ...

--- Outliers (best ANI < 90.0%) ---

  Fusarium_sp._NOVEL1.scaffolds.fa   Fusarium sp. NOVEL1   best_ANI=84.3%  best_match=Fusarium_sp._XYZ.scaffolds.fa
  ...
```

---

## Caching and Resume

The compare step (`SKANI_COMPARE` / `MASH_COMPARE` / `SOURMASH_COMPARE` /
`FASTANI_COMPARE`) uses Nextflow `storeDir`: if the output TSV already exists in the
results folder it is not recomputed, even across separate pipeline runs. Per-genome
sketches are cached independently under `--sketch_cache` (see above).

```bash
# Re-run safely; completed groups are skipped
nextflow run nextflow/compare_ANI.nf ... -resume

# Force recomputation of one group (e.g. after replacing genome files)
rm -rf results/ANI/GENUS/Fusarium/
nextflow run nextflow/compare_ANI.nf ... -resume
```

---

## Batching Large Groups

By default (`--ani_batch_size 0`) each taxonomic group runs as a single fastANI job
using all available threads. FastANI's sketch-based algorithm scales well; even
Fusarium (≈1941 genomes, ≈3.8M pairs) typically completes in a few hours with 8–16
threads.

For extremely large groups or strict wall-time limits, enable batching:

```bash
nextflow run nextflow/compare_ANI.nf ... --ani_batch_size 500
```

This splits the genome list into chunks of 500 and submits upper-triangle
batch-pair jobs (10 jobs for 1941 genomes) that run in parallel. Each job
handles at most 500×500 = 250k comparisons. Batch outputs are merged
automatically before reporting.

---

## Singularity

Every identity tool runs inside a BioContainers Singularity image — this is the
preferred path and no module load is needed for the compute tools. Images are
pulled once and cached in the shared Singularity cache
(`/bigdata/stajichlab/shared/lib/singularity_cache`), accessible from all nodes.

| Tool | Image (override with `--<tool>_container`) |
|---|---|
| skani | `skani:0.2.2--ha6fb395_2` |
| mash | `mash:2.3--hb105d93_10` |
| sourmash | `sourmash:4.9.4--hdfd78af_0` |
| fastANI | `fastani:1.34--hb66fcc3_5` |
| report / components | `python:3.12` |

All paths are prefixed with `https://depot.galaxyproject.org/singularity/`. The
report / component-finding step is pure-Python (stdlib only), so it runs in a plain
`python` image — **no `module load` anywhere in the pipeline**.

---

## Notes on ANI Thresholds

- **Cluster threshold (default 95%)**: A well-established species boundary for
  many bacteria; for fungi, intraspecific ANI is often >97%. Raise this to 98%
  for species-level clustering or lower it to 90% to identify genus-level groups.
- **Outlier threshold (default 90%)**: Genomes below this are flagged as potential
  misidentifications, novel lineages, or contaminants. They are listed separately
  for manual inspection, not automatically excluded.
- Both thresholds are reported in the output file header so results are
  self-documenting.
