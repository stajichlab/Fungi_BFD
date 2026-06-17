# PLAN_ANI.md — compare_ANI Workflow Implementation Plan

## Overview

A Nextflow DSL2 workflow (`nextflow/compare_ANI.nf`) that groups fungal genomes from
`samples.csv` by a user-specified taxonomic rank, runs fastANI all-vs-all within each
group via a Singularity container, then post-processes results into a human-readable
clustering and outlier report.

---

## Files to Create

| File | Purpose |
|---|---|
| `nextflow/compare_ANI.nf` | Main workflow |
| `nextflow/conf/profile_ANI.config` | SLURM resource labels + Singularity |
| `nextflow/bin/report_ani.py` | Python: raw fastANI TSV → clustering report |
| `nextflow/params_ani.yaml` | Example params file (new-style NF params) |
| `nextflow/run_ANI.sh` | SLURM wrapper launch script |
| `PLAN_ANI.md` | This document |
| `README_ANI.md` | User-facing usage and documentation |

**Modified:** `nextflow/nextflow.config` — add `ani` profile entry.

---

## Parameters

All overridable via `params_ani.yaml` or `--param value` on the CLI.

| Parameter | Default | Description |
|---|---|---|
| `params.samples` | `samples.csv` | Master metadata CSV (inherited from shared config) |
| `params.genome_dir` | `input/dna` | Directory containing `*.scaffolds.fa` genome files |
| `params.compare` | `GENUS` | Column to group by: `GENUS`, `FAMILY`, `ORDER`, `CLASS`, `PHYLUM` |
| `params.taxon` | `''` | Pre-filter rows before grouping (e.g. `PHYLUM:Ascomycota`) |
| `params.ani_cluster_threshold` | `95.0` | ANI% cutoff for forming clusters in the report |
| `params.ani_outlier_threshold` | `90.0` | ANI% below which a genome is flagged as an outlier |
| `params.min_group_size` | `2` | Skip groups with fewer genomes than this |
| `params.ani_batch_size` | `0` | Max genomes per fastANI job; **0 = disabled** (default — use threads) |
| `params.outdir` | `results/ANI` | Root output directory |
| `params.fastani_fraglen` | `3000` | fastANI fragment length |
| `params.fastani_kmer` | `16` | fastANI k-mer size |
| `params.n_test` | `0` | Limit number of **groups** processed (0 = all); applied after groupTuple |

---

## Channel Architecture

```
samples.csv
  ├─ splitCsv(header: true)
  ├─ taxonFilter (--taxon RANK:VALUE)
  ├─ map { row → (group_key, locustag, basename, genome_path) }
  ├─ filter { genome file exists }        [warn on missing]
  ├─ groupTuple by group_key              [collect all genomes per group]
  └─ filter { group.size() >= min_group_size }

For each group:
  ├─ if group.size() <= ani_batch_size (or batching disabled):
  │     ANI_COMPARE(group_name, all_genomes, all_genomes, "full")
  │     → REPORT_ANI
  └─ if group.size() > ani_batch_size:
        split into batches of ani_batch_size
        → upper-triangle batch pairs (i, j where i <= j)
        → ANI_COMPARE(group_name, batch_i, batch_j, "b${i}_b${j}")
        → collect all TSVs for group_name
        → MERGE_ANI_BATCHES(group_name, tsv_list)
        → REPORT_ANI
```

---

## Processes

### `ANI_COMPARE`

- **Label:** `fastani`
- **Container:** `https://depot.galaxyproject.org/singularity/fastani:1.34--hb66fcc3_5`
  (resolved via `cacheDir = '/bigdata/stajichlab/shared/lib/singularity_cache'`)
- **storeDir:** `${params.outdir}/${params.compare}/${group_name}/batches`
  — if the output TSV already exists there, the process is skipped on `-resume`
- **Input:** `tuple val(group_name), path(query_genomes, stageAs:'query/*'), path(ref_genomes, stageAs:'ref/*'), val(batch_tag), val(group_size)`
  — `stageAs` stages query and ref into separate subdirs so identical files don't collide when the same genome list is passed for both (single-job / diagonal-batch case)
- **Output:** `${group_name}.${batch_tag}.ani.tsv`
- **Script:**
  ```bash
  ls query/* > query_list.txt
  ls ref/*   > ref_list.txt
  fastANI --ql query_list.txt --rl ref_list.txt \
      -o ${group_name}.${batch_tag}.ani.tsv \
      --fragLen ${params.fastani_fraglen} \
      -k ${params.fastani_kmer} \
      -t ${task.cpus}
  ```
- **Dynamic resources** (based on total `group_size`, not batch size):
  - ≤ 200 genomes → 8 CPUs / 16 GB
  - ≤ 500 genomes → 24 CPUs / 48 GB
  - \> 500 genomes → 64 CPUs / 128 GB
- **Fixed resources:** 24 h wall time, `epyc` queue; Singularity image cached at shared cache dir (all nodes)

### `MERGE_ANI_BATCHES`

- **Label:** `report`
- **storeDir:** `${params.outdir}/${params.compare}/${group_name}`
- **Input:** `tuple val(group_name), path(batch_tsvs)`
- **Output:** `${group_name}.ani.tsv`
- **Script:** `cat` all batch TSVs into one merged file
- **Resources:** 1 CPU, 4 GB RAM, 1 h, `short` queue

### `REPORT_ANI`

- **Label:** `report`
- **publishDir:** `${params.outdir}/${params.compare}/${group_name}`, mode `copy`
- **Input:** `tuple val(group_name), path(ani_tsv), path(names_tsv)`
- **Output:** `${group_name}_ANI_report.txt`
- **Script:** calls `bin/report_ani.py`
- **Resources:** 1 CPU, 4 GB RAM, 1 h, `short` queue

---

## `bin/report_ani.py` Logic

1. Parse fastANI TSV: `query  reference  ANI  mapped_frags  total_frags`
2. Strip path prefixes (basename only); exclude self-comparisons
3. Deduplicate bidirectional pairs (keep max ANI per unordered pair)
4. **Union-Find clustering** at `--cluster-threshold` (default 95.0%)
5. Identify **outliers**: genomes whose max pairwise ANI to any other genome
   in the group is below `--outlier-threshold` (default 90.0%)
6. Write report with sections:
   - Header: group name, N genomes, date, thresholds used
   - Clusters (sorted by size descending): N members, median ANI, min ANI, member list
   - Outliers: genome filename + species name, best ANI found

---

## Output Structure

Single-job mode (default, `ani_batch_size = 0`):
```
results/ANI/
└── GENUS/
    ├── Fusarium/
    │   ├── batches/
    │   │   └── Fusarium.full.ani.tsv        # raw fastANI output (storeDir)
    │   ├── Fusarium_genome_names.tsv        # filename → species name mapping (publishDir)
    │   └── Fusarium_ANI_report.txt          # human-readable report (publishDir)
    └── Aspergillus/
        └── ...
```

Batched mode (`ani_batch_size > 0`):
```
results/ANI/
└── GENUS/
    ├── Fusarium/
    │   ├── batches/                         # per-batch TSVs (storeDir)
    │   │   ├── Fusarium.b0_b0.ani.tsv
    │   │   └── ...
    │   ├── Fusarium.ani.tsv                 # merged raw output (storeDir)
    │   ├── Fusarium_genome_names.tsv        # filename → species name mapping (publishDir)
    │   └── Fusarium_ANI_report.txt          # human-readable report (publishDir)
    └── Aspergillus/
        └── ...
```

---

## Batching Strategy

Triggered when `ani_batch_size > 0` and `group_size > ani_batch_size`.

For N genomes and batch size B:
- Number of batches: `K = ceil(N / B)`
- Number of jobs: `K*(K+1)/2` (upper triangle including diagonal)
- All jobs run in parallel on SLURM
- Maximum comparisons per job: `B × B`

Example with `ani_batch_size = 500`, Fusarium N=1941:
- 4 batches → 10 batch-pair jobs → max 250k comparisons per job
- Each job: ~5–10 min with 8 threads

Each batch job's output is cached independently in `batches/` via `storeDir`.
Interrupted runs resume from the last completed batch.

---

## Caching / Resume

- `ANI_COMPARE` uses `storeDir` so completed groups/batches are not recomputed
  even across separate `-resume` sessions
- If a genome FASTA changes, manually delete the group's storeDir folder
  (e.g. `results/ANI/GENUS/Fusarium/`) and re-run
- `REPORT_ANI` uses `publishDir` + standard Nextflow fingerprint caching

---

## Singularity

- Global `nextflow.config` already sets:
  - `singularity.autoMounts = true`
  - `singularity.cacheDir   = '/bigdata/stajichlab/shared/lib/singularity_cache'`
- The ANI profile enables singularity and sets the container on the `fastani` label
- The SIF image is pulled once to the shared cache and accessible from all nodes

---

## Two ANI Thresholds

| Parameter | Purpose |
|---|---|
| `ani_cluster_threshold` (95%) | Minimum ANI to place two genomes in the same cluster |
| `ani_outlier_threshold` (90%) | A genome is an outlier if its best hit in the group is below this |

Both appear in the report header so results are unambiguous.
