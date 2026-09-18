# BFD Functional Annotation Pipeline (`--pipeline BFD`)

Nextflow DSL2 pipeline that runs functional-annotation and genome-statistics tools across
every species in `samples.csv`, merges each tool's per-genome output into a master Parquet
table set (`tables/`), and (optionally) rebuilds the queryable `db/BFD.duckdb` from those
tables.

Run after `funannotate.nf` has produced `genome_annotation/<Species_Strain>/predict_results/`
(see [Relationship to funannotate.nf](#relationship-to-funannotatenf)).

---

## Where things live (quick map)

```
funannotate.nf  ──▶  genome_annotation/<name>/predict_results/
                              │
                              ▼
   BFD.nf  (--pipeline BFD, this pipeline)
     SETUP_INPUT → RUN_* (function) + genome-stats/telomeres/BUSCO → MERGE_*
                              │
                              ▼
                    tables/*.parquet   (master, unscoped, full-rebuild-per-run)
                              │
                    (opt-in: --run_build_duckdb true, or run separately)
                              ▼
        nextflow/bin/build_BFD_duckDB.sh  ──▶  db/BFD.duckdb   (authoritative schema)
                              │
                    ┌─────────┴─────────┐
                    ▼                   ▼
         MCP/BFD_mcp_server.py   scripts/extract_bfd_taxonomic_subset.py
         (query layer)           (post-hoc taxon-scoped DB extract)
```

`sql/schema.sql` exists as a hand-maintained schema reference but has drifted from what
`build_BFD_duckDB.sh` actually creates (see [Known gaps](#known-gaps-worth-fixing)) — treat
`build_BFD_duckDB.sh` as the authoritative schema source, not `sql/schema.sql`.

---

## Running it

Launch via the single entry point `nextflow/main.nf`, always with `--pipeline BFD`:

```bash
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile BFD --pipeline BFD \
    -resume
```

Or via the SLURM wrapper (recommended for production):

```bash
sbatch nextflow/run_functional.sh
```

**Booleans (`--run_*`, `--merge_all`, `--skip_merge`, `--run_build_duckdb`) must go in a
`-params-file`, not as CLI `--flag` overrides** — the installed nextflow/nf-schema version's
`validateParameters()` rejects CLI-passed booleans for `BFD.nf` (confirmed 2026-08-26). Edit
`nextflow/param_files/params_functional.yaml` (paired with `run_functional.sh`) to change which
tools run. Non-boolean overrides (`--n_test`, `--taxon`, `--asmid`) are still safe on the CLI:

```bash
# Restrict to a taxonomic subset
sbatch nextflow/run_functional.sh --taxon PHYLUM:Ascomycota
sbatch nextflow/run_functional.sh --taxon CLASS:Dothideomycetes --n_test 5

# Test: first 5 species only
sbatch nextflow/run_functional.sh --n_test 5
```

To change tool selection, edit the params file itself, e.g.:

```yaml
# nextflow/param_files/params_functional.yaml
pipeline: BFD
run_setup: true
run_pfam: false
run_cazy: true
run_merops: true
run_swissprot: true
run_signalp: true
run_tmhmm: false
run_targetp: false
run_idp: true
run_wolfpsort: true
run_predgpi: false
run_busco_genome: true
run_busco_pep: true
```

There is also `nextflow/param_files/params_functional_all.yaml` (paired with
`run_functional_all.sh`) for the "everything on" variant.

---

## Pipeline stages

```
samples.csv
    │
    ▼
SETUP_INPUT              Create symlinks in input/ pointing to funannotate predict_results/
    │                    Skips species whose symlinks already exist.
    │                    (--run_setup false to skip entirely if input/ is pre-populated)
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Per-species functional annotation (all tools run in parallel per        │
 │  species, results/function/<tool>/, storeDir)                           │
 │                                                                          │
 │  RUN_PFAM        hmmsearch --cut_ga vs Pfam-A                           │
 │  RUN_CAZY        dbcanlight cazyme + substrate                          │
 │  RUN_MEROPS      blastp vs MEROPS scan.lib                              │
 │  RUN_SWISSPROT   blastp vs SwissProt (BUILD_SWISSPROT_ANNOT)            │
 │  RUN_SIGNALP     SignalP 6 (GPU, fast mode)                             │
 │  RUN_TMHMM       TMHMM short format                                     │
 │  RUN_TARGETP     TargetP 2 non-plant                                    │
 │  RUN_IDP         AIUPred disorder prediction (GPU)                      │
 │  RUN_WOLFPSORT   runWolfPsortSummary fungi                              │
 │  RUN_PREDGPI     predgpi.py GFF3 mode                                   │
 │                                                                          │
 │  All per-species processes use storeDir: automatically skipped if all    │
 │  output files already exist on disk (no -resume required).               │
 └──────────────────────────────────────────────────────────────────────────┘
    │
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Per-genome statistics (results/genome_stats/<type>/, storeDir,          │
 │  hash-bucketed by ASMID/LOCUSTAG — see T-014 in .living/decisions.md)    │
 │                                                                          │
 │  BATCH_AA_FREQ, BATCH_CODON_FREQ   AA/codon frequency (batched, default  │
 │                                    freq_batch_size=50 genomes/job)       │
 │  CALC_INTERGENIC, CALC_GENE_STATS  intergenic distances, gene structure  │
 │  CALC_ASM_STATS                    assembly QC stats (keyed by ASMID)   │
 │  FIND_TELOMERES     (--run_telomeres, default false)                    │
 │  BUSCO_GENOME       (--run_busco_genome, default false, ASMID-keyed)    │
 │  BUSCO_PEP          (--run_busco_pep,    default false, LOCUSTAG-keyed,│
 │                       post-annotation QC — NOT used for representative- │
 │                       strain selection, see decision F)                 │
 └──────────────────────────────────────────────────────────────────────────┘
    │
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  MERGE steps (one per tool, controlled by --merge_all / --skip_merge)   │
 │  → writes/overwrites the master, unscoped tables/*.parquet              │
 │                                                                          │
 │  merge_all=true : glob ALL files in results/function/<tool>/ and        │
 │    results/genome_stats/<type>/, regardless of which species were run   │
 │    this session. A sync barrier ensures newly-run files are on disk     │
 │    before the glob fires. MERGE is silently skipped if no files exist.  │
 │  merge_all=false (default): merge only the files produced in this run.  │
 │  skip_merge=true (default): skip all MERGE steps unconditionally.       │
 │                                                                          │
 │  There is no more per-taxon tables/<Taxon>/ subset build — MERGE_SAMPLES│
 │  always builds one full master samples/species table (T-014 §D.2). A    │
 │  taxon-scoped view is a post-hoc extraction, see below.                 │
 └──────────────────────────────────────────────────────────────────────────┘
    │
    ▼ (opt-in, --run_build_duckdb true)
BUILD_DUCKDB   wraps bin/build_BFD_duckDB.sh — full rebuild of db/BFD.duckdb
               from tables/*.parquet every time it runs (no incremental upsert)
```

**Defaults that actually matter for incremental runs** (`nextflow.config` /
`conf/profile_BFD.config`, current as of this rewrite):

| Param | Default | Effect |
|---|---|---|
| `merge_all` | `false` | MERGE sees only this run's output, not the whole `results/` tree |
| `skip_merge` | `true` | MERGE steps don't run at all unless you set this `false` |
| `run_build_duckdb` | `false` | `db/BFD.duckdb` is not rebuilt automatically — separate opt-in step |
| `run_busco_genome` / `run_busco_pep` | `false` / `false` | BUSCO steps are off by default in the BFD pipeline |
| `run_telomeres` | `false` | Telomere finder is off by default |
| all functional `run_*` (pfam/cazy/merops/swissprot/signalp/tmhmm/targetp/idp/wolfpsort/predgpi) | `true` | On by default in `-profile BFD` |
| `run_aa_freq` / `run_codon_freq` / `run_intergenic` / `run_gene_stats` | `true` | On by default |

If you want a run that both merges everything on disk into `tables/` **and** rebuilds
`db/BFD.duckdb` in one shot, set in your params file:

```yaml
merge_all: true
skip_merge: false
run_build_duckdb: true
```

---

## Output structure

```
input/                          ← populated by SETUP_INPUT (symlinks into genome_annotation/)
  pep/    <name>.proteins.fa
  cds/    <name>.cds-transcripts.fa
  gff3/   <name>.gff3
  dna/    <name>.scaffolds.fa
  trna/   <name>.trna.gff3

results/function/                ← per-species per-tool outputs (storeDir), *.gz
  pfam_hmmscan/  cazy/<name>/  merops/  swissprot/  signalp/  tmhmm/  targetP/
  aiupred/  wolfpsort/  predgpi/

results/genome_stats/             ← per-genome outputs (storeDir), hash-bucketed by
  aa_freq/  codon_freq/  intergenic/  gene_stats/  asm_stats/  BUSCO_genome/
  BUSCO_protein/  telomeres/           ASMID or LOCUSTAG per T-014 (not by species name;
                                       a generated read-only symlink tree,
                                       results/genome_stats_by_name/, keeps human-browsable
                                       access without making it the canonical store)

tables/                          ← master merged Parquet (MERGE_* output), one
  species.parquet                  unpartitioned file per type, always rebuilt full
  samples.parquet
  asm_stats.parquet
  telomere_summary.parquet / telomere_tracts.parquet
  busco_genome.parquet
  pfam.parquet  cazy.parquet  cazy_overview.parquet  merops.parquet  swissprot.parquet
  signalp.parquet  tmhmm.parquet  targetp.parquet  idp.parquet  idp_summary.parquet
  wolfpsort.parquet  predgpi.parquet
  aa_freq.parquet  codon_freq.parquet  gene_intergenic_distances.parquet
  gene_info.parquet  gene_transcripts.parquet  gene_exons.parquet  gene_CDS.parquet
  gene_introns.parquet  gene_trnas.parquet  gene_proteins.parquet

db/BFD.duckdb                    ← rebuilt from tables/*.parquet by bin/build_BFD_duckDB.sh
                                    (see below); NOT written directly by BFD.nf unless
                                    --run_build_duckdb true

logs/nextflow/
  BFD_trace.txt
  BFD_report.html
  BFD_timeline.html
```

---

## Building / rebuilding `db/BFD.duckdb`

`db/BFD.duckdb` is built by **`nextflow/bin/build_BFD_duckDB.sh`**, either:

- automatically at the end of a `BFD.nf` run, if `--run_build_duckdb true` is set (wired via
  the `BUILD_DUCKDB` module in `subworkflows/local/BFD_MERGE.nf`, gated on all `MERGE_*`
  outputs so it only fires after every merge in that run has published), or
- manually, any time after `tables/` is populated:

  ```bash
  module load duckdb
  bash nextflow/bin/build_BFD_duckDB.sh
  # or: sbatch nextflow/bin/build_BFD_duckDB.sh   (has its own #SBATCH header)
  ```

It is a **full rebuild every time** (`rm -f db/BFD.duckdb`, then `CREATE TABLE ... AS SELECT
... FROM read_parquet(...)` per table) — there is no incremental upsert, by design (see T-014
§D.1 in `.living/decisions.md`), so it's opt-in rather than tacked onto every merge run.

It builds these tables/views (30 tables, 3 views — the authoritative list, read directly out
of the script rather than assumed):

| Group | Tables |
|---|---|
| Species/assembly | `species`, `asm_stats` |
| Telomeres | `telomere_summary`, `telomere_tracts`, view `v_telomere_both_ends` |
| BUSCO | `busco_genome` |
| Gene structure | `gene_info`, `gene_proteins`, `gene_transcripts`, `gene_exons`, `gene_CDS`, `gene_introns`, `gene_trna`, `gene_intergenic_distances` |
| Functional annotation | `pfam`, `cazy`, `cazy_overview`, `merops`, `signalp`, `tmhmm`, `targetp`, `idp`, `idp_summary`, `wolfpsort`, `predgpi` |
| Sequence composition | `aa_frequency`, `codon_frequency` |
| Comparative genomics (from `--pipeline paralogoscope`, not BFD) | `wgd_ks`, `wgd_ks_summary` |
| Views | `v_species_summary`, `v_protein_annotation`, `v_telomere_both_ends` |

`species` is the authoritative join key table (`LOCUSTAG`, `ASMID` both unique-indexed); every
other table joins to it by `LOCUSTAG` or `ASMID` (assembly-level: `asm_stats`,
`telomere_summary`, `busco_genome`) or by a derived `species_prefix` column (annotation-level
tables). **The master DB is only as complete as what's been computed so far** — if 136 of
13,580 genomes haven't finished `FUNANNOTATE_TRAIN`/annotation yet, those genomes simply have
rows in `species`/`asm_stats` (assembly-level, independent of annotation) but no matching rows
in the annotation-derived tables (`gene_*`, `pfam`, `cazy`, etc.) until they finish. This was
directly verified during `EXTRACT_BFD_TAXONOMIC_SUBSET` testing (`GENUS:Malassezia`, 81 genomes
in `species`, only 1 with computed function/gene data at the time) — partial completion is
expected steady state, not an error condition.

### Taxon-scoped views

Since `MERGE_SAMPLES` only builds one master, unscoped table set (T-014 §D.2 — per-taxon
`tables/<Taxon>/` re-merges were retired as a side effect of manifest-glob scoping, not a
deliberate feature), a taxon-scoped database is a **post-hoc extraction**, not a pipeline-time
artifact:

```bash
python3 scripts/extract_bfd_taxonomic_subset.py \
    --taxon GENUS:Malassezia \
    --master db/BFD.duckdb \
    -o Malassezia.duckdb
```

It filters `species` by the taxon predicate first, then `INNER JOIN`s every other table against
that filtered `species` set (never re-evaluating the predicate independently per table), mirrors
`build_BFD_duckDB.sh`'s indexes, and re-executes the two analytical views read live from the
master (`duckdb_views()`) so view definitions never drift out of sync. `--dry-run` reports the
match count without writing.

### Querying

```bash
module load duckdb
duckdb db/BFD.duckdb
```

Or via the MCP server (`MCP/BFD_mcp_server.py`), which defaults to `db/BFD.duckdb` and reads
`PROTEIN_DB_PATH` to point at a different DB (e.g. a taxon-scoped extract):

```bash
PROTEIN_DB_PATH=/path/to/Malassezia.duckdb python MCP/BFD_mcp_server.py
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` | Master species/genome table (`ASMID, SPECIES_IN, STRAIN, BIOPROJECT, NCBI_TAXONID, BUSCO_LINEAGE, PHYLUM, SUBPHYLUM, CLASS, SUBCLASS, ORDER, FAMILY, GENUS, SPECIES, TRANSL_TABLE, LOCUSTAG`) |
| `--genome_annotation` | `genome_annotation/` | funannotate output root (source for SETUP_INPUT) |
| `--outdir` | `results/function/` | Per-species per-tool result files (storeDir root) |
| `--genome_stats_outdir` | `results/genome_stats/` | Per-genome stats storeDir root |
| `--tables` | `tables/` | Merged Parquet output |
| `--taxon` | `""` (all) | Restrict to `RANK:VALUE`, e.g. `PHYLUM:Ascomycota` |
| `--asmid` | `""` (all) | Restrict to a single ASMID |
| `--n_test` | `0` (all) | Limit to first N samples after taxon filter |
| `--suppress` | `suppress.txt` | ASMID suppress list applied across BFD/funannotate/ANI |
| `--run_setup` | `true` | Run SETUP_INPUT symlink step |
| `--run_pfam` / `--run_cazy` / `--run_merops` / `--run_swissprot` / `--run_signalp` / `--run_tmhmm` / `--run_targetp` / `--run_idp` / `--run_wolfpsort` / `--run_predgpi` | `true` (each) | Per-tool functional-annotation toggles (params-file only, see [Running it](#running-it)) |
| `--run_aa_freq` / `--run_codon_freq` / `--run_intergenic` / `--run_gene_stats` | `true` (each) | Genome-stats toggles |
| `--run_busco_genome` / `--run_busco_pep` | `false` / `false` | BUSCO steps (assembly / protein set) |
| `--run_telomeres` | `false` | Telomere finder |
| `--freq_batch_size` | `50` | Genomes per SLURM job for AA/codon frequency batching |
| `--merge_all` | `false` | MERGE from all files in results/ (not just current run) |
| `--skip_merge` | `true` | Skip all MERGE steps |
| `--run_build_duckdb` | `false` | Rebuild `db/BFD.duckdb` from `tables/*.parquet` after MERGE completes |
| `--pfam_nodes` | `1` | SLURM nodes per hmmscan job (`>1` enables MPI mode) |

### Strain ID sanitizing

When a row is read from `samples.csv` the `STRAIN` field is normalised before use:

1. Leading/trailing whitespace stripped.
2. Quotes (`'` `"`) removed.
3. Only the first `;`-delimited token is kept (some entries list synonyms).
4. **Colons (`:`) replaced with a space** — e.g. `CBS:123` → `CBS 123`.

The cleaned strain is then joined to the species name (`{SPECIES}_{STRAIN}`) and remaining
whitespace and path-unsafe characters are collapsed to `_` to form the per-species `basename`
used for `results/function/` output paths. Assembly- and annotation-level storeDir paths
(`results/genome_stats/`) instead key on `ASMID`/`LOCUSTAG` directly (T-014) — see
[Output structure](#output-structure).

---

## Skip/cache behavior

Per-species RUN processes use `storeDir`, so Nextflow skips a species automatically if all
declared output files already exist on disk — even on a fresh run without `-resume`.

### Automatic re-run when proteins are updated

Before dispatching each protein-consuming process (`RUN_*`, `BATCH_AA_FREQ`,
`BATCH_CODON_FREQ`), the pipeline compares the modification time of the input proteins/CDS
FASTA against the genome's cached output file. If the output is **strictly older** than the
input, it is deleted so the genome is included in the next batch job. This means re-running
`funannotate.nf` for a genome automatically triggers re-annotation by every BFD functional tool
the next time `BFD.nf` is submitted — no manual cache-busting required.

Processes whose inputs are not proteins (`CALC_INTERGENIC`, `CALC_GENE_STATS`) are not affected
by this check; delete their storeDir outputs manually if a re-run is needed. To force
re-running a species for a specific tool regardless of timestamps, delete its output file(s)
from `results/function/<tool>/` before re-submitting.

> **Note on `-stub-run`:** the staleness check and deletions happen at workflow initialisation
> time (inside channel operators), so they execute even under `-stub-run`. To dry-test without
> deletions, disable the relevant tool flags (e.g. `run_pfam: false` in the params file).

---

## Adding new species to an existing dataset

```bash
# 1. Annotate new genomes first (funannotate.nf)
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Mucoromycota

# 2. Run functional annotation on new species only.
#    Existing species are automatically skipped (storeDir).
sbatch nextflow/run_functional.sh --taxon PHYLUM:Mucoromycota

# 3. Merge everything on disk (existing + new) into tables/, and rebuild the DB.
#    Set in params_functional.yaml (booleans can't go on the CLI):
#      merge_all: true
#      skip_merge: false
#      run_build_duckdb: true
sbatch nextflow/run_functional.sh --taxon PHYLUM:Mucoromycota
```

To run new species and defer the merge until later (e.g. when multiple batches are running in
parallel), set `skip_merge: true` in the params file for the per-batch runs, then do one final
run with `skip_merge: false` (and `run_build_duckdb: true` if you also want `db/BFD.duckdb`
rebuilt) once all batches are done.

---

## Required modules

| Tool | Module |
|---|---|
| Pfam HMM scan | `hmmer/3.4`, `db-pfam` |
| CAZyme annotation | `dbcanlight` |
| Protease families | `db-merops/124`, `ncbi-blast/2.16.0+` |
| SwissProt | `ncbi-blast/2.16.0+`, `db-swissprot` |
| Signal peptides | `signalp/6-gpu` (A100 GPU) |
| TM helices | `tmhmm` |
| Subcellular targeting | `targetp` |
| Disorder prediction | `aiupred` (A100 GPU) |
| Subcellular localization | `wolfpsort` |
| GPI anchors | `predgpi` |
| Gene/sequence statistics | `biopython` |
| DB build/query | `duckdb` |

---

## Relationship to funannotate.nf

`funannotate.nf` produces `genome_annotation/<Species_Strain>/predict_results/`. `BFD.nf` reads
from that directory via its `SETUP_INPUT` step (symlinks into `input/`). Always run
`funannotate.nf` first, then `BFD.nf`.

`interproscan6.nf` uses a filter to skip genomes whose `annotate_misc/iprscan.xml` already
exists, with the same timestamp-based staleness check as above: if
`predict_results/<name>.proteins.fa` is newer than the existing `iprscan.xml`, the genome is
re-queued and InterProScan 6 re-runs for it.

---

## Known gaps worth fixing

Flagged during this rewrite, not fixed here:

- **`sql/schema.sql` has drifted from `build_BFD_duckDB.sh`.** It documents several tables the
  build script does not create (`funguild`, `chrom_info`, `mmseqs_orthogroup_clusters`,
  `mmseqs_orthogroup_cluster_count`, `prosite`, `pfam_UoT`, `gene_pairwise_distances`) and is
  missing several the build script does create (`telomere_summary`, `telomere_tracts`,
  `busco_genome`, `wolfpsort`, `predgpi`, `gene_CDS`, `gene_intergenic_distances`, `wgd_ks`,
  `wgd_ks_summary`, and all 3 views). Until reconciled, treat `build_BFD_duckDB.sh` as ground
  truth for the current schema.
- **`swissprot.parquet` is produced by `MERGE_SWISSPROT` but `build_BFD_duckDB.sh` does not
  load it into `db/BFD.duckdb`** — `run_swissprot` defaults `true` and the table is merged, but
  there's no `swissprot` table in the DB yet.
