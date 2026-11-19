# Building the BFD database

How `tables/*.parquet` and `db/BFD.duckdb` are produced, what launches each
step, and the one ordering rule that must never be broken.

Audience: anyone rebuilding BFD or consuming it. If you only want to *query*
the database, you do not need this file — see `docs/BFD_Introduction.md`.

## 1. What is actually built

Two artefacts, in this order:

| artefact | where | what it is |
|---|---|---|
| `tables/*.parquet` | `Fungi_BFD_runs/tables/` | one unpartitioned Parquet file per table, the real data |
| `db/BFD.duckdb` | `Fungi_BFD_runs/BFD_function/db/` | **views** over those Parquet files, by absolute path |

The database is almost entirely views. Only five objects hold data:

- `swissprot_subject` (575,503 rows), materialized by the builder
- `ref_fungaltraits`, `ref_growth_form_override`, `ref_class_subphylum`,
  `ref_phylum_subphylum`, copied in from `lib/refs/` by `sql/bfd_traits.sql`

Two consequences worth internalising:

1. **A view reads its Parquet at query time.** Republish a Parquet file and
   every existing database serves the new rows immediately, with no rebuild.
   This is why the 2026-11-18 taxonomy republish reached both databases
   without either being rebuilt.
2. **The absolute path is baked into each view.** A database only works where
   `tables/` sits at the path it was built against. Moving `tables/` breaks
   every database built before the move.

## 2. The one rule

> **Never build the database while a `MERGE_*` step is writing `tables/`.**

`publishDir` copies in the background after a task ends, so a Parquet file can
exist at its final path while still being written. On 2026-11-06 `BUILD_DUCKDB`
read a half-copied `gene_proteins.parquet` and failed with "No magic bytes
found". `BUILD_DUCKDB` now waits for each published file to match the staged
size and end in `PAR1`, but that guard only covers files from *its own* run —
it cannot see a merge launched separately.

`rebuild_pipeline_duckdb.sh` enforces this by refusing to start when a Nextflow
run is live in the same directory. Nothing enforces it across directories;
that is on you.

## 3. Stage one — merging into `tables/`

Every table is produced by a `MERGE_*` Nextflow process that globs
`results/function/<tool>/` and writes one Parquet file.

### Control parameters

Set in `nextflow/nextflow.config`, overridden per run by a params file.

| parameter | meaning |
|---|---|
| `skip_merge` | `true` (default) skips every `MERGE_*` step |
| `merge_all` | `true` merges all files under `results/`, `false` only this run's |
| `merge_if_stale` | `true` runs a merge only when a result file is newer than its table, or the table is missing |
| `skip_merges` | comma-separated labels never to run |
| `force_merges` | comma-separated labels to run even when `merge_if_stale` says current — **use after changing a merge's code**, since its inputs are not newer |
| `coverage_min_frac` | a full merge refuses to publish a table covering less than this fraction (0.9) of the genomes in the table it replaces |
| `allow_shrink` | turns that check, and `MERGE_SAMPLES`' no-shrink check, off for a deliberate shrink |

### How a table gets published

Most merges publish by `publishDir` straight into `params.tables`.
`species.parquet` is different: `MERGE_SAMPLES` calls
`nextflow/bin/publish_table.py`, which writes the Parquet by atomic rename
under a per-table `flock` and then writes
`tables/_manifest/<table>.json` recording `row_count`, `built_at`, `built_by`,
`merge_run_id` and `parquet_sha256`.

**Consumers should gate on the manifest, not on file mtime.** `cp -p` preserves
mtimes, so a freshly installed file can look old.

### Launching a merge-only run

From `Fungi_BFD_runs/BFD_function/`:

```bash
sbatch run_BFD.sh          # with PARAMS_FILE pointing at a params file
# or directly:
nextflow run <checkout>/nextflow/main.nf -profile BFD \
    -params-file nextflow/param_files/params_merge_only.yaml -resume
```

`params_merge_only.yaml` turns every producer off and sets
`merge_all: true`, `merge_if_stale: true`,
`skip_merges: "MERGE_TARGETP,MERGE_INFERNAL"`, `run_build_duckdb: true`.

TargetP and Infernal are skipped because their results cover only 125 and 15
genomes; merging them would trip the coverage check.

**Known limitation:** deleting a genome's result files does not make the merge
rerun. Delete the affected table in `tables/` to force it.

## 4. Stage two — building `db/BFD.duckdb`

One script does the work in all cases: `nextflow/bin/build_BFD_duckDB.sh`.
It is a **full rebuild every time** (`rm -f "$DB"`); there is no incremental
upsert.

It builds in three passes:

1. **Catalog views** — `species`, `asm_stats`, `busco_genome` and the rest,
   generated from `sql/table_schema.json` by `bin/generate_bfd_catalog_sql.py`.
2. **`v_protein_annotation`** — the per-protein join. Note the
   `WHERE feature = 'GPI-anchor'` filter on `predgpi`: without it every protein
   PredGPI read is flagged `is_gpi_anchored` (223.9 M instead of 3.3 M).
3. **Compat, traits and helper SQL** — `sql/bfd_compat_views.sql`,
   `sql/bfd_traits.sql`, `sql/bfd_extra_views.sql`, applied in that order
   (`bfd_traits.sql` defines `v_genome_traits`, which `bfd_extra_views.sql`
   reads).

Pass 3 moved here from `BFD_dashboard` on 2026-11-19. Before that only the
dashboard wrapper applied it, so the pipeline database had 39 views and the
dashboard's had 59, and they drifted between builds.

### The line-1 path contract

Each of the three SQL files begins with its own default data path:

```sql
-- TABLES_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
-- REF_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs
```

That line is **not documentation — it is a substitution key.** The builder runs
`sed "s#<line-1 path>#<real path>#g"` over the file. A path in the body that
line 1 does not spell exactly survives the rewrite, and the database silently
reads whatever sits at the stale path.

`nextflow/tests/test_sql_path_headers.py` guards this. Run it after editing any
of the three files.

(`bfd_extra_views.sql` has no absolute read at all — every view is built on
other views — so its header is inert. The builder still requires it.)

### Environment contract

| variable | default | meaning |
|---|---|---|
| `TABLES` | `tables` | directory of Parquet files |
| `DBDIR` | `db` | where `BFD.duckdb` is written |
| `BIN_DIR` | *required* | absolute path to `nextflow/bin` |
| `SCHEMA` | *required* | absolute path to `sql/table_schema.json` |
| `SQL_DIR` | `$BIN_DIR/../../sql` | the three SQL files |
| `REF_DIR` | `$BIN_DIR/../../lib/refs` | FungalTraits reference data |

## 5. The three ways it gets launched

### A. Inside a pipeline run — the normal path

Set `run_build_duckdb: true`. The `BUILD_DUCKDB` process runs after every
`MERGE_*` step that ran this session, gated on their outputs, and publishes to
`params.db` (default `${launchDir}/db`).

Off by default, because it is a full rebuild every time.

Its Nextflow cache key covers `bin/build_BFD_duckDB.sh`,
`sql/table_schema.json`, the three SQL files and `lib/refs/`. **Change any of
them and the next `-resume` rebuilds.** Add anything new the builder reads to
that list in `modules/BFD/BUILD_DUCKDB/main.nf`, or `-resume` will serve a
stale database.

### B. Builder changed, data did not

From `Fungi_BFD_runs/BFD_function/`:

```bash
sbatch rebuild_pipeline_duckdb.sh
```

Runs the builder alone, no Nextflow, no merges. Builds into a staging directory
and moves the result into place at the end, so a failed build keeps the old
database. Refuses to run when a Nextflow run is live in the same directory.

This is the right tool after a change to the builder or the SQL.

### C. By hand

```bash
TABLES=/abs/tables DBDIR=/abs/db \
BIN_DIR=<checkout>/nextflow/bin \
SCHEMA=<checkout>/sql/table_schema.json \
    bash <checkout>/nextflow/bin/build_BFD_duckDB.sh
```

For a scratch copy. Do not point `DBDIR` at a published database directly —
use B, which stages.

## 6. Where things live

**`Fungi_BFD`** (this repo) — the pipeline, and everything version-controlled
that defines *how* a build happens: `nextflow/`, `sql/`, `lib/refs/`,
`nextflow/tests/`.

**`Fungi_BFD_runs`** — *this* instance of the build. Everything that says
*which* genomes and *which* run: `samples.csv`, `suppress.txt`, `tables/`,
`BFD_function/db/`, and the per-run `params_*.yaml` files. `BFD_function/`
links `samples.csv` and `suppress.txt` back up to the run root, so a release is
defined by what sits there.

To cut a different release, vary `samples.csv` / `suppress.txt` /
`params_*.yaml` in a run directory — not the pipeline.

**`BFD_dashboard`** — consumer only: Quarto figures and site pages, slimmed
WASM database subsets, atlas and summary tooling. It should not define schema.

## 7. Verifying a build

```sql
SELECT count(*) FROM duckdb_views() WHERE NOT internal;   -- expect 59
SELECT table_name FROM duckdb_tables() ORDER BY table_name;
-- expect: ref_class_subphylum, ref_fungaltraits,
--         ref_growth_form_override, ref_phylum_subphylum, swissprot_subject
```

Then check a view against its Parquet, and the manifest:

```sql
SELECT (SELECT count(*) FROM species) AS view_rows;
```

```bash
cat tables/_manifest/species.json      # row_count must match
```

`nextflow/bin/check_traits_vs_atlas.py` compares `v_genome_labels` traits
against an atlas Parquet and the FungalTraits CSV; it exits 1 on a mismatch.

Run the test suite after any change to the builder or SQL:

```bash
python3 -m pytest nextflow/tests/test_sql_path_headers.py \
                  nextflow/tests/test_traits_sql.py \
                  nextflow/tests/test_traits_refs.py -q
```

## 8. Troubleshooting

**"No magic bytes found"** — the builder read a Parquet file mid-copy. A merge
was running. See §2.

**A view's row count does not match its Parquet** — the database was built
against a different `tables/` directory. Check the path baked into the view:
`SELECT sql FROM duckdb_views() WHERE view_name = '<name>'`.

**`is_gpi_anchored` is true for almost every protein** — the database predates
the `feature = 'GPI-anchor'` filter. Rebuild with B.

**Traits columns are null, or `subphylum_source` is unexpected** — check
`lib/refs/` is present and `sql/bfd_traits.sql`'s line-1 `REF_DIR` matches the
paths in its body (§4).

**`-resume` served a stale database** — something the builder reads is missing
from `BUILD_DUCKDB`'s cache key (§5A).
