# BFD data-lake redesign: views-only `db/BFD.duckdb` over `tables/*.parquet`

Date: 2026-09-09
Status: approved design, not yet implemented
Reviewed by: Fable (HPC storage/DB-performance persona), a fork agent (nextflow/bin audit)

## Problem

`db/BFD.duckdb` is built by `nextflow/bin/build_BFD_duckDB.sh` as a full `rm -f` +
`CREATE TABLE ... AS SELECT ... FROM read_parquet(...)` copy of every table (30 tables, 3
views), every time it runs. This is expensive enough that it's opt-in (`--run_build_duckdb`,
default `false`) rather than run after every merge.

Meanwhile different tables in `tables/` are populated on independent schedules — assembly-level
tables (`asm_stats`, `busco_genome`) are already largely complete; annotation-derived tables
(`gene_*`, `pfam`, `cazy`, ...) lag behind whatever fraction of ~13,580 genomes have finished
`funannotate` prediction (136 genomes still pending as of 2026-09-09). A downstream consumer,
DeltaGain, needs current BUSCO completeness + assembly N50 without waiting for annotation to
finish everywhere, and already reads `tables/busco_genome.parquet` / `tables/asm_stats.parquet`
directly via `read_parquet()`, bypassing `BFD.duckdb` entirely (see
`DeltaGain/docs/superpowers/notes/busco-n50-datalake-requirements.md`).

## Decision

`db/BFD.duckdb` stops holding materialized copies. It becomes a thin, cheap-to-rebuild **catalog**
of `VIEW`s over `tables/*.parquet`, generated from a load-bearing schema contract, backed by a
per-table freshness manifest and atomic writes. Four parts:

### 1. Catalog: `VIEW`s, not copies

For each logical table, `nextflow/bin/build_BFD_duckDB.sh` (kept as the existing name/path —
every consumer already points at `db/BFD.duckdb`, no churn) emits:

```sql
CREATE VIEW asm_stats AS
SELECT sp.LOCUSTAG, a.ASMID, a.contig_count, a.total_length_bp AS TOTAL_LENGTH, ...
FROM read_parquet('<abs>/tables/asm_stats.parquet') AS a
LEFT JOIN read_parquet('<abs>/tables/species.parquet') AS sp USING (ASMID);
```

Two changes from the old materialized build:

- **`LEFT JOIN`, not `INNER JOIN`**, fact table (`asm_stats`/`busco_genome`) → `species`
  (dimension). The old `INNER JOIN` silently dropped assembly rows with no `species` match
  (Fable found 2 such rows in the live data) — `LEFT JOIN` keeps them (with `LOCUSTAG` null),
  and a per-table row-count check (view row count vs. that table's own manifest `row_count`,
  see part 2) surfaces the mismatch instead of hiding it.
- **Explicit typed columns from `table_schema.yaml`** (part 3), not `SELECT *`. DuckDB 1.1.3
  invalidates a view outright on any schema change to the underlying Parquet, even adding a
  column (`Binder Error: Contents of view were altered: types don't match!`, confirmed by
  Fable) — the view stays broken until the catalog is rebuilt. Driving explicit types from the
  schema contract, and rebuilding the catalog in the same step as any merge that changes a
  table, avoids ever shipping a broken view.
- **`species_prefix`/`LOCUSTAG` on annotation tables (`gene_*`, `pfam`, `cazy`, ...) become real
  materialized columns written at merge time**, not a computed `string_split()` expression in
  the view. A computed expression can't use Parquet row-group min/max statistics, so a filter on
  it forces a full scan — a real cost once `gene_intergenic_distances` is back to its full
  ~38M-row size (Fable, citing the old CSV-era DB). `MERGE_*` scripts additionally `ORDER BY`
  this column when writing, so row groups are sorted and pruning actually works.

Rebuilding the catalog is a fixed cost (~30 `CREATE VIEW` statements) independent of data
volume, so it is cheap enough to run after every `MERGE_*` completes. **Recommendation**: flip
`--run_build_duckdb`'s default from `false` to `true` now that it's cheap, so the catalog is
always current after any merge run — flagged here explicitly since it's a behavior change from
today, not something to change silently.

Existing consumers (MCP server, `pick_representative_strain.py` — already Parquet-native per
the fork audit, `sql`-level ad-hoc queries) need no path or table-name change: a `SELECT`
against a `VIEW` is transparent versus a `TABLE`. Per-taxon extraction
(`scripts/extract_bfd_taxonomic_subset.py`) is unaffected — it still produces a real
materialized DB on explicit request; this redesign only changes the always-current master
catalog.

Indexes are gone (a `VIEW` can't carry a `UNIQUE INDEX`); point lookups rely on DuckDB's
Parquet row-group statistics/predicate pushdown instead, which the sorted-write change above
restores for the columns that matter.

### 2. Freshness manifest: one sidecar per table

`tables/_manifest/<name>.json`, one file per table (not one shared file — Fable flagged that a
single `tables/_manifest.json` written by ~30 concurrent `MERGE_*` processes loses updates under
last-writer-wins). Written by a shared helper every `MERGE_*` script calls immediately after its
own Parquet write, using the same tmp+rename pattern as part 4:

```json
{
  "table": "busco_genome",
  "path": "tables/busco_genome.parquet",
  "row_count": 23125,
  "built_at": "2026-09-09T14:22:03Z",
  "built_by": "MERGE_BUSCO_GENOME",
  "merge_run_id": "<nextflow session/task hash>",
  "schema_version": 1,
  "columns": ["ASMID", "complete_pct", "single_pct", "duplicated_pct", "fragmented_pct",
              "missing_pct", "n_markers", "lineage"],
  "parquet_sha256": "<sha256>"
}
```

`merge_run_id` lets a consumer confirm that a set of related tables (e.g. all 7 `gene_*` tables)
came from the same generation rather than mixing an old and a refreshed table. The catalog
exposes the same sidecars as a `_table_manifest` view via `read_json_auto('tables/_manifest/*.json')`
for SQL-side consumers; a consumer that only reads Parquet directly (DeltaGain's current
pattern) reads the JSON sidecar with no DuckDB dependency at all.

`schema_version` is a bare per-table int, versioned in `table_schema.yaml` with a changelog
comment there — per-column version history was considered and rejected as over-engineering for
this workload (occasional bulk reads by a stable key, not point lookups or streaming).

### 3. Schema contract: `sql/table_schema.yaml` replaces `sql/schema.sql`

`sql/schema.sql` was confirmed drifted during the README rewrite (lists tables the build script
doesn't create — `funguild`, `chrom_info`, `mmseqs_orthogroup_clusters`,
`mmseqs_orthogroup_cluster_count`, `prosite`, `pfam_UoT`, `gene_pairwise_distances` — and is
missing several it does — `telomere_summary`, `telomere_tracts`, `busco_genome`, `wolfpsort`,
`predgpi`, `gene_CDS`, `gene_intergenic_distances`, `wgd_ks`, `wgd_ks_summary`, and all 3 views).
Replace it with one machine-readable file, per table:

```yaml
asm_stats:
  version: 1
  key: ASMID
  source_parquet: tables/asm_stats.parquet
  join: {table: species, on: ASMID, how: left}
  columns:
    - {name: ASMID, type: VARCHAR}
    - {name: LOCUSTAG, type: VARCHAR}       # from the species join
    - {name: TOTAL_LENGTH, type: BIGINT}    # aliased from total_length_bp
    # ...
```

This is load-bearing, used to (a) drive the typed `COPY ... (FORMAT PARQUET)` at merge time (no
more `read_csv_auto` inference drift between runs), (b) generate the `CREATE VIEW` catalog SQL,
(c) generate human-readable docs, (d) a lint step that runs `DESCRIBE` on each real Parquet file
and diffs it against this YAML to catch drift automatically going forward, rather than by manual
audit.

### 4. Atomic writes + integrity guards

**Atomic writes**: `MERGE_*` currently publishes via `publishDir mode: 'copy'`, which deletes an
existing target then streams the new file — not atomic (confirmed: Nextflow 26.04.6's copy
publisher). A concurrent `read_parquet()` during that window sees `ENOENT` or (per Fable, tested
against DuckDB 1.1.3) a loud `No magic bytes found at end of file` error — a transient failure,
not silent corruption, but avoidable. Storage is confirmed GPFS (not NFS), where `rename(2)` is
atomic. Fix: each `MERGE_*` writes to `tables/.tmp.<name>.<task_hash>.parquet` then `mv -f` into
place; the manifest sidecar is written the same way, sequenced after the Parquet rename
completes. The Nextflow work-dir copy is untouched (so `-resume` still works) — only the final
publish step changes.

**Duplicate-key guard**: neither `summarize_busco_stats.py` nor `summarize_asm_stats.py`
de-duplicates or fails on a repeated `ASMID` (confirmed — no `duplicated()`/`GROUP BY`/`DISTINCT`
anywhere in either). Under the views-only design there's no `CREATE TABLE ... UNIQUE INDEX` step
left to catch this. Add a shared `assert_unique_key(df, "ASMID")` to both, raising/exiting 1
before the Parquet write, matching this repo's fail-loud convention. Fable found a concrete path
to a real duplicate (not just theoretical): `BUSCO_GENOME`'s `storeDir` output is keyed by
`{ASMID}.BUSCO_summary.{lineage}.txt` — if a genome's lineage changes, the old file is never
cleaned up and both get merged. 0 duplicates on disk today (23,692 files / 23,125 rows — that
gap itself is unrelated staleness, not duplication), so this guard is preventive.

Additional assertions (Fable): `row_count > 0`, and for `ASMID`-keyed tables, `row_count >=`
the previous manifest's `row_count` unless an explicit override flag is passed — catches an
accidental partial-glob regression.

**Lower priority, same pattern**: `calculate_intergenic.py` and `build_genestats_table.py` only
`print()` a warning on a duplicate gene ID, never fail or dedup. Recommended upgrade to
fail-loud, lower priority than the `ASMID` guard since it's a per-gene key, not the genome-level
join key DeltaGain and the manifest depend on.

## `nextflow/bin` cleanup (audited separately, folded in here)

- `git rm nextflow/bin/collect_chrom_info.py nextflow/bin/rewrite_abinitio_param_paths.py` —
  both confirmed dead (zero references in any `.nf`/module/config, not documented as
  intentionally-retained tools the way `setup_singularity_cache.sh` is). Missed by the earlier
  `asm_reports`/`collect_asm_stats.py` cleanup pass.
- Every other script in `nextflow/bin/` (57 checked) writes an intermediate (stdout/csv/tsv)
  that its calling Nextflow module converts to Parquet via a `duckdb ... COPY` step — the
  established, correct pattern. No format-migration work needed beyond the guards above.
- `pick_representative_strain.py` already reads `tables/*.parquet` directly with zero
  `db/BFD.duckdb` references — already fully compatible, confirms decision F (2026-08-03) was
  carried through completely.
- `species_reuse_clusters.py` queries a separate `ani.duckdb` (from `compare_ANI.nf`) —
  unrelated to this redesign.

## Not risks (confirmed, not acted on)

- Per-query Parquet schema inference: a footer read, microseconds at this scale.
- DuckDB version skew between the CLI (1.1.3) building the catalog and the Python DuckDB (1.5.3)
  the MCP server runs: Parquet + views are version-neutral, and trivial-to-rebuild views make a
  mismatch a non-issue if it ever occurs — actually a benefit of dropping materialization.
- Stray files already in `tables/` (`telomeres.parquet.bak.*`, `BUSCO.csv.gz`): harmless once
  the catalog is generated from the YAML allowlist rather than a glob; add a lint warning for
  unrecognized files rather than treating it as blocking.

## Out of scope for this design

- Any change to `scripts/extract_bfd_taxonomic_subset.py` (still produces a real materialized,
  taxon-scoped DB on explicit request — unaffected).
- Materializing `v_species_summary` (which scans every gene table per query) — currently
  invisible cost with 1 annotated genome in the live data; revisit once gene tables are back to
  full scale, per Fable's "measure before deciding."
- Partitioning `tables/*.parquet` into multiple files per type — T-014 explicitly rejected this
  (small-file metadata cost), and nothing here reopens it.

## Handoff

A companion short doc answering DeltaGain's stated requirements point-by-point is at
`DeltaGain/docs/superpowers/notes/busco-n50-datalake-response.md`.
