# BFD Views-Only Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate `species`, `asm_stats`, and `busco_genome` — the three tables DeltaGain and the
representative-strain-selection tooling actually depend on — to the views-only, manifest-backed,
atomically-written design, proving the mechanism end-to-end before the remaining ~27 tables get
the same treatment in a fast-follow plan.

**Architecture:** A schema contract (`sql/table_schema.json`) drives both a codegen script that
emits `CREATE VIEW` SQL for `db/BFD.duckdb` (with explicit `CAST`s, not `SELECT *`) and a shared
publish helper (`nextflow/bin/publish_table.py`) that every migrated `MERGE_*` module calls to
atomically publish its Parquet file plus a per-table JSON freshness manifest. A shared
`assert_unique_key` check closes the duplicate-key gaps found in `summarize_busco_stats.py`,
`summarize_asm_stats.py`, and `build_species_table.py`.

**Tech Stack:** Python 3 stdlib only (deliberately no PyYAML — see revision note below), DuckDB
CLI (already `module load`ed throughout this pipeline), Nextflow DSL2.

**Spec:** `docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md`

## Revision note (2026-09-09, before any task was implemented)

This plan was reviewed twice before implementation started: once by a data-design/database
persona (Fable), once by a Python-code reviewer that actually materialized and ran the plan's
code in a scratch repo rather than only reading it. Both found real, concrete defects in the
first draft. All are fixed in the version below:

- **Switched the schema contract from YAML to JSON.** The Python reviewer executed the plan's own
  code and found PyYAML's `safe_load` (YAML 1.1 semantics) parses the bare dict key `on` in
  `join: {table: species, on: ASMID, how: left}` as the *boolean* `True`, not the string `"on"` —
  `join["on"]` then raises `KeyError` inside `generate_bfd_catalog_sql.py`, breaking both
  join-based views outright. Separately, Fable confirmed PyYAML isn't even installed in the
  `merge` Nextflow label's Python (`module load biopython` → no `yaml` module). JSON has neither
  problem and needs no new dependency.
- **Fixed the `species` schema entry** (Fable, confirmed against the real `build_species_table.py`
  and live `tables/species.parquet`): the draft included a `TRANSL_TABLE` column that
  `build_species_table.py` never actually writes (it's in `samples.csv` but not in that script's
  `REQUIRED_COLUMNS`/`OPTIONAL_COLUMNS`), and typed `NCBI_TAXONID` as `VARCHAR` when the live file
  has it as `BIGINT`. Both fixed below.
- **Added explicit `CAST` in the generated `CREATE VIEW` SQL.** The draft's codegen read a `type`
  field from the schema but never used it — no actual type enforcement, so the design's stated
  protection against silent type drift didn't exist in the code. Fixed in Task 4.
- **Added existence guards to the codegen.** Fable confirmed `CREATE VIEW` binds (and errors)
  eagerly in DuckDB 1.1.3 if the underlying Parquet file or `read_json_auto` glob doesn't exist —
  and `busco_genome.parquet` won't exist on a fresh run since `run_busco_genome` defaults `false`,
  nor will `tables/_manifest/` before this plan's changes have ever run once. Fixed in Task 4:
  skip (with a `-- SKIP:` comment) any table or the manifest view whose source doesn't exist yet.
- **Fixed absolute-path resolution for `build_BFD_duckDB.sh`.** The draft's `BIN_DIR`/`SCHEMA` env
  vars were never actually passed by the `BUILD_DUCKDB` module that invokes this script in
  production, and defaulted to paths relative to whatever directory Nextflow's work-dir happens
  to be — which is never the repo root. Fixed by adding a step to Task 10 that updates
  `BUILD_DUCKDB/main.nf` to pass both as absolute, `${projectDir}`-derived paths.
- **Fixed a real false-positive in the row-count-shrink check.** Fable traced actual pipeline
  semantics: `merge_all` defaults `false`, meaning `MERGE_ASM_STATS`/`MERGE_BUSCO_GENOME` normally
  see only the current run's (small) manifest, not a full glob — so the shrink check as
  originally written would reject nearly every default/`--taxon`/`--n_test` run. Fixed in Task 3:
  the check is now opt-in (`--enforce-no-shrink`), passed only when `params.merge_all` is true.
- **Added a duplicate-key guard to `build_species_table.py`** (not in the original scope note's
  exclusion list — Fable pointed out `species(ASMID)`/`species(LOCUSTAG)` uniqueness was enforced
  by the old materialized table's `UNIQUE INDEX` and is now enforced by nothing). Added as a step
  in Task 9.
- **Manifest `columns` now reflects the real Parquet file** (via `DESCRIBE`), not the schema
  contract's intended columns — otherwise the manifest can't actually detect the drift it exists
  to catch. Fixed in Task 3.
- **Closed a TOCTOU race** the Python reviewer found between the shrink-check read and the
  publish write (two concurrent `publish_table.py` runs for the same table could interleave).
  Fixed in Task 3 with a `flock`-based per-table lock (GPFS supports real POSIX locking, unlike
  NFS — confirmed by Fable's storage check during the design review).
- **Fixed `assert_unique_key`'s `KeyError`-on-missing-column** into a clean `sys.exit` message,
  and the generated-SQL test's fragile raw-CSV-substring assertions into proper CSV parsing (both
  flagged by the Python reviewer). Fixed in Tasks 2 and 4.
- **Confirmed NOT a bug** (Python reviewer tested directly against DuckDB 1.1.3): `t.ORDER`
  (alias-qualified) works fine; only a bare unqualified `SELECT ORDER` would fail. The codegen
  always qualifies column references, so no fix was needed there.

## Global Constraints

- Every new/modified Python script must fail loudly (non-zero exit, message to stderr) on bad
  input — never warn-and-continue. This repo's established convention (`.living/learnings.md`).
- `sql/table_schema.json` paths embedded in generated `CREATE VIEW` SQL must be **absolute** — a
  view is evaluated at query time by whatever process opens `db/BFD.duckdb` (e.g. the MCP server,
  possibly from a different working directory), so a relative path would silently break for any
  consumer not running from the repo root.
- No `${BASH_SOURCE[0]}`/`dirname` path resolution in any script that could run under SLURM —
  resolve paths via caller-supplied arguments/env vars instead (this repo's own documented
  HPCC gotcha).
- No new non-stdlib Python dependencies — the pipeline's various `module load`ed Python
  environments are not guaranteed to have anything beyond stdlib available (confirmed: no
  PyYAML in the `merge` label's environment).
- Match this repo's existing test style: standalone `python3 nextflow/tests/test_X.py` scripts
  with plain `assert` (see `nextflow/tests/test_hash_bucket_parity.py`) — no pytest, no fixtures.

## Scope note (read before starting)

This plan covers `species`, `asm_stats`, `busco_genome`, plus the two confirmed dead scripts and
the two lower-priority warn→fail-loud upgrades. It deliberately does **not** touch:
- `samples.parquet` (still published via the old `publishDir mode:'copy'` path — not part of the
  catalog today, no consumer depends on its atomicity yet).
- Any of the other ~27 tables (`gene_*`, `pfam`, `cazy`, ...) — same mechanism, mechanical
  repetition, tracked as a separate fast-follow plan once this one is proven working.
- `subset_samples.py`'s duplicate-key guard — `samples.parquet` is out of scope per above, so its
  producer script is too.
- Fixing the `-resume`-doesn't-restore-a-deleted-published-file interaction that removing
  `publishDir` introduces (Fable, Task 7/8/9) — a real, disclosed semantic change (a cached task
  no longer re-invokes `publish_table.py` on resume), not something this plan solves. Recovery if
  it ever matters: manually re-run `publish_table.py` against the cached task's work-dir file
  (found via `nextflow log -f workdir`).

---

### Task 1: Schema contract — `sql/table_schema.json`

**Files:**
- Create: `sql/table_schema.json`

**Interfaces:**
- Produces: a JSON mapping `{table_name: {version, key, source_parquet, join?, columns: [{name, type, from?}], changelog}}`, loaded by `nextflow/bin/bfd_common.load_schema()` (Task 2) and consumed by `nextflow/bin/generate_bfd_catalog_sql.py` (Task 4) and `nextflow/bin/publish_table.py` (Task 3).

- [ ] **Step 1: Write the file**

```json
{
  "species": {
    "version": 1,
    "key": "LOCUSTAG",
    "source_parquet": "species.parquet",
    "columns": [
      {"name": "ASMID", "type": "VARCHAR"},
      {"name": "LOCUSTAG", "type": "VARCHAR"},
      {"name": "PHYLUM", "type": "VARCHAR"},
      {"name": "SUBPHYLUM", "type": "VARCHAR"},
      {"name": "CLASS", "type": "VARCHAR"},
      {"name": "ORDER", "type": "VARCHAR"},
      {"name": "FAMILY", "type": "VARCHAR"},
      {"name": "GENUS", "type": "VARCHAR"},
      {"name": "SPECIES", "type": "VARCHAR"},
      {"name": "BUSCO_LINEAGE", "type": "VARCHAR"},
      {"name": "SUBCLASS", "type": "VARCHAR"},
      {"name": "NCBI_TAXONID", "type": "BIGINT"},
      {"name": "STRAIN", "type": "VARCHAR"}
    ],
    "changelog": [
      "v1 (2026-09-09): matches build_species_table.py's REQUIRED_COLUMNS + OPTIONAL_COLUMNS output exactly (confirmed against live tables/species.parquet during design review -- do not add TRANSL_TABLE, samples.csv has it but build_species_table.py never copies it)"
    ]
  },
  "asm_stats": {
    "version": 1,
    "key": "ASMID",
    "source_parquet": "asm_stats.parquet",
    "join": {"table": "species", "on": "ASMID", "how": "left"},
    "columns": [
      {"name": "LOCUSTAG", "type": "VARCHAR", "from": "species.LOCUSTAG"},
      {"name": "ASMID", "type": "VARCHAR"},
      {"name": "contig_count", "type": "BIGINT"},
      {"name": "TOTAL_LENGTH", "type": "BIGINT", "from": "total_length_bp"},
      {"name": "min_contig_bp", "type": "BIGINT"},
      {"name": "max_contig_bp", "type": "BIGINT"},
      {"name": "median_contig_bp", "type": "DOUBLE"},
      {"name": "mean_contig_bp", "type": "DOUBLE"},
      {"name": "L50", "type": "BIGINT"},
      {"name": "N50_bp", "type": "BIGINT"},
      {"name": "L90", "type": "BIGINT"},
      {"name": "N90_bp", "type": "BIGINT"},
      {"name": "GC_PERCENT", "type": "DOUBLE", "from": "gc_pct"},
      {"name": "n_gap_count", "type": "BIGINT"},
      {"name": "total_n_bases", "type": "BIGINT"},
      {"name": "masked_bases", "type": "BIGINT"},
      {"name": "masked_pct", "type": "DOUBLE"},
      {"name": "t2t_scaffolds", "type": "BIGINT"},
      {"name": "telomere_fwd", "type": "BIGINT"},
      {"name": "telomere_rev", "type": "BIGINT"}
    ],
    "changelog": [
      "v1 (2026-09-09): matches build_BFD_duckDB.sh's pre-migration asm_stats CREATE TABLE"
    ]
  },
  "busco_genome": {
    "version": 1,
    "key": "ASMID",
    "source_parquet": "busco_genome.parquet",
    "join": {"table": "species", "on": "ASMID", "how": "left"},
    "columns": [
      {"name": "LOCUSTAG", "type": "VARCHAR", "from": "species.LOCUSTAG"},
      {"name": "ASMID", "type": "VARCHAR"},
      {"name": "complete_pct", "type": "DOUBLE"},
      {"name": "single_pct", "type": "DOUBLE"},
      {"name": "duplicated_pct", "type": "DOUBLE"},
      {"name": "fragmented_pct", "type": "DOUBLE"},
      {"name": "missing_pct", "type": "DOUBLE"},
      {"name": "n_markers", "type": "BIGINT"},
      {"name": "lineage", "type": "VARCHAR"}
    ],
    "changelog": [
      "v1 (2026-09-09): initial schema-contract entry"
    ]
  }
}
```

(Field-level documentation lives in this plan and the design doc, not inline in the JSON, which
can't carry comments — `version` bumps on any column add/rename/remove, logged in that table's
`changelog`; `key` is the column `assert_unique_key()` enforces before write; `join.how` is a SQL
join type used verbatim; a column's `from` is either a same-table source name being renamed, or
an already-qualified `"<jointable>.<col>"` reference for a joined-in column — omit `from` when the
output name matches the source name.)

- [ ] **Step 2: Validate it parses and has the expected shape**

Run: `python3 -c "import json; d = json.load(open('sql/table_schema.json')); assert set(d) == {'species','asm_stats','busco_genome'}, d.keys(); assert d['asm_stats']['join']['on'] == 'ASMID'; print('OK', list(d))"`
Expected: `OK ['species', 'asm_stats', 'busco_genome']`

- [ ] **Step 3: Commit**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
git add sql/table_schema.json
git commit -m "Add schema contract for species/asm_stats/busco_genome

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 2: `nextflow/bin/bfd_common.py` — shared schema loader + duplicate-key guard

**Files:**
- Create: `nextflow/bin/bfd_common.py`
- Test: `nextflow/tests/test_bfd_common.py`

**Interfaces:**
- Produces: `load_schema(path) -> dict`, `assert_unique_key(rows: list[dict], key_field: str, table_name: str) -> None` (exits 1 on duplicate or on a row missing `key_field`). Consumed by Tasks 3, 4, 5, 6.

- [ ] **Step 1: Write the failing test**

```python
#!/usr/bin/env python3
"""test_bfd_common.py -- unit tests for nextflow/bin/bfd_common.py.

Usage: python3 nextflow/tests/test_bfd_common.py
"""
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent / "bin"))

from bfd_common import assert_unique_key, load_schema  # noqa: E402


def test_load_schema_reads_real_contract():
    repo_root = TESTS_DIR.parent.parent
    schema = load_schema(str(repo_root / "sql" / "table_schema.json"))
    assert set(schema) == {"species", "asm_stats", "busco_genome"}, schema.keys()
    assert schema["asm_stats"]["key"] == "ASMID"
    assert schema["asm_stats"]["join"] == {"table": "species", "on": "ASMID", "how": "left"}


def test_assert_unique_key_passes_on_unique_rows():
    rows = [{"ASMID": "A1", "x": 1}, {"ASMID": "A2", "x": 2}]
    assert_unique_key(rows, "ASMID", "test_table")  # must not raise/exit


def _run_snippet(code):
    bin_dir = str(TESTS_DIR.parent / "bin")
    return subprocess.run(
        [sys.executable, "-c", f"import sys; sys.path.insert(0, {bin_dir!r}); {code}"],
        capture_output=True, text=True,
    )


def test_assert_unique_key_fails_loud_on_duplicate():
    result = _run_snippet(
        "from bfd_common import assert_unique_key; "
        "assert_unique_key([{'ASMID':'A1'},{'ASMID':'A1'}], 'ASMID', 'test_table')"
    )
    assert result.returncode != 0, "expected non-zero exit on duplicate key"
    assert "duplicate" in result.stderr.lower(), result.stderr


def test_assert_unique_key_fails_loud_on_missing_key_field():
    result = _run_snippet(
        "from bfd_common import assert_unique_key; "
        "assert_unique_key([{'OTHER':'x'}], 'ASMID', 'test_table')"
    )
    assert result.returncode != 0, "expected non-zero exit, not an uncaught KeyError"
    assert "missing key field" in result.stderr.lower(), result.stderr
    assert "Traceback" not in result.stderr, "should be a clean sys.exit, not a raw traceback"


if __name__ == "__main__":
    test_load_schema_reads_real_contract()
    test_assert_unique_key_passes_on_unique_rows()
    test_assert_unique_key_fails_loud_on_duplicate()
    test_assert_unique_key_fails_loud_on_missing_key_field()
    print("All tests passed.")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 nextflow/tests/test_bfd_common.py`
Expected: `ModuleNotFoundError: No module named 'bfd_common'`

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""bfd_common.py -- shared helpers for BFD merge/summarize/catalog scripts.

Used by summarize_busco_stats.py, summarize_asm_stats.py, publish_table.py,
and generate_bfd_catalog_sql.py. sql/table_schema.json is the schema
contract this module loads -- see
docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md.

Deliberately JSON, not YAML: PyYAML's YAML-1.1 boolean coercion silently
turns a bare `on` dict key into the boolean True (confirmed the hard way --
see this plan's "Revision note"), and PyYAML isn't installed in every
Python environment this pipeline's Nextflow labels use anyway.
"""
import json
import sys


def load_schema(path):
    """Load sql/table_schema.json into a dict keyed by table name."""
    with open(path) as fh:
        schema = json.load(fh)
    if not isinstance(schema, dict):
        sys.exit(f"ERROR: {path} did not parse to a mapping of table -> definition")
    return schema


def assert_unique_key(rows, key_field, table_name):
    """Fail loudly (exit 1) if `key_field` repeats across `rows`, or if any
    row is missing `key_field` entirely.

    `rows` is a list of dicts -- the in-memory row-accumulation pattern used
    by summarize_busco_stats.py/summarize_asm_stats.py before writing CSV.
    """
    seen = {}
    dupes = set()
    for row in rows:
        if key_field not in row:
            sys.exit(f"ERROR: {table_name}: row missing key field '{key_field}': {row}")
        key = row[key_field]
        if key in seen:
            dupes.add(key)
        else:
            seen[key] = row
    if dupes:
        sys.exit(
            f"ERROR: {table_name}: duplicate {key_field} value(s) found, "
            f"refusing to write ({len(dupes)} distinct duplicated key(s)). "
            f"First few: {sorted(dupes)[:5]}"
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 nextflow/tests/test_bfd_common.py`
Expected: `All tests passed.`

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/bfd_common.py nextflow/tests/test_bfd_common.py
git commit -m "Add bfd_common: schema loader + duplicate-key guard

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 3: `nextflow/bin/publish_table.py` — atomic publish + manifest sidecar

**Files:**
- Create: `nextflow/bin/publish_table.py`
- Test: `nextflow/tests/test_publish_table.py`

**Interfaces:**
- Consumes: `bfd_common.load_schema(path)` (Task 2).
- Produces: a CLI (`--table`, `--local-parquet`, `--tables-dir`, `--built-by`, `--merge-run-id`, `--schema`, `--enforce-no-shrink`, `--allow-shrink`) that atomically publishes `<tables-dir>/<table>.parquet` and `<tables-dir>/_manifest/<table>.json`, serialized per-table via `flock`. Consumed by Tasks 7, 8, 9 (Nextflow module changes).

- [ ] **Step 1: Write the failing test**

```python
#!/usr/bin/env python3
"""test_publish_table.py -- unit tests for nextflow/bin/publish_table.py.

Usage: python3 nextflow/tests/test_publish_table.py
Requires: `duckdb` CLI on PATH (module load duckdb on HPCC).
"""
import csv
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
BIN_DIR = TESTS_DIR.parent / "bin"
REPO_ROOT = TESTS_DIR.parent.parent
SCHEMA_PATH = str(REPO_ROOT / "sql" / "table_schema.json")


def make_fixture_parquet(path, rows_tsv):
    tsv = path.with_suffix(".tsv")
    tsv.write_text(rows_tsv)
    subprocess.run(
        ["duckdb", "-c",
         f"COPY (SELECT * FROM read_csv_auto('{tsv}', delim='\\t', sample_size=-1)) "
         f"TO '{path}' (FORMAT PARQUET);"],
        check=True, capture_output=True,
    )


def run_publish(local, tables_dir, run_id, extra_args=()):
    return subprocess.run(
        [sys.executable, str(BIN_DIR / "publish_table.py"),
         "--table", "busco_genome",
         "--local-parquet", str(local),
         "--tables-dir", str(tables_dir),
         "--built-by", "TEST",
         "--merge-run-id", run_id,
         "--schema", SCHEMA_PATH,
         *extra_args],
        capture_output=True, text=True,
    )


def test_publish_writes_parquet_and_manifest_atomically():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        local = td / "busco_genome.parquet"
        make_fixture_parquet(local, "ASMID\tcomplete_pct\nGCA_1\t99.0\n")
        tables_dir = td / "tables"
        tables_dir.mkdir()

        result = run_publish(local, tables_dir, "test-run-1")
        assert result.returncode == 0, result.stderr

        published = tables_dir / "busco_genome.parquet"
        manifest_file = tables_dir / "_manifest" / "busco_genome.json"
        assert published.is_file()
        assert manifest_file.is_file()
        manifest = json.loads(manifest_file.read_text())
        assert manifest["table"] == "busco_genome"
        assert manifest["row_count"] == 1
        assert manifest["built_by"] == "TEST"
        assert manifest["merge_run_id"] == "test-run-1"
        assert manifest["schema_version"] == 1
        assert "ASMID" in manifest["columns"]  # derived from the real Parquet file via DESCRIBE
        assert len(manifest["parquet_sha256"]) == 64

        leftovers = list(tables_dir.glob(".tmp.*")) + list((tables_dir / "_manifest").glob(".tmp.*"))
        assert leftovers == [], leftovers


def test_publish_rejects_empty_parquet():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        local = td / "busco_genome.parquet"
        make_fixture_parquet(local, "ASMID\tcomplete_pct\n")  # header only, 0 rows
        tables_dir = td / "tables"
        tables_dir.mkdir()

        result = run_publish(local, tables_dir, "test-run-2")
        assert result.returncode != 0
        assert "0 rows" in result.stderr or "0 rows" in result.stdout


def test_publish_shrink_check_is_opt_in_and_overridable():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        tables_dir = td / "tables"
        tables_dir.mkdir()

        big = td / "big.parquet"
        make_fixture_parquet(big, "ASMID\tcomplete_pct\nA1\t1.0\nA2\t2.0\n")
        small = td / "small.parquet"
        make_fixture_parquet(small, "ASMID\tcomplete_pct\nA1\t1.0\n")

        assert run_publish(big, tables_dir, "run-a", ["--enforce-no-shrink"]).returncode == 0

        # Without --enforce-no-shrink (the default, matching merge_all=false runs),
        # a smaller table is accepted -- this is normal behavior for a
        # current-run-only merge, not a data-integrity violation.
        result_default = run_publish(small, tables_dir, "run-b")
        assert result_default.returncode == 0, result_default.stderr

        # Re-publish the big one, then confirm --enforce-no-shrink actually
        # rejects a real shrink when asked to.
        assert run_publish(big, tables_dir, "run-c", ["--enforce-no-shrink"]).returncode == 0
        result_shrink = run_publish(small, tables_dir, "run-d", ["--enforce-no-shrink"])
        assert result_shrink.returncode != 0

        result_override = run_publish(
            small, tables_dir, "run-e", ["--enforce-no-shrink", "--allow-shrink"]
        )
        assert result_override.returncode == 0, result_override.stderr


if __name__ == "__main__":
    test_publish_writes_parquet_and_manifest_atomically()
    test_publish_rejects_empty_parquet()
    test_publish_shrink_check_is_opt_in_and_overridable()
    print("All tests passed.")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 nextflow/tests/test_publish_table.py`
Expected: `FileNotFoundError` (`publish_table.py` doesn't exist yet)

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""publish_table.py -- atomically publish a MERGE_* Parquet output into
tables/, and write a matching per-table freshness-manifest sidecar into
tables/_manifest/.

Usage:
    python3 publish_table.py --table busco_genome \\
        --local-parquet busco_genome.parquet \\
        --tables-dir /abs/path/tables \\
        --built-by MERGE_BUSCO_GENOME \\
        --merge-run-id <workflow.sessionId> \\
        --schema /abs/path/sql/table_schema.json \\
        [--enforce-no-shrink] [--allow-shrink]

Writes (in this order, under a per-table flock so two concurrent runs for
the SAME table can't interleave their read-check-write sequence):
    <tables-dir>/<table>.parquet          (atomic rename into place)
    <tables-dir>/_manifest/<table>.json   (atomic rename, written AFTER
                                            the parquet rename above)

Fails loudly (non-zero exit) if the local Parquet file is missing, empty,
or has 0 rows. --enforce-no-shrink additionally fails if the row count
dropped versus the previous manifest -- pass it only when this run's
manifest is known to represent the FULL current dataset (i.e.
params.merge_all == true), never on a current-run-only partial merge,
where a smaller table is expected and correct. --allow-shrink overrides
--enforce-no-shrink for a deliberate, known-smaller republish (e.g. after
suppress.txt removes genomes).
"""
import argparse
import fcntl
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bfd_common import load_schema  # noqa: E402


def sha256_file(path, chunk_size=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def duckdb_row_count(parquet_path):
    out = subprocess.run(
        ["duckdb", "-csv", "-noheader", "-c",
         f"SELECT COUNT(*) FROM read_parquet('{parquet_path}');"],
        capture_output=True, text=True, check=True,
    )
    return int(out.stdout.strip())


def duckdb_describe_columns(parquet_path):
    """Real column names of the Parquet file being published, via DESCRIBE --
    NOT copied from the schema contract, so the manifest can actually reveal
    drift between what the contract expects and what got written."""
    out = subprocess.run(
        ["duckdb", "-csv", "-c", f"DESCRIBE SELECT * FROM read_parquet('{parquet_path}');"],
        capture_output=True, text=True, check=True,
    )
    lines = out.stdout.strip().splitlines()[1:]  # skip "column_name,column_type,..." header
    return [line.split(",")[0] for line in lines if line]


def atomic_write_bytes_from_file(src_path, final_path):
    final_path = Path(final_path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = final_path.parent / f".tmp.{final_path.name}.{os.getpid()}.{int(time.time() * 1000)}"
    with open(src_path, "rb") as src, open(tmp_path, "wb") as dst:
        while True:
            chunk = src.read(1 << 20)
            if not chunk:
                break
            dst.write(chunk)
    os.replace(tmp_path, final_path)  # atomic rename, same filesystem (GPFS)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--table", required=True)
    p.add_argument("--local-parquet", required=True)
    p.add_argument("--tables-dir", required=True)
    p.add_argument("--built-by", required=True)
    p.add_argument("--merge-run-id", required=True)
    p.add_argument("--schema", required=True)
    p.add_argument("--enforce-no-shrink", action="store_true",
                    help="only pass this on a full/merge_all=true run")
    p.add_argument("--allow-shrink", action="store_true")
    args = p.parse_args()

    local_parquet = Path(args.local_parquet)
    if not local_parquet.is_file() or local_parquet.stat().st_size == 0:
        sys.exit(f"ERROR: {local_parquet} is missing or empty; refusing to publish {args.table}")

    schema = load_schema(args.schema)
    entry = schema.get(args.table)
    if entry is None:
        sys.exit(f"ERROR: table '{args.table}' has no entry in {args.schema}")

    row_count = duckdb_row_count(str(local_parquet))
    if row_count == 0:
        sys.exit(f"ERROR: {local_parquet} has 0 rows; refusing to publish {args.table}")

    tables_dir = Path(args.tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir = tables_dir / "_manifest"
    final_parquet = tables_dir / f"{args.table}.parquet"
    final_manifest = manifest_dir / f"{args.table}.json"

    # Serialize the whole read-check-write sequence per table: two concurrent
    # publishes for the same table (e.g. an accidental overlapping resume)
    # must not interleave, or the final manifest could describe a different
    # parquet file than the one actually on disk.
    lock_path = tables_dir / f".lock.{args.table}"
    with open(lock_path, "w") as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)

        if args.enforce_no_shrink and not args.allow_shrink and final_manifest.is_file():
            try:
                prev_count = json.loads(final_manifest.read_text()).get("row_count", 0)
            except (json.JSONDecodeError, OSError):
                prev_count = 0
            if row_count < prev_count:
                sys.exit(
                    f"ERROR: {args.table} row_count dropped from {prev_count} to {row_count} "
                    f"-- refusing to publish (looks like a partial-glob regression); pass "
                    f"--allow-shrink to override if this is expected"
                )

        checksum = sha256_file(local_parquet)
        columns = duckdb_describe_columns(str(local_parquet))

        atomic_write_bytes_from_file(local_parquet, final_parquet)

        manifest = {
            "table": args.table,
            "path": f"tables/{args.table}.parquet",
            "row_count": row_count,
            "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "built_by": args.built_by,
            "merge_run_id": args.merge_run_id,
            "schema_version": entry["version"],
            "columns": columns,
            "parquet_sha256": checksum,
        }
        manifest_dir.mkdir(parents=True, exist_ok=True)
        tmp_manifest = manifest_dir / f".tmp.{args.table}.json.{os.getpid()}.{int(time.time() * 1000)}"
        tmp_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        os.replace(tmp_manifest, final_manifest)

        fcntl.flock(lock_fh, fcntl.LOCK_UN)

    print(f"Published {args.table}: {row_count} rows -> {final_parquet}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 nextflow/tests/test_publish_table.py`
Expected: `All tests passed.`

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/publish_table.py nextflow/tests/test_publish_table.py
git commit -m "Add publish_table.py: atomic Parquet publish + freshness manifest

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 4: `nextflow/bin/generate_bfd_catalog_sql.py` — codegen `CREATE VIEW` SQL

**Files:**
- Create: `nextflow/bin/generate_bfd_catalog_sql.py`
- Test: `nextflow/tests/test_generate_bfd_catalog_sql.py`

**Interfaces:**
- Consumes: `bfd_common.load_schema(path)` (Task 2).
- Produces: a CLI (`--schema`, `--tables-dir`) printing `CREATE VIEW` SQL to stdout — one typed, `CAST`-enforced view per schema table (skipped with a `-- SKIP:` comment if its source Parquet, or a joined table's source Parquet, doesn't exist yet) plus `_table_manifest` (skipped the same way if no manifest JSON exists yet). Consumed by Task 10 (`build_BFD_duckDB.sh`).

- [ ] **Step 1: Write the failing test**

```python
#!/usr/bin/env python3
"""test_generate_bfd_catalog_sql.py -- generates the catalog SQL against a
tiny fixture Parquet dataset, feeds it into an in-memory duckdb, and checks
the resulting views return the expected joined/aliased/typed rows -- and
that a missing source file is skipped instead of crashing the whole build.

Usage: python3 nextflow/tests/test_generate_bfd_catalog_sql.py
Requires: `duckdb` CLI on PATH.
"""
import csv
import io
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
BIN_DIR = TESTS_DIR.parent / "bin"
REPO_ROOT = TESTS_DIR.parent.parent
SCHEMA_PATH = str(REPO_ROOT / "sql" / "table_schema.json")


def write_parquet(path, tsv_text):
    tsv = path.with_suffix(".tsv")
    tsv.write_text(tsv_text)
    subprocess.run(
        ["duckdb", "-c",
         f"COPY (SELECT * FROM read_csv_auto('{tsv}', delim='\\t', sample_size=-1)) "
         f"TO '{path}' (FORMAT PARQUET);"],
        check=True, capture_output=True,
    )


def query_csv_rows(db_path, sql):
    result = subprocess.run(["duckdb", "-csv", str(db_path), "-c", sql],
                             capture_output=True, text=True, check=True)
    reader = csv.DictReader(io.StringIO(result.stdout))
    return list(reader)


def test_catalog_views_join_cast_and_skip_missing_correctly():
    with tempfile.TemporaryDirectory() as td:
        tables_dir = Path(td) / "tables"
        tables_dir.mkdir()
        (tables_dir / "_manifest").mkdir()

        write_parquet(tables_dir / "species.parquet",
                      "LOCUSTAG\tASMID\tSPECIES\tSTRAIN\tPHYLUM\tSUBPHYLUM\tCLASS\tSUBCLASS\tORDER\tFAMILY\tGENUS\tNCBI_TAXONID\tBUSCO_LINEAGE\n"
                      "LOC1\tGCA_1\tfoo\tbar\tAscomycota\tNA\tNA\tNA\tNA\tNA\tFusarium\t123\tfungi_odb10\n")
        write_parquet(tables_dir / "asm_stats.parquet",
                      "ASMID\tcontig_count\ttotal_length_bp\tmin_contig_bp\tmax_contig_bp\tmedian_contig_bp\tmean_contig_bp\tL50\tN50_bp\tL90\tN90_bp\tgc_pct\tn_gap_count\ttotal_n_bases\tmasked_bases\tmasked_pct\tt2t_scaffolds\ttelomere_fwd\ttelomere_rev\n"
                      "GCA_1\t10\t1000000\t500\t200000\t50000\t100000\t3\t150000\t8\t50000\t48.5\t0\t0\t10000\t1.0\t2\t2\t2\n")
        # busco_genome.parquet deliberately NOT created -- run_busco_genome
        # defaults false in production, so its view must be skipped, not
        # crash the whole catalog build.
        (tables_dir / "_manifest" / "asm_stats.json").write_text(
            '{"table": "asm_stats", "row_count": 1, "built_at": "2026-09-09T00:00:00Z", '
            '"built_by": "TEST", "merge_run_id": "r1", "schema_version": 1, '
            '"columns": ["ASMID"], "parquet_sha256": "x"}'
        )

        gen = subprocess.run(
            [sys.executable, str(BIN_DIR / "generate_bfd_catalog_sql.py"),
             "--schema", SCHEMA_PATH, "--tables-dir", str(tables_dir)],
            capture_output=True, text=True, check=True,
        )
        assert "CREATE VIEW busco_genome" not in gen.stdout, gen.stdout
        assert "SKIP: busco_genome" in gen.stdout, gen.stdout

        db_path = Path(td) / "test.duckdb"
        subprocess.run(["duckdb", str(db_path)], input=gen.stdout, text=True, check=True)

        rows = query_csv_rows(db_path, "SELECT LOCUSTAG, ASMID, TOTAL_LENGTH, GC_PERCENT FROM asm_stats;")
        assert len(rows) == 1, rows
        assert rows[0]["LOCUSTAG"] == "LOC1"
        assert rows[0]["ASMID"] == "GCA_1"
        assert rows[0]["TOTAL_LENGTH"] == "1000000"
        assert rows[0]["GC_PERCENT"] == "48.5"

        manifest_rows = query_csv_rows(db_path, 'SELECT "table", row_count FROM _table_manifest;')
        assert manifest_rows == [{"table": "asm_stats", "row_count": "1"}], manifest_rows


if __name__ == "__main__":
    test_catalog_views_join_cast_and_skip_missing_correctly()
    print("All tests passed.")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 nextflow/tests/test_generate_bfd_catalog_sql.py`
Expected: `FileNotFoundError` (`generate_bfd_catalog_sql.py` doesn't exist yet)

- [ ] **Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""generate_bfd_catalog_sql.py -- emit CREATE VIEW statements for
db/BFD.duckdb from sql/table_schema.json: one typed, CAST-enforced VIEW per
table over read_parquet(tables/<t>.parquet) (joined per the schema's `join`
entry), plus a _table_manifest VIEW over the per-table JSON manifest
sidecars. A table (or the manifest view) whose source file(s) don't exist
yet is skipped with a `-- SKIP:` comment instead of crashing the whole
build -- DuckDB's CREATE VIEW binds eagerly and errors immediately if the
underlying Parquet/glob is missing, and several tables (e.g. busco_genome,
since run_busco_genome defaults false) legitimately won't exist on a fresh
or partial run.

Usage:
    python3 generate_bfd_catalog_sql.py --schema sql/table_schema.json \\
        --tables-dir /abs/path/tables > catalog.sql
    duckdb db/BFD.duckdb < catalog.sql
"""
import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bfd_common import load_schema  # noqa: E402


def quote(path):
    return path.replace("'", "''")


def build_view_sql(table_name, entry, schema, tables_dir):
    source_path = f"{tables_dir}/{entry['source_parquet']}"
    join = entry.get("join")

    if not os.path.isfile(source_path):
        return None, f"-- SKIP: {table_name} (source parquet not found: {source_path})"
    if join:
        joined_entry = schema[join["table"]]
        joined_path = f"{tables_dir}/{joined_entry['source_parquet']}"
        if not os.path.isfile(joined_path):
            return None, f"-- SKIP: {table_name} (join target parquet not found: {joined_path})"

    select_cols = []
    for col in entry["columns"]:
        col_type = col["type"]
        src_expr = col.get("from")
        if src_expr and "." in src_expr:
            raw = src_expr
        elif src_expr:
            raw = f"t.{src_expr}"
        else:
            raw = f"t.{col['name']}"
        select_cols.append(f"CAST({raw} AS {col_type}) AS {col['name']}")
    select_clause = ",\n    ".join(select_cols)

    if join:
        joined_path = f"{tables_dir}/{schema[join['table']]['source_parquet']}"
        on_col = join["on"]
        how = join["how"].upper()
        sql = (
            f"CREATE VIEW {table_name} AS\n"
            f"SELECT\n    {select_clause}\n"
            f"FROM read_parquet('{quote(source_path)}') AS t\n"
            f"{how} JOIN read_parquet('{quote(joined_path)}') AS {join['table']} "
            f"ON t.{on_col} = {join['table']}.{on_col};"
        )
    else:
        sql = (
            f"CREATE VIEW {table_name} AS\n"
            f"SELECT\n    {select_clause}\n"
            f"FROM read_parquet('{quote(source_path)}') AS t;"
        )
    return sql, None


def build_manifest_view_sql(tables_dir):
    manifest_dir = f"{tables_dir}/_manifest"
    if not glob.glob(f"{manifest_dir}/*.json"):
        return None, f"-- SKIP: _table_manifest (no JSON files found under {manifest_dir})"
    manifest_glob = f"{manifest_dir}/*.json"
    sql = f"CREATE VIEW _table_manifest AS\nSELECT * FROM read_json_auto('{quote(manifest_glob)}');"
    return sql, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema", required=True)
    ap.add_argument("--tables-dir", required=True)
    args = ap.parse_args()

    schema = load_schema(args.schema)
    tables_dir = args.tables_dir.rstrip("/")

    print("-- Auto-generated by generate_bfd_catalog_sql.py -- do not hand-edit.")
    print("-- Source of truth: sql/table_schema.json\n")
    for table_name, entry in schema.items():
        sql, skip_comment = build_view_sql(table_name, entry, schema, tables_dir)
        print(sql if sql else skip_comment)
        print()
    sql, skip_comment = build_manifest_view_sql(tables_dir)
    print(sql if sql else skip_comment)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 nextflow/tests/test_generate_bfd_catalog_sql.py`
Expected: `All tests passed.`

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/generate_bfd_catalog_sql.py nextflow/tests/test_generate_bfd_catalog_sql.py
git commit -m "Add generate_bfd_catalog_sql.py: schema-driven, type-safe CREATE VIEW codegen

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 5: Add `assert_unique_key` to `summarize_busco_stats.py`

**Files:**
- Modify: `nextflow/bin/summarize_busco_stats.py`

**Interfaces:**
- Consumes: `bfd_common.assert_unique_key(rows, key_field, table_name)` (Task 2).

- [ ] **Step 1: Find the row-accumulation point**

Run: `grep -n "rows = \[\]\|rows.append\|def main\|write" nextflow/bin/summarize_busco_stats.py`

Confirm there's a `rows = []` list built via `.append(...)` and a final write step (matches
Fable's confirmed finding — no dict-by-key, no dedup).

- [ ] **Step 2: Add the guard immediately before the write**

Add near the top of the file:
```python
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bfd_common import assert_unique_key  # noqa: E402
```
(add `from pathlib import Path` to the existing imports if not already present)

Immediately before the code that writes `rows` out to CSV/TSV, add:
```python
assert_unique_key(rows, "ASMID", "busco_genome")
```

- [ ] **Step 3: Verify manually with a duplicate-row fixture**

Construct a two-file manifest where both report the same ASMID (copy an existing
`*.BUSCO_summary.*.txt` fixture to a second lineage-suffixed name for the same ASMID) and run
the script against it; confirm it exits non-zero with a message containing "duplicate".

- [ ] **Step 4: Run the existing repo test suite to confirm no regression**

Run: `python3 nextflow/tests/test_hash_bucket_parity.py` (unrelated but the fastest existing
smoke check that the `nextflow/bin` import path still works)
Expected: passes as before.

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/summarize_busco_stats.py
git commit -m "summarize_busco_stats.py: fail loud on duplicate ASMID

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 6: Add `assert_unique_key` to `summarize_asm_stats.py`

**Files:**
- Modify: `nextflow/bin/summarize_asm_stats.py`

**Interfaces:**
- Consumes: `bfd_common.assert_unique_key(rows, key_field, table_name)` (Task 2).

- [ ] **Step 1: Find the row-accumulation point**

Run: `grep -n "rows = \[\]\|rows.append\|def main\|write" nextflow/bin/summarize_asm_stats.py`

- [ ] **Step 2: Add the guard, same pattern as Task 5**

```python
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bfd_common import assert_unique_key  # noqa: E402
```
Immediately before the write step:
```python
assert_unique_key(rows, "ASMID", "asm_stats")
```

- [ ] **Step 3: Verify manually with a duplicate-row fixture**, same approach as Task 5 Step 3,
using a duplicated `asm_stats` manifest entry.

- [ ] **Step 4: Run the existing repo test suite to confirm no regression**

Run: `python3 nextflow/tests/test_hash_bucket_parity.py`
Expected: passes as before.

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/summarize_asm_stats.py
git commit -m "summarize_asm_stats.py: fail loud on duplicate ASMID

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 7: `MERGE_BUSCO_GENOME/main.nf` — atomic publish via `publish_table.py`

**Files:**
- Modify: `nextflow/modules/BFD/MERGE_BUSCO_GENOME/main.nf`

**Interfaces:**
- Consumes: `nextflow/bin/publish_table.py` (Task 3).

- [ ] **Step 1: Remove `publishDir` and add the publish call to `script:`**

Replace the whole file with:
```groovy
include { tablesDir } from '../../common/utils.nf'

process MERGE_BUSCO_GENOME {
    label      'merge'

    input:
    path manifest

    output:
    path "busco_genome.parquet", emit: parquet

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/summarize_busco_stats.py \\
        --manifest ${manifest} \\
        -o         busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         busco_genome \\
        --local-parquet busco_genome.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_BUSCO_GENOME \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    printf 'ASMID\\tcomplete_pct\\tsingle_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\tn_markers\\tlineage\\nSTUB_ASM\\t99.0\\t98.0\\t1.0\\t0.5\\t0.5\\t758\\tfungi_odb10\\n' | gzip > busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         busco_genome \\
        --local-parquet busco_genome.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_BUSCO_GENOME \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
```

Two things changed from the pre-review draft: (1) the stub now emits one synthetic row instead of
a header-only (0-row) file, since `publish_table.py` rejects 0-row Parquet files; (2)
`--enforce-no-shrink` is only passed on the real `script:` block, and only when
`params.merge_all` is true — under the default `merge_all=false`, a run's manifest legitimately
covers only the current run's genomes, so a smaller table is expected and must not be rejected
(Fable's finding — see this plan's Revision note).

- [ ] **Step 2: Smoke-test the module in isolation via `-stub-run`**

Run: `cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_annotation_asco && nextflow run ../Fungi_BFD/nextflow/main.nf -c ../Fungi_BFD/nextflow/nextflow.config -profile test --pipeline BFD -stub-run -resume 2>&1 | tail -40`
Expected: run completes, `MERGE_BUSCO_GENOME` task shows `COMPLETED`, and
`tables/busco_genome.parquet` + `tables/_manifest/busco_genome.json` exist in the test profile's
output directory (`nextflow/tests/output/tables/`, per `test.config` — confirm this path in
`nextflow/conf/test.config` if unsure).

- [ ] **Step 3: Commit**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
git add nextflow/modules/BFD/MERGE_BUSCO_GENOME/main.nf
git commit -m "MERGE_BUSCO_GENOME: atomic publish via publish_table.py

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 8: `MERGE_ASM_STATS/main.nf` — atomic publish via `publish_table.py`

**Files:**
- Modify: `nextflow/modules/BFD/MERGE_ASM_STATS/main.nf`

**Interfaces:**
- Consumes: `nextflow/bin/publish_table.py` (Task 3).

- [ ] **Step 1: Remove `publishDir` and add the publish call**

Replace the whole file with:
```groovy
include { tablesDir } from '../../common/utils.nf'

process MERGE_ASM_STATS {
    label      'merge'

    input:
    path manifest
    path samples

    output:
    path "asm_stats.parquet", emit: parquet

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/summarize_asm_stats.py \\
        --manifest ${manifest} \\
        --samples  ${samples} \\
        -o         asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         asm_stats \\
        --local-parquet asm_stats.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_ASM_STATS \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    printf 'ASMID\\tSPECIES\\tSTRAIN\\tcontig_count\\ttotal_length_bp\\tmin_contig_bp\\tmax_contig_bp\\tmedian_contig_bp\\tmean_contig_bp\\tL50\\tN50_bp\\tL90\\tN90_bp\\tgc_pct\\tn_gap_count\\ttotal_n_bases\\tmasked_bases\\tmasked_pct\\tt2t_scaffolds\\ttelomere_fwd\\ttelomere_rev\\nSTUB_ASM\\tFoo bar\\tCBS 1\\t10\\t1000000\\t500\\t200000\\t50000\\t100000\\t3\\t150000\\t8\\t50000\\t48.5\\t0\\t0\\t10000\\t1.0\\t2\\t2\\t2\\n' | gzip > asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         asm_stats \\
        --local-parquet asm_stats.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_ASM_STATS \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
```

Same 0-row-stub fix and same `merge_all`-conditional `--enforce-no-shrink` as Task 7.

- [ ] **Step 2: Smoke-test via `-stub-run`**

Run: same command as Task 7 Step 2.
Expected: `MERGE_ASM_STATS` task `COMPLETED`, `tables/asm_stats.parquet` +
`tables/_manifest/asm_stats.json` present.

- [ ] **Step 3: Commit**

```bash
git add nextflow/modules/BFD/MERGE_ASM_STATS/main.nf
git commit -m "MERGE_ASM_STATS: atomic publish via publish_table.py

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 9: `MERGE_SAMPLES/main.nf` + `build_species_table.py` — atomic publish + dedup for `species.parquet`

**Files:**
- Modify: `nextflow/modules/BFD/MERGE_SAMPLES/main.nf`
- Modify: `nextflow/bin/build_species_table.py`

**Interfaces:**
- Consumes: `nextflow/bin/publish_table.py` (Task 3).
- Note: `samples.parquet` is intentionally left on the old `publishDir mode:'copy'` path — see
  "Scope note" at the top of this plan.

- [ ] **Step 1: Add a duplicate-key guard to `build_species_table.py`**

This script streams rows straight to the CSV writer rather than accumulating a list, so it can't
reuse `bfd_common.assert_unique_key` as-is; track seen keys inline instead. In the `main()` loop
(the `for row in reader:` block), change:
```python
for row in reader:
    if keep is not None and (row.get(args.key) or "").strip() not in keep:
        continue
    writer.writerow({c: (row.get(c) or "") for c in out_columns})
    n += 1
```
to:
```python
seen_keys = set()
seen_asmids = set()
for row in reader:
    if keep is not None and (row.get(args.key) or "").strip() not in keep:
        continue
    row_key = (row.get(args.key) or "").strip()
    row_asmid = (row.get("ASMID") or "").strip()
    if row_key in seen_keys:
        sys.exit(f"ERROR: duplicate {args.key} value '{row_key}' in {args.samples} -- "
                  f"refusing to write species table")
    if row_asmid in seen_asmids:
        sys.exit(f"ERROR: duplicate ASMID value '{row_asmid}' in {args.samples} -- "
                  f"refusing to write species table")
    seen_keys.add(row_key)
    seen_asmids.add(row_asmid)
    writer.writerow({c: (row.get(c) or "") for c in out_columns})
    n += 1
```
(`sys` is already imported at the top of this file — confirm with `grep -n "^import sys" nextflow/bin/build_species_table.py` before assuming so.)

- [ ] **Step 2: Restrict `publishDir` to `samples.parquet` only, publish `species.parquet` explicitly**

Replace the file's `process` block with:
```groovy
include { tablesDir } from '../../common/utils.nf'

// Per T-014 §D.2: always builds the single, full/unscoped master samples+species
// table, regardless of --taxon. taxonRowFilter() already restricts which genomes
// are *processed* under --taxon (at the base genome channel); MERGE_SAMPLES no
// longer also restricts the table it writes, so a --taxon run no longer produces
// a separate tables/<Taxon>/ subset -- see also the retired `matched` input this
// process used to take.
//
// species.parquet is published atomically with a freshness manifest (see
// docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md); samples.parquet
// stays on the old publishDir copy path -- it isn't part of the view catalog yet.
process MERGE_SAMPLES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy', pattern: 'samples.parquet'

    input:
    path(samples)

    output:
    path "samples.parquet", emit: samples
    path "species.parquet", emit: species

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${projectDir}/bin/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         species \\
        --local-parquet species.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_SAMPLES \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    python3 ${projectDir}/bin/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${projectDir}/bin/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         species \\
        --local-parquet species.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_SAMPLES \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
```

`subset_samples.py`/`build_species_table.py` already produce non-empty output against the test
profile's `samples.csv` fixture, so — unlike Tasks 7/8 — no stub-fixture row-count fix should be
needed here; verify this assumption in Step 3 and add a synthetic row to the stub only if it
turns out to produce 0 rows.

- [ ] **Step 3: Smoke-test via `-stub-run`, including the new dedup guard**

Run: same command as Task 7 Step 2.
Expected: `MERGE_SAMPLES` task `COMPLETED`, `tables/species.parquet` + `tables/_manifest/species.json`
present, `tables/samples.parquet` also present (via the unchanged `publishDir` path).

Then verify the new guard fires: copy the test profile's `samples.csv` fixture, duplicate one
data row (same `LOCUSTAG` and `ASMID`), run
`python3 nextflow/bin/build_species_table.py --samples <dup-copy> --key LOCUSTAG -o /tmp/species.csv.gz`
directly, confirm non-zero exit with an "ERROR: duplicate" message.

- [ ] **Step 4: Commit**

```bash
git add nextflow/modules/BFD/MERGE_SAMPLES/main.nf nextflow/bin/build_species_table.py
git commit -m "MERGE_SAMPLES: atomic publish for species.parquet + duplicate-key guard

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 10: Rewrite `build_BFD_duckDB.sh` — views for species/asm_stats/busco_genome

**Files:**
- Modify: `nextflow/bin/build_BFD_duckDB.sh`
- Modify: `nextflow/modules/BFD/BUILD_DUCKDB/main.nf`

**Interfaces:**
- Consumes: `nextflow/bin/generate_bfd_catalog_sql.py` (Task 4).

- [ ] **Step 1: Pass absolute paths from the Nextflow module (fixes the pre-review draft's broken path resolution)**

In `nextflow/modules/BFD/BUILD_DUCKDB/main.nf`, change:
```groovy
    TABLES=${params.tables} DBDIR=. bash ${projectDir}/bin/build_BFD_duckDB.sh
```
to:
```groovy
    TABLES=${params.tables} DBDIR=. \\
    BIN_DIR=${projectDir}/bin \\
    SCHEMA=${projectDir}/../sql/table_schema.json \\
    bash ${projectDir}/bin/build_BFD_duckDB.sh
```
This is required, not optional: the standalone script's own `${BIN_DIR:-nextflow/bin}` /
`${SCHEMA:-sql/table_schema.json}` fallbacks (Step 2 below) are relative to whatever directory
the script happens to run in — for this module invocation that's a Nextflow work directory, where
neither relative path resolves to anything. Confirmed by Fable during plan review: without this
change, every real (non-manual) catalog rebuild fails.

- [ ] **Step 2: Replace the `species`/`asm_stats`/`busco_genome` `CREATE TABLE` blocks with a codegen call**

In `nextflow/bin/build_BFD_duckDB.sh`, find the `# ─── species ───` and `# ─── asm_stats ───`
blocks and the `busco_genome` block further down (search `CREATE TABLE busco_genome`). Replace
all three with one block placed where `species` currently is:

```bash
# ─── species / asm_stats / busco_genome: schema-contract-driven VIEWs ───────
# These three are VIEWs over tables/*.parquet, not materialized copies -- see
# sql/table_schema.json and
# docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md. Everything
# below this block is still a materialized CREATE TABLE copy (fast-follow
# migrates the rest to the same pattern).
BIN_DIR="${BIN_DIR:?BIN_DIR must be set to an absolute path to nextflow/bin}"
SCHEMA="${SCHEMA:?SCHEMA must be set to an absolute path to sql/table_schema.json}"
# A VIEW's read_parquet('...') path is evaluated at *query* time by whatever
# process opens db/BFD.duckdb (e.g. the MCP server, possibly from a
# different cwd), so it must be absolute -- resolve $SRC defensively even
# though the production caller (BUILD_DUCKDB/main.nf) already passes
# params.tables, which is itself already absolute.
SRC_ABS="$(cd "$SRC" && pwd)"
CATALOG_SQL="$(mktemp)"
python3 "$BIN_DIR/generate_bfd_catalog_sql.py" --schema "$SCHEMA" --tables-dir "$SRC_ABS" > "$CATALOG_SQL"
duckdb "$DB" < "$CATALOG_SQL"
rm -f "$CATALOG_SQL"
```

Note the deliberate `${VAR:?message}` (not `${VAR:-default}`) for `BIN_DIR`/`SCHEMA` — a silent
wrong-default here would be a repeat of the exact bug this step exists to fix; failing loudly if
the caller (Step 1) forgot to set them is the correct behavior.

Remove the old `species`/`asm_stats`/`busco_genome` `CREATE TABLE ... CREATE UNIQUE INDEX ...`
blocks entirely (the `busco_genome` removal also means deleting its
`if [ -f ... busco_genome.parquet ]` guard — the codegen's own existence check, Task 4, now
covers this by skipping the view instead).

- [ ] **Step 3: Run the full build script against real data and verify the views work**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
module load duckdb
BIN_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin \
SCHEMA=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/sql/table_schema.json \
bash /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin/build_BFD_duckDB.sh
duckdb db/BFD.duckdb -c "SELECT COUNT(*) FROM species;"
duckdb db/BFD.duckdb -c "SELECT COUNT(*) FROM asm_stats;"
duckdb db/BFD.duckdb -c "SELECT LOCUSTAG, ASMID, TOTAL_LENGTH, GC_PERCENT FROM asm_stats LIMIT 5;"
```
Expected: both `COUNT(*)` queries return rows without error; `TOTAL_LENGTH`/`GC_PERCENT` (the
aliased, `CAST`-typed columns) are present and populated. If `tables/busco_genome.parquet`
doesn't exist yet in this environment (likely — `run_busco_genome` defaults `false`), confirm the
script still completes successfully and prints a `-- SKIP: busco_genome` line rather than
aborting; `SELECT * FROM busco_genome` should then fail with a clear "Table with name
busco_genome does not exist" error, not a build-time crash. Same check for `_table_manifest` if
`tables/_manifest/` doesn't exist yet in this environment before any Task 7/8/9 module has run.

- [ ] **Step 4: Confirm a downstream consumer still works unmodified**

```bash
python3 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/bin/pick_representative_strain.py --help
```
(Already Parquet-native per the earlier audit — this just confirms the script still imports/runs
after the surrounding changes; no code change expected here.)

- [ ] **Step 5: Commit**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
git add nextflow/bin/build_BFD_duckDB.sh nextflow/modules/BFD/BUILD_DUCKDB/main.nf
git commit -m "build_BFD_duckDB.sh: species/asm_stats/busco_genome become VIEWs

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 11: Warn→fail-loud upgrade for duplicate gene IDs (lower priority)

**Files:**
- Modify: `nextflow/bin/calculate_intergenic.py`
- Modify: `nextflow/bin/build_genestats_table.py`

**Interfaces:** none beyond stdlib `sys.exit`.

- [ ] **Step 1: Locate the existing warning**

Run: `grep -n "WARNING: Duplicate gene ID" nextflow/bin/calculate_intergenic.py nextflow/bin/build_genestats_table.py`

- [ ] **Step 2: Change `print(...)` to `sys.exit(...)` in both files**

For each match, replace:
```python
print(f"WARNING: Duplicate gene ID {gene_id} in {gff}")
```
with:
```python
sys.exit(f"ERROR: Duplicate gene ID {gene_id} in {gff} -- refusing to continue")
```
(add `import sys` at the top of either file if not already imported)

- [ ] **Step 3: Verify manually against a GFF3 fixture with a repeated gene ID**

Construct or reuse an existing test GFF3, duplicate one `ID=gene-` line, run each script against
it, confirm non-zero exit with the new "ERROR: Duplicate gene ID" message instead of a printed
warning and continued execution.

- [ ] **Step 4: Run the existing repo test suite to confirm no regression**

Run: `python3 nextflow/tests/test_hash_bucket_parity.py`
Expected: passes as before (unrelated script, but confirms nothing in `nextflow/bin` broke on
import).

- [ ] **Step 5: Commit**

```bash
git add nextflow/bin/calculate_intergenic.py nextflow/bin/build_genestats_table.py
git commit -m "calculate_intergenic.py, build_genestats_table.py: fail loud on duplicate gene ID

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 12: Remove dead scripts

**Files:**
- Delete: `nextflow/bin/collect_chrom_info.py`
- Delete: `nextflow/bin/rewrite_abinitio_param_paths.py`

- [ ] **Step 1: Re-confirm zero references before deleting**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
grep -rn "collect_chrom_info\|rewrite_abinitio_param_paths" nextflow/ --include="*.nf" --include="*.config" | grep -v /work/
```
Expected: no output (confirms the fork audit's finding still holds).

- [ ] **Step 2: Remove the files**

```bash
git rm nextflow/bin/collect_chrom_info.py nextflow/bin/rewrite_abinitio_param_paths.py
```

- [ ] **Step 3: Commit**

```bash
git commit -m "Remove dead scripts: collect_chrom_info.py, rewrite_abinitio_param_paths.py

Both confirmed unreferenced by any .nf/module/config file, missed by the
earlier asm_reports/collect_asm_stats.py dead-code cleanup pass.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

### Task 13: End-to-end smoke test on the real (non-stub) pipeline

**Files:** none (verification only).

- [ ] **Step 1: Run a small real `--pipeline BFD` slice**

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_annotation_asco
sbatch ../Fungi_BFD/nextflow/run_functional.sh --n_test 3
```
(with `merge_all: true`, `skip_merge: false`, `run_build_duckdb: true` set in a **scratch copy**
of `nextflow/param_files/params_functional.yaml` for this run only — don't edit the checked-in
file, since this plan doesn't change those pipeline-wide defaults)

- [ ] **Step 2: After the job completes, verify the manifest and views reflect the real run**

```bash
cat tables/_manifest/asm_stats.json
cat tables/_manifest/busco_genome.json
cat tables/_manifest/species.json
module load duckdb
duckdb db/BFD.duckdb -c "SELECT * FROM _table_manifest;"
duckdb db/BFD.duckdb -c "SELECT LOCUSTAG, ASMID, complete_pct FROM busco_genome LIMIT 5;"
```
Expected: manifests show recent `built_at` timestamps and `built_by` matching the module names;
the views return real rows for the 3 test genomes.

- [ ] **Step 3: Confirm no leftover temp or lock files anywhere under `tables/`**

```bash
find tables/ -name ".tmp.*" -o -name ".lock.*"
```
Expected: empty output — a leftover `.tmp.*` would indicate a crash mid-publish; a leftover
`.lock.*` file itself is harmless (an empty file `flock` reuses next time) but worth noting if
present.

No commit for this task — it's a verification pass, not a code change.
