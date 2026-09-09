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
