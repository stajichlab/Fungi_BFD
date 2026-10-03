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
