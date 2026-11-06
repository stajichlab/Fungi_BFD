#!/usr/bin/env python3
"""manifest_to_parquet.py -- build one Parquet table from many per-genome CSV
files (plain or .gz) listed in a merge manifest, with DuckDB reading the files
directly.

Usage:
    manifest_to_parquet.py --manifest gene_stats.manifest.txt \\
        --suffix .gene_exons.csv.gz --out gene_exons.parquet \\
        [--types "exon_id:VARCHAR,order:BIGINT,..."] [--delim ,] \\
        [--threads 4] [--memory-limit 18GB] [--temp-dir $SCRATCH] [--sample-files 300]

Manifest: one line per file, TAB-separated, path in field 1 (toManifest() in
modules/common/utils.nf). --suffix keeps only paths ending in it.

Column types: detected by DuckDB on an evenly spaced sample of --sample-files
files (full scan of each sampled file), then fixed for the full read; --types
overrides individual columns. A value that does not fit the fixed type in a
file outside the sample fails the build loudly (no silent coercion).

Why (profiled 2026-11-06, 1,020 genomes of gene_exons, 4 CPUs): the former
zcat-concatenate-then-read_csv_auto(sample_size=-1) path took 116 s and 3.4 GB
RSS and, at ~23k genomes, failed on a single .csv.gz with "Trying to seek a
compressed CSV File". Reading the file list directly with fixed types took
26 s and 0.86 GB, with identical rows (count and checksum). Writing the
Parquet to $SCRATCH first was no faster.

The file list goes into a SQL file passed on stdin: one shell argument is
limited to 128 KB.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time


def sql_str(s):
    return "'" + s.replace("'", "''") + "'"


def file_list_sql(paths):
    return "[" + ",".join(sql_str(p) for p in paths) + "]"


def run_duckdb(sql, csv_out=False):
    with tempfile.NamedTemporaryFile("w", suffix=".sql", delete=False) as fh:
        fh.write(sql)
        path = fh.name
    try:
        cmd = ["duckdb"] + (["-csv", "-noheader"] if csv_out else [])
        with open(path) as stdin:
            res = subprocess.run(cmd, stdin=stdin, capture_output=True, text=True)
        if res.returncode != 0 or "Error" in res.stderr:
            sys.exit(f"ERROR: duckdb failed:\n{res.stderr.strip()[:2000]}")
        return res.stdout
    finally:
        os.unlink(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--types", default="", help="col:TYPE,... overrides")
    ap.add_argument("--delim", default=",")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("NCPU", "4")))
    ap.add_argument("--memory-limit", default=None, help="DuckDB memory_limit, e.g. 18GB")
    ap.add_argument("--temp-dir", default=None, help="DuckDB spill dir (node-local scratch)")
    ap.add_argument("--sample-files", type=int, default=300)
    a = ap.parse_args()

    paths = []
    for line in open(a.manifest):
        p = line.rstrip("\n").split("\t")[0]
        if p and p.endswith(a.suffix):
            paths.append(p)
    if not paths:
        sys.exit(f"ERROR: no files in {a.manifest} end with '{a.suffix}'")
    step = max(1, len(paths) // a.sample_files)
    sample = paths[::step][: a.sample_files]
    delim = sql_str(a.delim)

    settings = [f"SET threads = {a.threads};"]
    if a.memory_limit:
        settings.append(f"SET memory_limit = {sql_str(a.memory_limit)};")
    if a.temp_dir:
        settings.append(f"SET temp_directory = {sql_str(a.temp_dir)};")
    settings = "\n".join(settings) + "\n"

    t0 = time.time()
    desc = run_duckdb(settings + f"SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM "
                      f"read_csv({file_list_sql(sample)}, header=true, delim={delim}, "
                      f"union_by_name=true, sample_size=-1));\n", csv_out=True)
    cols = [tuple(l.split(",", 1)) for l in desc.strip().splitlines()]
    override = dict(kv.split(":", 1) for kv in a.types.split(",") if kv.strip())
    unknown = set(override) - {c for c, _ in cols}
    if unknown:
        sys.exit(f"ERROR: --types names columns not in the files: {sorted(unknown)}")
    cols = [(c, override.get(c, t)) for c, t in cols]
    t_detect = time.time() - t0

    colspec = "{" + ", ".join(f"{sql_str(c)}: {sql_str(t)}" for c, t in cols) + "}"
    tmp_out = a.out + ".tmp"
    t1 = time.time()
    run_duckdb(settings + f"COPY (SELECT * FROM read_csv({file_list_sql(paths)}, header=true, "
               f"delim={delim}, columns={colspec})) TO {sql_str(tmp_out)} (FORMAT PARQUET);\n")
    os.replace(tmp_out, a.out)
    n = run_duckdb(f"SELECT count(*) FROM read_parquet({sql_str(a.out)});\n", csv_out=True).strip()
    print(f"{a.out}: {n} rows from {len(paths)} files; types from {len(sample)} files "
          f"({t_detect:.0f} s), build {time.time() - t1:.0f} s; "
          + ", ".join(f"{c}:{t}" for c, t in cols))


if __name__ == "__main__":
    main()
