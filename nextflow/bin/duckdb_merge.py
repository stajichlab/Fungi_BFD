#!/usr/bin/env python3
"""duckdb_merge.py -- run a DuckDB SQL merge template over a list of per-genome
result files and write Parquet.

Usage:
    duckdb_merge.py --sql nextflow/sql/merge/pfam.sql --files list.txt \\
        --out pfam.parquet [--threads 4] [--memory-limit 18GB] [--temp-dir $SCRATCH] \\
        [--param KEY=VALUE ...]

--files: one path per line, or a merge manifest (path in TAB field 1).
The template is plain DuckDB SQL with these placeholders:
    {files}      SQL list literal of the input paths: ['a.gz','b.gz',...]
    {out}        output Parquet path (quoted SQL string)
    {KEY}        any --param KEY=VALUE (inserted verbatim)
It must end in a COPY ... TO {out} (FORMAT PARQUET ...) statement.

Why: the per-tool Python parsers (merge_*.py) wrote one large CSV that DuckDB
then re-read with read_csv_auto(sample_size=-1). On 1,018 genomes the
SwissProt Python step took 1,266 s of 1,507 s; DuckDB parsing the DIAMOND
tables directly took 189 s with the same rows (2026-11-06). The SQL goes to
duckdb on stdin: one shell argument is limited to 128 KB.
"""
import argparse, os, subprocess, sys, time


def sql_str(s):
    return "'" + s.replace("'", "''") + "'"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sql", required=True)
    ap.add_argument("--files", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("NCPU", "4")))
    ap.add_argument("--memory-limit", default=None)
    ap.add_argument("--temp-dir", default=None)
    ap.add_argument("--param", action="append", default=[])
    a = ap.parse_args()

    body = open(a.sql).read()
    paths = [l.rstrip("\n").split("\t")[0] for l in open(a.files)]
    paths = [p for p in paths if p]
    if not paths and "{files}" in body:
        sys.exit(f"ERROR: no input files in {a.files}")
    tmp_out = a.out + ".tmp"
    subs = {"files": "[" + ",".join(sql_str(p) for p in paths) + "]", "out": sql_str(tmp_out)}
    for kv in a.param:
        k, v = kv.split("=", 1)
        subs[k] = v
    for k, v in subs.items():
        body = body.replace("{" + k + "}", v)
    # insertion order is kept: inputs arrive sorted by genome, which keeps the
    # Parquet row groups prunable by protein_id / species_prefix.
    head = [f"SET threads = {a.threads};"]
    if a.memory_limit:
        head.append(f"SET memory_limit = {sql_str(a.memory_limit)};")
    if a.temp_dir:
        head.append(f"SET temp_directory = {sql_str(a.temp_dir)};")
    sql = "\n".join(head) + "\n" + body
    t0 = time.time()
    res = subprocess.run(["duckdb"], input=sql, capture_output=True, text=True)
    if res.returncode != 0 or "Error" in res.stderr:
        sys.exit(f"ERROR: duckdb failed:\n{res.stderr.strip()[:3000]}")
    os.replace(tmp_out, a.out)
    n = subprocess.run(["duckdb", "-csv", "-noheader", "-c",
                        f"SELECT count(*) FROM read_parquet({sql_str(a.out)})"],
                       capture_output=True, text=True, check=True).stdout.strip()
    print(f"{a.out}: {n} rows from {len(paths)} files in {time.time() - t0:.0f} s ({os.path.basename(a.sql)})")


if __name__ == "__main__":
    main()
