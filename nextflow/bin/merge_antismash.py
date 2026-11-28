#!/usr/bin/env python3
"""Merge per-genome antiSMASH-lite tables into four zstd Parquet tables (schema v2).

Input: any mix of <LOCUSTAG>.antismash.{regions,cds,known,known_pairs}.tsv.gz written by
modules/BFD/ANTISMASH (via bin/extract_antismash_subset.py --minimal). The `genome` column of
every file is the LOCUSTAG, so it is used as `locustag` directly.
Output (in --outdir): bgc_cluster, bgc_cluster_gene, bgc_mibig_hit, bgc_mibig_pair (.parquet).
Join key between the tables: (locustag, record, region). gene_id joins to gene_info.gene_id.
The tables are defined in nextflow/sql/merge/antismash.sql (schema v2, 2026-11-22); this script
only collects the files and fills the placeholders. No sequences are stored.

Header-only files (genomes with no cluster) are fine. Files with no header at all (empty
gzip, e.g. from the stub block) are skipped. Hits above --topn are dropped (default 5).

Usage: merge_antismash.py [-o OUTDIR] [--topn 5] FILE [FILE ...]
An argument "@LIST" reads one path per line from the file LIST (needed above ~20,000 genomes,
where the paths no longer fit in the 2 MB argument limit).
"""
import argparse, gzip, os, subprocess, sys

KINDS = ("regions", "cds", "known", "known_pairs")


def sql_str(s):
    return "'" + s.replace("'", "''") + "'"


def usable(path):
    try:
        with gzip.open(path, "rb") as fh:
            return fh.read(1) != b""
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--topn", type=int, default=5)
    ap.add_argument("--sql", default=None, help="SQL template (default: nextflow/sql/merge/antismash.sql)")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("NCPU", "4")))
    a = ap.parse_args()
    files = []
    for f in a.files:
        if f.startswith("@"):
            files += [l.strip() for l in open(f[1:]) if l.strip()]
        else:
            files.append(f)
    a.files = files

    by = {k: [] for k in KINDS}
    skipped = 0
    for f in a.files:
        for k in KINDS:
            if f.endswith(f".antismash.{k}.tsv.gz"):
                if usable(f):
                    by[k].append(f)
                else:
                    skipped += 1
                break
        else:
            sys.exit(f"ERROR: not an antiSMASH table: {f}")
    if not by["regions"]:
        sys.exit("ERROR: no usable regions files")
    for k in KINDS:   # file-name order = LOCUSTAG order
        by[k].sort(key=lambda p: (os.path.basename(p), p))
    sql_path = a.sql or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sql", "merge", "antismash.sql")
    body = open(sql_path).read()
    os.makedirs(a.outdir, exist_ok=True)
    tables = ("bgc_cluster", "bgc_cluster_gene", "bgc_mibig_hit", "bgc_mibig_pair")
    subs = {"regions": by["regions"], "cds": by["cds"], "known": by["known"], "pairs": by["known_pairs"]}
    subs = {k: "[" + ",".join(sql_str(p) for p in v) + "]" for k, v in subs.items()}
    subs["topn"] = str(a.topn)
    for t in tables:
        subs["out_" + t] = sql_str(os.path.join(a.outdir, t + ".parquet.tmp"))
    for k, v in subs.items():
        body = body.replace("{" + k + "}", v)
    sql = f"SET threads = {a.threads};\n" + body
    res = subprocess.run(["duckdb", "-bail"], input=sql, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit("ERROR: duckdb failed:\n" + res.stderr.strip()[:3000])
    for t in tables:
        final = os.path.join(a.outdir, t + ".parquet")
        os.replace(final + ".tmp", final)
        cnt = subprocess.run(["duckdb", "-csv", "-noheader", "-c", f"SELECT count(*) FROM read_parquet({sql_str(final)})"],
                             capture_output=True, text=True, check=True).stdout.strip()
        print(f"{final}: {cnt} rows", file=sys.stderr)
    print(f"genomes: {len(by['regions'])} regions files; skipped {skipped} empty files", file=sys.stderr)


if __name__ == "__main__":
    main()
