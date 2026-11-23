#!/usr/bin/env python3
"""Merge per-genome antiSMASH-lite tables into four zstd Parquet tables (schema v2).

Input: any mix of <LOCUSTAG>.antismash.{regions,cds,known,known_pairs}.tsv.gz written by
modules/BFD/ANTISMASH (via bin/extract_antismash_subset.py --minimal). The `genome` column of
every file is the LOCUSTAG, so it is used as `locustag` directly.
Output (in --outdir): bgc_cluster, bgc_cluster_gene, bgc_mibig_hit, bgc_mibig_pair (.parquet).
Join key between the tables: (locustag, record, region). gene_id joins to gene_info.gene_id.
Same columns, types and sort order as DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/
build_bgc_parquet.sql (schema v2, 2026-11-22). No sequences are stored.

Header-only files (genomes with no cluster) are fine. Files with no header at all (empty
gzip, e.g. from the stub block) are skipped. Hits above --topn are dropped (default 5).

Usage: merge_antismash.py [-o OUTDIR] [--topn 5] FILE [FILE ...]
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
    ap.add_argument("--threads", type=int, default=int(os.environ.get("NCPU", "4")))
    a = ap.parse_args()

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
    L = {k: "[" + ",".join(sql_str(p) for p in v) + "]" for k, v in by.items()}
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda n: sql_str(os.path.join(a.outdir, n + ".parquet.tmp"))
    n = a.topn
    cp = "(FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000)"
    sql = f"""
SET threads = {a.threads};
CREATE TEMP TABLE r AS SELECT * FROM read_csv({L['regions']}, delim='\\t', header=true, auto_detect=true, types={{'contig_edge':'BOOLEAN'}});
CREATE TEMP TABLE c AS SELECT * FROM read_csv({L['cds']}, delim='\\t', header=true, auto_detect=true, types={{'locus_tag':'VARCHAR','protein_id':'VARCHAR'}});
CREATE TEMP TABLE k AS SELECT * FROM read_csv({L['known']}, delim='\\t', header=true, auto_detect=true);
CREATE TEMP TABLE p AS SELECT * FROM read_csv({L['known_pairs']}, delim='\\t', header=true, auto_detect=true, types={{'query_locus_tag':'VARCHAR','mibig_protein':'VARCHAR'}});
CREATE TEMP TABLE best AS
  SELECT genome, record, region, arg_max(mibig_accession, (similarity, -rank)) AS best_mibig, max(similarity) AS best_similarity,
         arg_max(rank, (similarity, -rank)) AS best_rank
  FROM k WHERE rank <= {n} GROUP BY genome, record, region;
COPY (
  SELECT r.genome AS locustag, r.record, r.region::USMALLINT AS region, r.start::UINTEGER AS start, r."end"::UINTEGER AS "end", r.length::UINTEGER AS length,
         r.products, r.categories, r.contig_edge, r.core_start::UINTEGER AS core_start, r.core_end::UINTEGER AS core_end,
         r.candidate_kinds, r.n_cds::USMALLINT AS n_cds, r.n_biosynthetic::USMALLINT AS n_biosynthetic, r.n_additional::USMALLINT AS n_additional,
         b.best_mibig, b.best_similarity::UTINYINT AS best_similarity, b.best_rank::UTINYINT AS best_rank
  FROM r LEFT JOIN best b ON b.genome=r.genome AND b.record=r.record AND b.region=r.region
  ORDER BY r.genome, r.record, r.region) TO {out('bgc_cluster')} {cp};
COPY (
  SELECT c.genome AS locustag, c.record, c.region::USMALLINT AS region, c.locus_tag AS gene_id, c.start::UINTEGER AS start, c."end"::UINTEGER AS "end",
         (c.strand='-') AS is_minus, c.gene_kind, c.functions
  FROM c ORDER BY c.genome, c.record, c.region, c.start) TO {out('bgc_cluster_gene')} {cp};
COPY (
  SELECT k.genome AS locustag, k.record, k.region::USMALLINT AS region, k.rank::UTINYINT AS rank, k.mibig_accession, k.hit_description AS compounds, k.hit_cluster_type,
         k.similarity::UTINYINT AS similarity, k.core_gene_hits::USMALLINT AS core_gene_hits, k.gene_hits::USMALLINT AS gene_hits, k.n_entry_proteins::USMALLINT AS n_entry_proteins,
         k.synteny_score::SMALLINT AS synteny_score, k.blast_score::FLOAT AS blast_score
  FROM k WHERE k.rank <= {n} ORDER BY k.genome, k.record, k.region, k.rank) TO {out('bgc_mibig_hit')} {cp};
COPY (
  SELECT p.genome AS locustag, p.record, p.region::USMALLINT AS region, p.rank::UTINYINT AS rank, p.mibig_accession, p.query_locus_tag AS gene_id,
         p.mibig_gene, p.mibig_protein, p.perc_ident::UTINYINT AS perc_ident, p.perc_coverage::FLOAT AS perc_coverage, p.blast_score::FLOAT AS blast_score, p.evalue::DOUBLE AS evalue
  FROM p WHERE p.rank <= {n} ORDER BY p.genome, p.record, p.region, p.rank, p.query_locus_tag) TO {out('bgc_mibig_pair')} {cp};
"""
    res = subprocess.run(["duckdb"], input=sql, capture_output=True, text=True)
    if res.returncode != 0 or "Error" in res.stderr:
        sys.exit("ERROR: duckdb failed:\n" + res.stderr.strip()[:3000])
    for t in ("bgc_cluster", "bgc_cluster_gene", "bgc_mibig_hit", "bgc_mibig_pair"):
        final = os.path.join(a.outdir, t + ".parquet")
        os.replace(final + ".tmp", final)
        cnt = subprocess.run(["duckdb", "-csv", "-noheader", "-c", f"SELECT count(*) FROM read_parquet({sql_str(final)})"],
                             capture_output=True, text=True, check=True).stdout.strip()
        print(f"{final}: {cnt} rows", file=sys.stderr)
    print(f"genomes: {len(by['regions'])} regions files; skipped {skipped} empty files", file=sys.stderr)


if __name__ == "__main__":
    main()
