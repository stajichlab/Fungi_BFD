#!/usr/bin/env python3
"""Merge per-species wgd ksd ks.tsv files into tables-loadable Parquet.

Follows the BFD MERGE_* convention: one unpartitioned Parquet per table type
in tables/ (here wgd.ks.parquet + wgd.ks.summary.parquet), keyed by
species_prefix (LOCUSTAG, first token of the wgd gene id) and genome
(sampletag, the ks.tsv filename stem).

The ks.tsv files are globbed off disk (from params.outdir/wgd_ksd) rather
than taken from a live Nextflow channel, so a -resume run that only
recomputes part of the dataset still merges *everything* published so far --
the same stale-storeDir-wave rationale as gatedGlobIn in modules/common/utils.nf.

Raw wgd columns are pruned to the analysis-relevant subset; numeric columns
(dN/dS etc.) are wgd "NaN" strings and become NULL.
"""

import argparse
import glob
import os
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

KS_SUFFIX = ".tsv.ks.tsv"

# Raw ks.tsv columns kept in the merged pairs table (plus genome/species_prefix).
# N, S, dN, dS, dN/dS, alignmentlength, t are floats (may be NULL when the
# pair had no valid alignment); pair/family/g1/g2/gene1/gene2 are strings.
# node is the wgd tree-node id the pair was dated at (numeric in the raw file,
# kept as string); wgd mix/peak node-average pairs by (family, node), so the
# peak-summary framework (build_wgd_ksd_summary.py) needs it.
STR_COLS = ["pair", "family", "g1", "g2", "gene1", "gene2", "node"]
FLOAT_COLS = ["N", "S", "dN", "dN/dS", "dS", "alignmentlength", "t"]


def species_prefix_of(gene):
    s = str(gene).split("_", 1)
    return s[0] if s else ""


def genome_of(path):
    # <sampletag>.cds-transcripts.fa.tsv.ks.tsv -> sampletag (== makeSampleTag
    # SPECIES_STRAIN, the BFD primary key; species_prefix/LOCUSTAG is derived
    # separately from the wgd gene ids).
    stem = os.path.basename(path)[: -len(KS_SUFFIX)]
    return stem.removesuffix(".cds-transcripts.fa")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", nargs="+", required=True,
                    help="ks.tsv glob pattern(s); matched on disk")
    ap.add_argument("-o", "--out", default="wgd.ks.parquet",
                    help="output merged pairs parquet")
    ap.add_argument("-s", "--summary", default="wgd.ks.summary.parquet",
                    help="output per-genome summary parquet")
    args = ap.parse_args()

    files = []
    for pat in args.glob:
        files.extend(glob.glob(pat))
    files.sort()
    if not files:
        print(f"ERROR: no *.tsv.ks.tsv files matched {args.glob}; "
              f"refusing to write an empty merged table", file=sys.stderr)
        return 1
    print(f"Merging {len(files)} ks.tsv files", file=sys.stderr)

    string_type = pa.string()
    float_type = pa.float64()
    schema = pa.schema(
        [("genome", string_type)] + [("species_prefix", string_type)] +
        [(c, string_type) for c in STR_COLS] +
        [(c, float_type) for c in FLOAT_COLS]
    )

    summaries = []  # (genome, species_prefix, n_pairs, n_pairs_with_ds, n_families)
    n_rows = 0
    with pq.ParquetWriter(args.out, schema, compression="zstd",
                          use_dictionary=True, version="2.6", data_page_version="2.0") as w:
        for path in files:
            genome = genome_of(path)
            try:
                df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                                 na_values=["NaN", "nan", ""])
            except pd.errors.EmptyDataError:
                df = pd.DataFrame(columns=STR_COLS + FLOAT_COLS)
            if df.empty:
                # keep schema stable; write an empty row group
                table = schema.empty_table()
            else:
                for c in STR_COLS:
                    if c not in df.columns:
                        df[c] = ""
                for c in FLOAT_COLS:
                    if c not in df.columns:
                        df[c] = float("nan")
                    df[c] = pd.to_numeric(df[c], errors="coerce")
                df["genome"] = genome
                # species_prefix: LOCUSTAG from the first token of gene1
                df["species_prefix"] = df["gene1"].map(species_prefix_of)
                table = pa.Table.from_pandas(df[schema.names], schema=schema,
                                             preserve_index=False)
            w.write_table(table)
            n_rows += table.num_rows
            if not df.empty:
                summaries.append((
                    genome,
                    str(df["species_prefix"].iloc[0]) if len(df) else "",
                    len(df),
                    int(df["dS"].notna().sum()),
                    int(df["family"].nunique()),
                ))
            print(f"  {genome}: {len(df):,} rows", file=sys.stderr)

    # per-genome summary
    sum_schema = pa.schema([
        ("genome", string_type),
        ("species_prefix", string_type),
        ("n_pairs", pa.int64()),
        ("n_pairs_with_ds", pa.int64()),
        ("n_families", pa.int64()),
    ])
    if summaries:
        sum_df = pd.DataFrame(summaries, columns=[f.name for f in sum_schema])
        sum_tbl = pa.Table.from_pandas(sum_df, schema=sum_schema, preserve_index=False)
    else:
        sum_tbl = sum_schema.empty_table()
    pq.write_table(sum_tbl, args.summary, compression="zstd",
                   use_dictionary=True, version="2.6")

    print(f"Done: {args.out} ({n_rows:,} merged rows), {args.summary} "
          f"({sum_tbl.num_rows} genomes)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
