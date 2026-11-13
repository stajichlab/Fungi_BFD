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
STR_COLS = ["pair", "family", "g1", "g2", "gene1", "gene2", "node", "status"]
# alignmentcoverage/alignmentidentity/strippedalignmentlength are kept: dS is only
# reliable when enough codons survive wgd's family-wide gap stripping
# (strippedalignmentlength is in nucleotides; codons = value / 3).
FLOAT_COLS = ["N", "S", "dN", "dN/dS", "dS", "alignmentlength",
              "strippedalignmentlength", "alignmentcoverage", "alignmentidentity", "t",
              "pair_codons", "ng_dS", "ng_codons", "copies1", "copies2", "npairs"]
# copies1/copies2/npairs (collapsed identical sequences) come from wgd_pair_ks.py: npairs = number of gene
# pairs the row represents (copies1 x copies2; n(n-1)/2 for status 'identical').
# pair_codons = (N+S)/3 is added by this script. It is the number of codons codeml
# used for the pair. Under wgd ksd --pairwise the strippedalignmentlength and
# alignmentlength columns are family-level constants, so this column is the only
# per-pair alignment-length measure.


NODE_INT = ["node_pairs", "node_rep_pairs", "node_pass", "node_sampled", "node_gated", "n_valid"]
NODE_FLOAT = ["dS_mean", "dS_min", "dS_max", "identity_mean"]


def node_table(path, df, genome, prefix):
    """One row per (family, node): counts from <path>.nodes.tsv.gz plus dS statistics of the kept sampled pairs.
    dS statistics cover the rows that passed the merge filters (--max-ds, --min-pair-codons, ...).
    Identical-copy nodes (kind = identical) have dS = 0."""
    npath = path + ".nodes.tsv.gz"
    if not os.path.exists(npath):
        return None
    nd = pd.read_csv(npath, sep="\t", dtype={"family": str, "node": str, "kind": str})
    if "node_gated" not in nd.columns:
        nd["node_gated"] = 0
    if len(df):
        d = df[["family", "node", "dS", "alignmentidentity"]].copy()
        d = d[d["dS"].notna()]
        g = d.groupby(["family", "node"]).agg(
            n_valid=("dS", "size"), dS_mean=("dS", "mean"), dS_min=("dS", "min"),
            dS_max=("dS", "max"), identity_mean=("alignmentidentity", "mean")).reset_index()
        nd = nd.merge(g, on=["family", "node"], how="left")
    else:
        for c in ["n_valid"] + NODE_FLOAT:
            nd[c] = float("nan")
    nd["n_valid"] = nd["n_valid"].fillna(0)
    nd.insert(0, "species_prefix", prefix)
    nd.insert(0, "genome", genome)
    return nd


def identical_stats(path):
    """Collapsed identical-CDS groups of a genome, from <path>.nodes.tsv.gz (kind = identical).
    Returns (groups, genes in groups, largest group, identical gene pairs / all gene pairs in the families).
    A group of n copies stands for n(n-1)/2 gene pairs, so n is recovered from the pair count."""
    npath = path + ".nodes.tsv.gz"
    if not os.path.exists(npath):
        return 0, 0, 0, float("nan")
    nd = pd.read_csv(npath, sep="\t", usecols=["kind", "node_pairs"])
    tot = float(nd["node_pairs"].sum())
    idn = nd[nd["kind"] == "identical"]["node_pairs"].astype("float64")
    if idn.empty:
        return 0, 0, 0, 0.0 if tot else float("nan")
    copies = ((1.0 + (1.0 + 8.0 * idn) ** 0.5) / 2.0).round().astype("int64")
    return int(len(idn)), int(copies.sum()), int(copies.max()), float(idn.sum() / tot) if tot else float("nan")


def species_prefix_of(gene):
    s = str(gene).split("_", 1)
    return s[0] if s else ""


def genome_of(path):
    # <sampletag>.cds-transcripts.fa.tsv.ks.tsv -> sampletag (== makeSampleTag
    # SPECIES_STRAIN, the BFD primary key; species_prefix/LOCUSTAG is derived
    # separately from the wgd gene ids).
    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    stem = base[: -len(KS_SUFFIX)]
    return stem.removesuffix(".cds-transcripts.fa")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", nargs="+", required=True,
                    help="ks.tsv glob pattern(s); matched on disk")
    ap.add_argument("-o", "--out", default="wgd.ks.parquet",
                    help="output merged pairs parquet")
    ap.add_argument("-s", "--summary", default="wgd.ks.summary.parquet",
                    help="output per-genome summary parquet")
    ap.add_argument("--max-ds", type=float, default=None,
                    help="drop pairs with dS above this value (saturated; wgd "
                         "plots 0-5). Default: keep all")
    ap.add_argument("--min-codons", type=int, default=None,
                    help="drop pairs with fewer than this many codons left after "
                         "wgd's gap stripping (strippedalignmentlength/3). Pairs "
                         "with a missing value are dropped. Default: keep all")
    ap.add_argument("--min-pair-codons", type=int, default=None,
                    help="drop pairs with (N+S)/3 below this value (per-pair codons "
                         "used by codeml; use with wgd ksd --pairwise). Pairs with a "
                         "missing value are dropped. Default: keep all")
    ap.add_argument("--failed-glob", default=None,
                    help="glob of *.failed.txt files written by WGD_PAIR_KS; the merge is written, "
                         "then the script exits 3 when any exist")
    ap.add_argument("--nodes-out", default=None,
                    help="write the per-node table (one row per genome, family, node) to this parquet; "
                         "needs the <ks table>.nodes.tsv.gz files from wgd_pair_ks.py --node-sample")
    ap.add_argument("--check-out", default=None,
                    help="write the problem list (or OK) to this file and exit 0 even when problems "
                         "exist, so a Nextflow process still publishes the tables; the workflow "
                         "logs an error when the file is not OK")
    ap.add_argument("--expected-genomes", type=int, default=0,
                    help="genomes selected in this run; the script exits 3 (after writing) when "
                         "fewer ks tables exist")
    ap.add_argument("--drop-null-ds", action="store_true",
                    help="drop pairs with no valid dS estimate")
    args = ap.parse_args()

    files = []
    for pat in args.glob:
        files.extend(glob.glob(pat))
    # side files of wgd_pair_ks.py live in the same folder and match *.ks.tsv*
    files = sorted(f for f in files if not f.endswith((".nodes.tsv.gz", ".groups.tsv.gz", ".status.txt", ".failed.txt")))
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

    summaries = []  # (genome, species_prefix, n_pairs_raw, n_pairs, n_pairs_with_ds, n_families, n_nodes, n_nodes_with_ds, hybrid, identical stats x4)
    n_rows = 0
    node_schema = pa.schema(
        [("genome", string_type), ("species_prefix", string_type), ("family", string_type), ("node", string_type),
         ("kind", string_type)] + [(c, pa.int64()) for c in NODE_INT[:5]] + [("n_valid", pa.int64())] +
        [(c, float_type) for c in NODE_FLOAT])
    nodes_writer = (pq.ParquetWriter(args.nodes_out, node_schema, compression="zstd", use_dictionary=True,
                                     version="2.6", data_page_version="2.0") if args.nodes_out else None)
    n_nodes_total = 0
    no_nodes_file = 0
    with pq.ParquetWriter(args.out, schema, compression="zstd",
                          use_dictionary=True, version="2.6", data_page_version="2.0") as w:
        for path in files:
            genome = genome_of(path)
            try:
                df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                                 na_values=["NaN", "nan", ""])
            except pd.errors.EmptyDataError:
                df = pd.DataFrame(columns=STR_COLS + FLOAT_COLS)
            n_raw = len(df)
            if not df.empty and "N" in df.columns and "S" in df.columns:
                df["pair_codons"] = (pd.to_numeric(df["N"], errors="coerce")
                                     + pd.to_numeric(df["S"], errors="coerce")) / 3.0
                if "status" in df.columns and "strippedalignmentlength" in df.columns:
                    # identical-copy rows (collapsed groups) have no N/S; codons = CDS length / 3
                    ident = df["status"] == "identical"
                    df.loc[ident, "pair_codons"] = pd.to_numeric(
                        df.loc[ident, "strippedalignmentlength"], errors="coerce") / 3.0
            prefix0 = species_prefix_of(df["gene1"].iloc[0]) if n_raw and "gene1" in df.columns else ""
            if not df.empty and (args.max_ds is not None or args.drop_null_ds
                                 or args.min_codons is not None
                                 or args.min_pair_codons is not None):
                ds = pd.to_numeric(df["dS"], errors="coerce") if "dS" in df.columns \
                    else pd.Series(float("nan"), index=df.index)
                keep = pd.Series(True, index=df.index)
                if args.drop_null_ds:
                    keep &= ds.notna()
                if args.max_ds is not None:
                    keep &= ~(ds > args.max_ds)  # NaN is kept unless --drop-null-ds
                if args.min_codons is not None:
                    sl = pd.to_numeric(df["strippedalignmentlength"], errors="coerce") \
                        if "strippedalignmentlength" in df.columns \
                        else pd.Series(float("nan"), index=df.index)
                    keep &= (sl / 3) >= args.min_codons  # NaN -> False
                if args.min_pair_codons is not None:
                    pc = df["pair_codons"] if "pair_codons" in df.columns \
                        else pd.Series(float("nan"), index=df.index)
                    keep &= pc >= args.min_pair_codons  # NaN -> False
                df = df[keep].copy()
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
            n_nodes = n_nodes_ds = 0
            if nodes_writer is not None:
                nt = node_table(path, df, genome, prefix0)
                if nt is None:
                    no_nodes_file += 1
                else:
                    n_nodes = len(nt)
                    n_nodes_ds = int((nt["n_valid"] > 0).sum())
                    n_nodes_total += n_nodes
                    for c in NODE_INT[:5] + ["n_valid"]:
                        nt[c] = pd.to_numeric(nt[c], errors="coerce").fillna(0).astype("int64")
                    for c in NODE_FLOAT:
                        nt[c] = pd.to_numeric(nt[c], errors="coerce")
                    nodes_writer.write_table(pa.Table.from_pandas(nt[node_schema.names], schema=node_schema,
                                                                  preserve_index=False))
            if n_raw:
                summaries.append((
                    genome,
                    prefix0,
                    n_raw,
                    len(df),
                    int(df["dS"].notna().sum()) if len(df) else 0,
                    int(df["family"].nunique()) if len(df) else 0,
                    n_nodes,
                    n_nodes_ds,
                    int(" x " in genome.replace("_", " ")),
                ) + identical_stats(path))
            print(f"  {genome}: {n_raw:,} raw -> {len(df):,} rows", file=sys.stderr)

    # per-genome summary
    sum_schema = pa.schema([
        ("genome", string_type),
        ("species_prefix", string_type),
        ("n_pairs_raw", pa.int64()),   # pairs in ks.tsv before filtering
        ("n_pairs", pa.int64()),       # pairs kept in wgd.ks.parquet
        ("n_pairs_with_ds", pa.int64()),
        ("n_families", pa.int64()),
        ("n_nodes", pa.int64()),          # nodes of the family trees (0 when there is no nodes file)
        ("n_nodes_with_ds", pa.int64()),  # nodes with at least one kept dS
        ("hybrid", pa.int64()),           # 1 when the genome name has " x " (interspecific hybrid)
        ("n_identical_groups", pa.int64()),            # groups of genes with an identical CDS inside a family
        ("n_genes_in_identical_groups", pa.int64()),
        ("largest_identical_group", pa.int64()),
        ("identical_pair_share", pa.float64()),        # identical gene pairs / all gene pairs of the families
    ])
    if summaries:
        sum_df = pd.DataFrame(summaries, columns=[f.name for f in sum_schema])
        sum_tbl = pa.Table.from_pandas(sum_df, schema=sum_schema, preserve_index=False)
    else:
        sum_tbl = sum_schema.empty_table()
    pq.write_table(sum_tbl, args.summary, compression="zstd",
                   use_dictionary=True, version="2.6")
    if nodes_writer is not None:
        nodes_writer.close()
        print(f"Nodes: {args.nodes_out} ({n_nodes_total:,} nodes; {no_nodes_file} genomes without a nodes file)",
              file=sys.stderr)

    print(f"Done: {args.out} ({n_rows:,} merged rows), {args.summary} "
          f"({sum_tbl.num_rows} genomes)", file=sys.stderr)
    rc = 0
    problems = []
    if args.failed_glob:
        failed = sorted(glob.glob(args.failed_glob))
        have = {genome_of(x) for x in files}
        if failed:
            print(f"Checking {len(failed)} failure list(s) from WGD_PAIR_KS:", file=sys.stderr)
            for f in failed:
                for line in open(f):
                    gid = line.split(":", 1)[0].strip()
                    if gid in have:
                        continue  # a later run produced the table; the marker is stale
                    print(f"  {os.path.basename(f)}: {line.rstrip()}", file=sys.stderr)
                    problems.append(f"failed: {os.path.basename(f)}: {line.rstrip()}")
                    rc = 3
    if args.expected_genomes and len(files) < args.expected_genomes:
        rc = 3
        print(f"ERROR: {len(files)} ks tables found, {args.expected_genomes} genomes expected; "
              f"genomes are missing (dropped tasks, or no CDS file)", file=sys.stderr)
        problems.append(f"missing genomes: {len(files)} tables, {args.expected_genomes} expected")
    if args.nodes_out and no_nodes_file:
        rc = 3
        print(f"ERROR: {no_nodes_file} genomes have no nodes file", file=sys.stderr)
        problems.append(f"genomes without a nodes file: {no_nodes_file}")
    if args.check_out:
        with open(args.check_out, "w") as fh:
            fh.write("\n".join(problems) + "\n" if problems else "OK\n")
        return 0
    return rc


if __name__ == "__main__":
    sys.exit(main())
