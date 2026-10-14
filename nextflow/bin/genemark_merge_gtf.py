#!/usr/bin/env python3
"""Concatenate chunk GTFs in the order given and renumber gene/transcript IDs.

gmes_petap.pl numbers genes 1..N contiguously inside each run, so a running offset per
chunk reproduces the numbering of a monolithic run. Fails (exit 1) if a chunk's IDs are
not contiguous or if the merged gene count differs from the sum over chunks.
Plain or .gz input. Writes the merged GTF to --out.
"""
import argparse
import gzip
import re
import sys

G = re.compile(r'"(\d+)_g"')
T = re.compile(r'"(\d+)_t"')


def open_in(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("gtfs", nargs="+")
    a = ap.parse_args()
    g_off = t_off = 0
    expected_genes = 0
    written_gene_lines = 0
    with open(a.out, "w") as out:
        for p in a.gtfs:
            g_ids, t_ids = set(), set()
            with open_in(p) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    m = G.search(line)
                    if m:
                        g_ids.add(int(m.group(1)))
                        line = G.sub(lambda x: '"%d_g"' % (int(x.group(1)) + g_off), line, count=1)
                    m = T.search(line)
                    if m:
                        t_ids.add(int(m.group(1)))
                        line = T.sub(lambda x: '"%d_t"' % (int(x.group(1)) + t_off), line, count=1)
                    if line.split("\t")[2:3] == ["gene"]:
                        written_gene_lines += 1
                    out.write(line)
            if g_ids and (min(g_ids) != 1 or max(g_ids) != len(g_ids)):
                sys.exit("ERROR: %s: gene IDs are not contiguous 1..N (n=%d, max=%d)" % (p, len(g_ids), max(g_ids)))
            if t_ids and (min(t_ids) != 1 or max(t_ids) != len(t_ids)):
                sys.exit("ERROR: %s: transcript IDs are not contiguous 1..N" % p)
            expected_genes += len(g_ids)
            g_off += len(g_ids)
            t_off += len(t_ids)
            print("%s\tgenes=%d\ttranscripts=%d" % (p, len(g_ids), len(t_ids)), file=sys.stderr)
    if written_gene_lines != expected_genes:
        sys.exit("ERROR: merged gene lines %d != sum of chunk genes %d" % (written_gene_lines, expected_genes))
    print("total genes=%d transcripts=%d" % (g_off, t_off), file=sys.stderr)


main()
