#!/usr/bin/env python3
"""Split a FASTA into chunks of whole contigs for scatter/gather GeneMark prediction.

Contigs stay in genome order, so chunk GTFs concatenated in chunk order match the
monolithic gmes_petap.pl order. Chunks are sized by contig count (measured run time
follows contig count, not bases). A chunk smaller than --min-chunk-bp is merged into the
previous chunk, because gmes_petap.pl exits with an error on input under 100,000 bp.
Output: <outdir>/chunk_NN.fa.gz (gzip). Prints one line per chunk to stdout.
"""
import argparse
import gzip
import math
import os
import subprocess
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        return subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True).stdout
    return open(path)


def read_fasta(path):
    name, seq_len, lines = None, 0, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, seq_len, lines
                name, seq_len, lines = line, 0, []
            else:
                lines.append(line)
                seq_len += len(line.strip())
    if name is not None:
        yield name, seq_len, lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--contigs-per-chunk", type=int, required=True)
    ap.add_argument("--min-chunk-bp", type=int, default=100000)
    a = ap.parse_args()

    recs = list(read_fasta(a.input))
    n = len(recs)
    if n == 0:
        sys.exit("ERROR: no sequences in " + a.input)
    k = max(1, math.ceil(n / a.contigs_per_chunk))
    base, extra = divmod(n, k)
    sizes = [base + (1 if i < extra else 0) for i in range(k)]
    bounds, start = [], 0
    for s in sizes:
        bounds.append([start, start + s])
        start += s

    def chunk_bp(b):
        return sum(r[1] for r in recs[b[0]:b[1]])

    # fold undersized chunks into their predecessor (or successor for the first)
    i = 0
    while len(bounds) > 1 and i < len(bounds):
        if chunk_bp(bounds[i]) < a.min_chunk_bp:
            if i > 0:
                bounds[i - 1][1] = bounds[i][1]
            else:
                bounds[1][0] = bounds[0][0]
            del bounds[i]
        else:
            i += 1

    os.makedirs(a.outdir, exist_ok=True)
    for idx, (lo, hi) in enumerate(bounds):
        path = os.path.join(a.outdir, "chunk_%04d.fa.gz" % idx)
        with gzip.open(path, "wt") as out:
            for name, _len, lines in recs[lo:hi]:
                out.write(name)
                out.writelines(lines)
        print("chunk_%04d\tcontigs=%d\tbp=%d" % (idx, hi - lo, chunk_bp([lo, hi])))


main()
