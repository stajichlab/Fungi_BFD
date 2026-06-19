#!/usr/bin/env python3
"""
sourmash_matrix_to_long.py — Convert a `sourmash compare --ani --csv` square
matrix into the long `query<TAB>reference<TAB>ANI` form consumed by report_ani.py.

`sourmash compare --csv` writes:
  - a header row of N signature labels
  - N rows of N similarity values, in the same order as the header
    (there is no row-label column)

With --ani the values are ANI estimates in the range 0..1; we multiply by 100 to
match the percentage scale used by fastANI/skani/mash. Self-comparisons (the
diagonal) and the lower triangle are skipped — only i<j pairs are emitted.

Signature labels are basename-d so they match the genome filenames in the
names TSV (compare_ANI.nf sketches each genome with --name <genome_filename>).
"""

import argparse
import csv
import sys


def basename(label):
    return label.rsplit('/', 1)[-1].strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True, help='sourmash compare --csv matrix')
    ap.add_argument('--output', required=True, help='output long-form TSV (q r ANI)')
    args = ap.parse_args()

    with open(args.input, newline='') as fh:
        reader = csv.reader(fh)
        rows = [row for row in reader if row]

    if not rows:
        open(args.output, 'w').close()
        print("sourmash_matrix_to_long: empty matrix", file=sys.stderr)
        return

    labels = [basename(x) for x in rows[0]]
    data = rows[1:]
    n = len(labels)

    n_pairs = 0
    with open(args.output, 'w') as out:
        for i, row in enumerate(data):
            if i >= n:
                break
            for j in range(i + 1, min(len(row), n)):
                try:
                    val = float(row[j])
                except ValueError:
                    continue
                ani = val * 100.0 if val <= 1.0 else val
                out.write(f"{labels[i]}\t{labels[j]}\t{ani:.4f}\n")
                n_pairs += 1

    print(f"sourmash_matrix_to_long: {n} genomes -> {n_pairs} pairs",
          file=sys.stderr)


if __name__ == '__main__':
    main()
