#!/usr/bin/env python3
"""Merge per-species TargetP summary results into a tables-loadable CSV."""

import csv
import gzip
import re
import sys
import argparse

CS_RE = re.compile(r'CS pos:\s+(\d+)-(\d+)\.\s+(\S+)\.\s+Pr:\s+(\S+)')

def main():
    parser = argparse.ArgumentParser(description="Merge TargetP results")
    parser.add_argument("summaries", nargs="+",
                        help="*_summary.targetp2.gz files")
    parser.add_argument("-o", "--outfile", default="targetP.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'prediction', 'probability',
                    'cleavage_position_start', 'cleavage_position_end',
                    'cleavage_probability', 'motif'])
        for f in sorted(args.summaries):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t", 2)
                    if len(parts) < 3:
                        continue
                    pid, prediction, rest = parts
                    if prediction == "noTP":
                        continue
                    fields = rest.split("\t")
                    SP  = fields[1] if len(fields) > 1 else ""
                    mTP = fields[2] if len(fields) > 2 else ""
                    CS  = fields[3] if len(fields) > 3 else ""
                    prob = SP if prediction == "SP" else (mTP if prediction == "mTP" else "")
                    m = CS_RE.match(CS)
                    cs_start, cs_end, motif, cs_prob = "0", "0", "", "0.0"
                    if m:
                        cs_start, cs_end, motif, cs_prob = m.groups()
                    w.writerow([pid.split("_")[0], pid, prediction, prob,
                                cs_start, cs_end, cs_prob, motif])

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
