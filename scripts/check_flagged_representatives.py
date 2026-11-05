#!/usr/bin/env python3
"""
Cross-reference find_ani_label_mismatches.py's flagged ASMIDs against the
current representative picks (repr_assignments.tsv) -- the targeted-triage
guard used in place of full Option 1 (DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md)
for pre-launch checks. A flagged genome that never wins representative rank
cannot repeat the KCTC_13826BP failure mode (a mislabeled genome promoted to
representative, poisoning its species' shared Trinity / ab-initio reuse), so
this narrows the 172+/241+ dataset-wide mismatch flags down to the small set
that actually matters before an RNA-seq/train launch.

Usage:
    python scripts/check_flagged_representatives.py \\
        --mismatches ani_label_mismatches.tsv \\
        --predict-input genome_annotation/_reuse_assignments/predict_input_for_ani.tsv \\
        --repr-assignments genome_annotation/_reuse_assignments/repr_assignments.tsv

Exit code is the number of hits (flagged ASMID that IS a representative).
0 = clean to proceed; >0 = must resolve (relabel like KCTC, or implement
Option 1) before launching.
"""

import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mismatches", required=True,
                   help="output TSV from find_ani_label_mismatches.py")
    p.add_argument("--predict-input", required=True,
                   help="predict_input_for_ani.tsv (asmid -> out mapping)")
    p.add_argument("--repr-assignments", required=True,
                   help="repr_assignments.tsv (out, is_representative, ...)")
    p.add_argument("--outfile", default=None, help="write hits TSV here instead of stdout")
    return p.parse_args()


def main():
    args = parse_args()

    flagged = {}
    with open(args.mismatches, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for asmid in row["minority_asmids"].split(","):
                asmid = asmid.strip()
                if asmid:
                    flagged[asmid] = (row["species"], row["confidence"], row["max_ani_minority_to_majority"])

    asmid_to_out = {}
    with open(args.predict_input, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            a = (row.get("asmid") or "").strip()
            o = (row.get("out") or "").strip()
            if a and o:
                asmid_to_out[a] = o

    reps = set()
    with open(args.repr_assignments, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["is_representative"] == "True":
                reps.add(row["out"])

    hits = []
    for asmid, (species, conf, maxani) in flagged.items():
        out = asmid_to_out.get(asmid)
        if out and out in reps:
            hits.append((asmid, out, species, conf, maxani))

    out_fh = open(args.outfile, "w", newline="") if args.outfile else sys.stdout
    writer = csv.writer(out_fh, delimiter="\t")
    writer.writerow(["asmid", "out", "species", "confidence", "max_ani_minority_to_majority"])
    for h in hits:
        writer.writerow(h)
    if args.outfile:
        out_fh.close()

    print(f"[info] {len(flagged)} flagged ASMIDs, {len(hits)} currently a representative",
          file=sys.stderr)
    return len(hits)


if __name__ == "__main__":
    sys.exit(main())
