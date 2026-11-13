#!/usr/bin/env python3
"""merge_filter.py -- which result files belong to genomes in samples.csv
minus suppress.txt. Python twin of loadMergeAllowedKeys()/mergeFileAllowed()
in nextflow/modules/common/utils.nf; keep the two in sync (tests:
nextflow/tests/test_merge_filter.py).

As a script it is a report-only orphan scan: it lists files under the given
result folders that a merge would drop, and never moves or deletes anything.

Usage:
    merge_filter.py --samples samples.csv --suppress suppress.txt \\
        RESULTS_DIR [RESULTS_DIR ...] [-o orphans.tsv]

Output TSV: tool, file name, path, bytes.
"""
import argparse
import csv
import os
import re
import sys

TARGETP_SUFFIX = "_summary.targetp2.gz"


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_representative_strain import make_sample_tag  # noqa: E402  (mirror of makeSampleTag)


def load_suppress(path):
    if not path or not os.path.isfile(path):
        return set()
    out = set()
    for line in open(path):
        k = line.strip().split(",")[0].strip()
        if k and not k.startswith("#"):
            out.add(k)
    return out


def allowed_keys(samples, suppress_path):
    sup = load_suppress(suppress_path)
    keys = set()
    for r in csv.DictReader(open(samples)):
        asmid = (r.get("ASMID") or "").strip()
        if not asmid or asmid in sup:
            continue
        keys.add(asmid)
        tag = (r.get("LOCUSTAG") or "").strip()
        if tag:
            keys.add(tag)
        keys.add(make_sample_tag((r.get("SPECIES") or "").strip(), (r.get("STRAIN") or "").strip()))
    return keys


def file_allowed(keys, name):
    """True when a prefix of `name` ending just before a '.' is an allowed key,
    or name is <key>_summary.targetp2.gz. A '_' boundary is not enough."""
    for m in re.finditer(r"\.", name):
        i = m.start()
        if i > 0 and name[:i] in keys:
            return True
    return name.endswith(TARGETP_SUFFIX) and name[: -len(TARGETP_SUFFIX)] in keys


def main():
    ap = argparse.ArgumentParser(description="Report result files a merge would drop (report only).")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--suppress", default=None)
    ap.add_argument("dirs", nargs="+", help="e.g. results/function results/genome_stats")
    ap.add_argument("-o", "--out", default="-")
    a = ap.parse_args()
    keys = allowed_keys(a.samples, a.suppress)
    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["tool", "file", "path", "bytes"])
    counts = {}
    for d in a.dirs:
        for tool in sorted(os.listdir(d)):
            top = os.path.join(d, tool)
            if not os.path.isdir(top):
                continue
            for root, _dirs, files in os.walk(top):
                for f in files:
                    if not file_allowed(keys, f):
                        p = os.path.join(root, f)
                        w.writerow([tool, f, p, os.path.getsize(p)])
                        counts[tool] = counts.get(tool, 0) + 1
    for tool, n in sorted(counts.items()):
        print(f"{tool}: {n} files not in samples.csv or suppressed", file=sys.stderr)
    if not counts:
        print("no orphan files", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
