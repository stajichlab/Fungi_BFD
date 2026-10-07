#!/usr/bin/env python3
"""List genomes whose last task attempt FAILED, from Nextflow trace file(s).

Used with errorStrategy 'ignore' after the last retry (conf/profile_funannotate.config): a
failed genome is dropped and the wave continues, so failures are only visible in the trace.
A (process, tag) pair counts as failed when its latest trace row is FAILED (a later COMPLETED
or CACHED row, for example after a -resume, clears it). Reads plain or .gz traces.

Usage: list_failed_genomes.py logs/nextflow/funannotate_trace.*.txt > failed_genomes.tsv
"""
import csv
import gzip
import sys


def rows(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main(paths):
    last = {}
    for path in paths:  # later files and later rows win
        for r in rows(path):
            proc = r["name"].split(" (")[0].split(":")[-1]
            last[(proc, r.get("tag") or r["name"])] = r
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["process", "genome", "status", "exit", "task_id", "hash"])
    n = 0
    for (proc, tag), r in sorted(last.items()):
        if r["status"] == "FAILED":
            w.writerow([proc, tag, r["status"], r["exit"], r["task_id"], r["hash"]])
            n += 1
    print(f"{n} failed genome/process pairs", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
