#!/usr/bin/env python3
"""Census of genemark_store sizes across a runs tree (read-only).

Finds <runs>/genemark_store and <runs>/*/genemark_store. For each store it counts species
folders and sums bytes per file type: .genemark.gtf, .genemark.gtf.gz, .genemark.mod,
.other.gff3, and anything else. A store reached by two paths (symlink) is counted once, by
real path. Sizes are stat() sizes (apparent size), not disk blocks.

Output: TSV to stdout, one row per store plus a TOTAL row. Progress goes to stderr.
Usage: genemark_store_census.py --runs /bigdata/.../Fungi_BFD_runs [--workers 16] > census.tsv
"""
import argparse
import concurrent.futures as cf
import glob
import os
import sys
import time

KINDS = ("gtf", "gtf_gz", "mod", "other_gff3", "other")


def kind(name):
    if name.endswith(".genemark.gtf"):
        return "gtf"
    if name.endswith(".genemark.gtf.gz"):
        return "gtf_gz"
    if name.endswith(".genemark.mod"):
        return "mod"
    if name.endswith(".other.gff3"):
        return "other_gff3"
    return "other"


def scan_species(path):
    """bytes and counts per kind for one species folder; symlinks are counted, not followed."""
    b = dict.fromkeys(KINDS, 0)
    n = dict.fromkeys(KINDS, 0)
    links = 0
    try:
        with os.scandir(path) as it:
            for e in it:
                if e.is_symlink():
                    links += 1
                    continue
                if e.is_file(follow_symlinks=False):
                    k = kind(e.name)
                    b[k] += e.stat(follow_symlinks=False).st_size
                    n[k] += 1
    except OSError as err:
        print(f"WARN {path}: {err}", file=sys.stderr)
    return b, n, links


def scan_store(store, workers):
    species = [e.path for e in os.scandir(store) if e.is_dir(follow_symlinks=False)]
    tb = dict.fromkeys(KINDS, 0)
    tn = dict.fromkeys(KINDS, 0)
    links = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        for b, n, l in ex.map(scan_species, species):
            for k in KINDS:
                tb[k] += b[k]
                tn[k] += n[k]
            links += l
    return len(species), tb, tn, links


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()

    stores = sorted(set(glob.glob(os.path.join(a.runs, "genemark_store")) + glob.glob(os.path.join(a.runs, "*", "genemark_store"))))
    seen, todo = {}, []
    for s in stores:
        rp = os.path.realpath(s)
        if rp in seen:
            print(f"SKIP {s}: same real path as {seen[rp]}", file=sys.stderr)
            continue
        seen[rp] = s
        if os.path.isdir(rp):
            todo.append(s)

    cols = ["store", "species_dirs", "symlinks"] + [f"bytes_{k}" for k in KINDS] + [f"n_{k}" for k in KINDS]
    print("\t".join(cols))
    tot = {"species": 0, "links": 0, "b": dict.fromkeys(KINDS, 0), "n": dict.fromkeys(KINDS, 0)}
    for s in todo:
        t0 = time.time()
        ns, tb, tn, links = scan_store(s, a.workers)
        print(f"{s}: {ns} species dirs in {time.time() - t0:.0f}s", file=sys.stderr)
        print("\t".join([os.path.relpath(s, a.runs), str(ns), str(links)] + [str(tb[k]) for k in KINDS] + [str(tn[k]) for k in KINDS]), flush=True)
        tot["species"] += ns
        tot["links"] += links
        for k in KINDS:
            tot["b"][k] += tb[k]
            tot["n"][k] += tn[k]
    print("\t".join(["TOTAL", str(tot["species"]), str(tot["links"])] + [str(tot["b"][k]) for k in KINDS] + [str(tot["n"][k]) for k in KINDS]))


if __name__ == "__main__":
    main()
