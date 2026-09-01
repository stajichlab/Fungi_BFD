#!/usr/bin/env python3
"""Collect WGD per-step wall-time measurements from a paralogoscope (resume) run.

Reads Nextflow work-dir metadata (.command.begin, .command.log mtime, .exitcode)
plus the published dmd families TSVs, and emits wgd_perf_measurements.tsv.

Usage:
    collect_measurements.py RUN_DIR OUT_TSV [--cds-dir DIR]

RUN_DIR  : the top-level dir holding work_r3/ and out_r3/ (the -w work dir
           parent and the --outdir of the test run, respectively).
OUT_TSV  : path to write the measurements table.
"""

import argparse
import glob
import os
import re
import sys

def parse_begin(p):
    try:
        with open(p) as fh:
            c = fh.read().strip()
        if c:
            return c
        import datetime
        return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return None

def log_mtime(d):
    p = os.path.join(d, ".command.log")
    try:
        import datetime
        return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return None

def family_stats(fams_tsv):
    """One line = one family: GFxxxx \\t gene1, gene2, ..."""
    sizes = []
    if not os.path.exists(fams_tsv):
        return None
    for line in open(fams_tsv):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2 or parts[0] == "":
            continue
        sizes.append(len([x for x in parts[1].split(",") if x.strip()]))
    if not sizes:
        return None
    return {
        "n_families": len(sizes),
        "n_multi": sum(1 for s in sizes if s >= 2),
        "max_family_size": max(sizes),
        "genes_in_families": sum(sizes),
    }

def walk_work(workdir):
    """Return: list of (dir, script_text) for tasks that have a .command.sh."""
    found = []
    for cmd in glob.glob(os.path.join(workdir, "*", "*", ".command.sh")):
        try:
            txt = open(cmd, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        found.append((os.path.dirname(cmd), txt))
    return found

def classify(script, ksd_log):
    if "wgd ksd" in script:
        return "ksd"
    if "wgd dmd" in script:
        return "dmd"
    if "wgd syn" in script:
        return "syn"
    return "?"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("out_tsv")
    ap.add_argument("--cds-dir", default=None, help="input/cds dir to count genes per genome")
    args = ap.parse_args()

    rows = []
    for d, script in walk_work(os.path.join(args.run_dir, "work_r3")):
        step = classify(script, None)
        if step == "?":
            continue
        beg = parse_begin(os.path.join(d, ".command.begin"))
        mlog = log_mtime(d)
        try:
            rc = open(os.path.join(d, ".exitcode")).read().strip()
        except OSError:
            rc = ""
        # sample id: search for <id>.cds-transcripts.fa in script or dir list
        samp = None
        m = re.search(r"([A-Za-z0-9._-]+)\.cds-transcripts\.fa", script)
        if m:
            samp = m.group(1)
        # ksd progress from log
        n_fam_log = ksd_families(d) if step == "ksd" else None
        # published outputs
        outs = [f for f in os.listdir(d) if not f.startswith(".")]
        rows.append({
            "task_dir": d, "step": step, "sample": samp, "exitcode": rc,
            "begin": beg, "end_log_mtime": mlog, "ksd_families_logged": n_fam_log,
            "outputs": ";".join(sorted(outs[:6])),
        })

    # attach family stats + gene counts from out_r3 publish dirs
    outdir = os.path.join(args.run_dir, "out_r3")
    fam_by_samp = {}
    for fam in glob.glob(os.path.join(outdir, "*", "wgd_dmd", "*.tsv")):
        samp = os.path.basename(os.path.dirname(os.path.dirname(fam)))
        fs = family_stats(fam)
        if fs:
            fam_by_samp[samp] = fs

    genes_by_samp = {}
    if args.cds_dir:
        wanted = {r["sample"] for r in rows if r["sample"]}
        for base in wanted:
            fa = os.path.join(args.cds_dir, base + ".cds-transcripts.fa")
            if not os.path.exists(fa):
                continue
            n = 0
            for line in open(fa):
                if line.startswith(">"):
                    n += 1
            genes_by_samp[base] = n

    with open(args.out_tsv, "w") as fh:
        cols = ["sample", "step", "task_dir", "exitcode", "begin", "end_log_mtime",
                "n_genes", "n_families", "n_multi_families", "max_family_size",
                "genes_in_families", "ksd_families_logged", "outputs"]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["step"], r["begin"] or "")):
            samp = r["sample"] or ""
            fs = fam_by_samp.get(samp, {})
            vals = [samp, r["step"], r["task_dir"].replace(args.run_dir + "/", ""),
                    r["exitcode"], r["begin"], r["end_log_mtime"],
                    genes_by_samp.get(samp, ""), fs.get("n_families", ""),
                    fs.get("n_multi", ""), fs.get("max_family_size", ""),
                    fs.get("genes_in_families", ""), r["ksd_families_logged"] or "",
                    r["outputs"]]
            fh.write("\t".join(str(v) for v in vals) + "\n")
    print(f"wrote {args.out_tsv}")

def ksd_families(d):
    p = os.path.join(d, ".command.log")
    n = 0
    try:
        for line in open(p, errors="replace"):
            if "Analysing family" in line:
                n += 1
    except OSError:
        pass
    return n

if __name__ == "__main__":
    main()
