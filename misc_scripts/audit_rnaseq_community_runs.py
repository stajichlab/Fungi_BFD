#!/usr/bin/env python3
"""Flag species whose selected RNA-seq runs are community (metatranscriptome) samples.

Why this exists
---------------
SRA_QUERY / SRA_QUERY_BATCH pick RNA-seq runs by TaxID + LibraryStrategy and
rank the newest first. Some studies file one run per microbe taxon under that
taxon's name even though every run is a whole-community sample. Confirmed
2026-09-24: two Sichuan Agricultural University goat-rumen studies
(SRP500439/PRJNA1097788, SRP606484/PRJNA1301529; "Microbe sample from
<species>", SampleName like H.RUEMN.5.1) list 12 taxa each, 6 of them fungi.
In the Funannotate benchmark, Rhizopus microsporus got 100% rumen reads
(26 of 49.2M reads mapped; Trinity-GG 0 transcripts) and Batrachochytrium
dendrobatidis 86%.

The query modules now drop these runs (LibrarySource METATRANSCRIPTOMIC /
METAGENOMIC, rumen keywords, study denylist). But the sra_query/*.csv caches
are storeDir outputs written before that fix, and SRA keeps growing. This
script is the monitoring check: re-run it after any SRA_QUERY refresh.

What it does
------------
1. Fetches from ENA every run with tax_tree(4751) (Fungi), library_strategy
   RNA-Seq and library_source METATRANSCRIPTOMIC or METAGENOMIC, plus every
   run of the denylisted studies.
2. Reads every <species>.sra_query.csv cache and reports each selected run
   that is in that set, with its rank in the species' list, study, center
   and whether the species' _norm_R1.fastq.gz read file is non-empty (reads
   were already downloaded and may already feed Trinity/PASA).

Output: a TSV (default rnaseq_community_run_audit.tsv) and a summary on
stderr. Exit code 0 always; the TSV is the monitoring artifact.

Usage
-----
    python3 audit_rnaseq_community_runs.py \\
        --sra-query-dir <launchDir>/rnaseq_reads/sra_query \\
        --reads-dir     <launchDir>/rnaseq_reads \\
        [--count-reads --threads 8] --out rnaseq_community_run_audit.tsv

    Fungi_BFD: <launchDir> = /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
    nf_funannotate1 / benchmark cells: <launchDir> = the run directory.

Kept byte-identical in Fungi_BFD (misc_scripts/) and nf_funannotate1
(scripts/); change both together.
"""
from __future__ import annotations

import argparse
import csv
import glob
import io
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import urllib.parse
import urllib.request

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,study_accession,secondary_study_accession,library_source,center_name,scientific_name,sample_title"

# Studies confirmed to file community samples under individual taxon names.
DENYLIST_STUDIES = {
    "PRJNA1097788": "SRP500439 Sichuan Agricultural University goat rumen metatranscriptome",
    "PRJNA1301529": "SRP606484 Sichuan Agricultural University goat rumen metatranscriptome",
}


def ena_query(query: str) -> list[dict]:
    params = urllib.parse.urlencode({"result": "read_run", "query": query, "fields": FIELDS,
                                     "format": "tsv", "limit": 0})
    with urllib.request.urlopen(f"{ENA_SEARCH}?{params}", timeout=300) as r:
        text = r.read().decode()
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def flagged_runs() -> dict[str, dict]:
    flagged = {}
    for src in ("METATRANSCRIPTOMIC", "METAGENOMIC"):
        rows = ena_query(f'tax_tree(4751) AND library_strategy="RNA-Seq" AND library_source="{src}"')
        print(f"[INFO] ENA fungal RNA-Seq runs with library_source={src}: {len(rows)}", file=sys.stderr)
        for r in rows:
            r["reason"] = f"library_source={src}"
            flagged[r["run_accession"]] = r
    for study, note in DENYLIST_STUDIES.items():
        rows = ena_query(f'study_accession="{study}"')
        print(f"[INFO] denylisted study {study}: {len(rows)} runs", file=sys.stderr)
        for r in rows:
            prev = flagged.get(r["run_accession"])
            r["reason"] = (prev["reason"] + ";" if prev else "") + f"denylist:{study}"
            flagged[r["run_accession"]] = r
    return flagged


def run_read_counts(r1: str) -> dict[str, int]:
    """Count reads per run accession in a _norm_R1.fastq.gz (header '@SRR123.45/1')."""
    cmd = f"zcat '{r1}' | awk 'NR%4==1{{split(substr($1,2),a,\".\");c[a[1]]++}} END{{for(k in c)print k\"\\t\"c[k]}}'"
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True).stdout
    return {k: int(v) for k, v in (l.split("\t") for l in out.splitlines() if l)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sra-query-dir", required=True, help="directory of <species>.sra_query.csv caches")
    ap.add_argument("--reads-dir", default=None, help="directory of <species>_norm_R1.fastq.gz (optional)")
    ap.add_argument("--out", default="rnaseq_community_run_audit.tsv")
    ap.add_argument("--count-reads", action="store_true",
                    help="for flagged species, count reads per run in the R1 file (full zcat; slow)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    flagged = flagged_runs()
    if not flagged:
        print("[ERROR] ENA returned no flagged runs -- API problem? Refusing to report a clean audit.", file=sys.stderr)
        return 1

    caches = sorted(glob.glob(os.path.join(args.sra_query_dir, "*.sra_query.csv")))
    print(f"[INFO] sra_query caches: {len(caches)}", file=sys.stderr)
    if not caches:
        print(f"[ERROR] no *.sra_query.csv under {args.sra_query_dir}", file=sys.stderr)
        return 1

    hits = []
    n_runs = 0
    for path in caches:
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        n_runs += len(rows)
        for rank, row in enumerate(rows, 1):
            acc = row.get("sra_accession", "")
            if acc not in flagged:
                continue
            f = flagged[acc]
            reads = ""
            if args.reads_dir:
                r1 = os.path.join(args.reads_dir, f"{row['species_tag']}_norm_R1.fastq.gz")
                reads = "downloaded" if os.path.exists(r1) and os.path.getsize(r1) > 0 else "not_downloaded"
            hits.append({
                "species_tag": row["species_tag"], "sra_accession": acc, "rank_in_cache": rank,
                "n_runs_in_cache": len(rows), "reason": f["reason"], "study": f["study_accession"],
                "srp": f["secondary_study_accession"], "center": f["center_name"],
                "ena_scientific_name": f["scientific_name"], "sample_title": f["sample_title"],
                "reads_file": reads,
            })

    if args.count_reads and args.reads_dir:
        todo = sorted({h["species_tag"] for h in hits if h["reads_file"] == "downloaded"})
        print(f"[INFO] counting reads per run for {len(todo)} species", file=sys.stderr)
        paths = {s: os.path.join(args.reads_dir, f"{s}_norm_R1.fastq.gz") for s in todo}
        with ThreadPoolExecutor(args.threads) as ex:
            counts = dict(zip(todo, ex.map(run_read_counts, [paths[s] for s in todo])))
        for h in hits:
            c = counts.get(h["species_tag"])
            if c is None:
                continue
            tot = sum(c.values())
            h["reads_from_run"] = c.get(h["sra_accession"], 0)
            h["reads_total"] = tot
            h["pct_reads_from_run"] = f"{100 * h['reads_from_run'] / tot:.1f}" if tot else ""
            h["reads_file"] = "contains_run" if h["reads_from_run"] else "run_not_in_file"

    cols = ["species_tag", "sra_accession", "rank_in_cache", "n_runs_in_cache", "reason", "study", "srp",
            "center", "ena_scientific_name", "sample_title", "reads_file", "reads_from_run", "reads_total",
            "pct_reads_from_run"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(hits)

    species = sorted({h["species_tag"] for h in hits})
    dl = sorted({h["species_tag"] for h in hits if h["reads_file"] in ("downloaded", "contains_run")})
    print(f"[INFO] runs checked: {n_runs} in {len(caches)} caches", file=sys.stderr)
    print(f"[RESULT] flagged runs: {len(hits)}; species affected: {len(species)}"
          + (f"; with reads already downloaded: {len(dl)}" if args.reads_dir else ""), file=sys.stderr)
    for s in species:
        accs = ",".join(h["sra_accession"] for h in hits if h["species_tag"] == s)
        pct = sum(float(h.get("pct_reads_from_run") or 0) for h in hits if h["species_tag"] == s)
        tag = (f"\t{pct:.1f}% of R1 reads from flagged runs" if args.count_reads else "\tDOWNLOADED") if s in dl else ""
        print(f"  {s}\t{accs}{tag}", file=sys.stderr)
    print(f"[INFO] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
