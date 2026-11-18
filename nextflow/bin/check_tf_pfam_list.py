#!/usr/bin/env python3
"""Check assets/fungal_tf_pfam_domains.csv against a Pfam release.

For each CSV row, look up the accession in Pfam-A.hmm.dat (the `#=GF ID` and
`#=GF AC` lines). The version suffix is removed (PF00172.25 -> PF00172).
The check fails (exit 1) when the accession is absent from the dat file
or when `pfam_name` differs from the Pfam ID.

A missing --pfam-dat path given on the command line is an error (exit 2).
A missing default path is logged and skipped (exit 0). The default can be
absent in containers.

Pfam-A.hmm.dat may be plain text or .gz.
"""
import argparse
import csv
import gzip
import os
import sys

DEFAULT_DAT = "/srv/projects/db/pfam/current/Pfam-A.hmm.dat"
DEFAULT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           os.pardir, "assets", "fungal_tf_pfam_domains.csv")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, newline="")


def load_pfam_ids(dat_path):
    """Return {accession without version: Pfam ID}."""
    ids = {}
    cur_id = None
    with open_text(dat_path) as fh:
        for line in fh:
            if line.startswith("#=GF ID"):
                cur_id = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC") and cur_id is not None:
                acc = line.split(None, 2)[2].strip().split(".")[0]
                ids[acc] = cur_id
                cur_id = None
    return ids


def check(csv_path, ids):
    """Return (report lines, number of failures)."""
    lines = []
    failures = 0
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            acc = row["pfam_acc"].split(".")[0]
            name = row["pfam_name"]
            dat_id = ids.get(acc)
            if dat_id is None:
                status = "FAIL absent from Pfam dat"
                failures += 1
            elif dat_id != name:
                status = "FAIL name differs"
                failures += 1
            else:
                status = "ok"
            lines.append("\t".join([acc, name, dat_id or "-", status]))
    return lines, failures


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--csv", default=DEFAULT_CSV, help="TF domain list CSV")
    ap.add_argument("--pfam-dat", default=None,
                    help="Pfam-A.hmm.dat (default: %s)" % DEFAULT_DAT)
    args = ap.parse_args(argv)

    dat = args.pfam_dat
    if dat is None:
        dat = DEFAULT_DAT
        if not os.path.exists(dat):
            print("check_tf_pfam_list: default Pfam dat not found (%s); skipped" % dat,
                  file=sys.stderr)
            return 0
    elif not os.path.exists(dat):
        print("check_tf_pfam_list: ERROR: --pfam-dat not found: %s" % dat, file=sys.stderr)
        return 2

    ids = load_pfam_ids(dat)
    lines, failures = check(args.csv, ids)
    print("pfam_acc\tcsv_name\tdat_id\tstatus")
    print("\n".join(lines))
    if failures:
        print("check_tf_pfam_list: %d row(s) failed" % failures, file=sys.stderr)
        return 1
    print("check_tf_pfam_list: %d rows ok" % len(lines), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
