#!/usr/bin/env python3
"""Merge per-species SwissProt homology (DIAMOND blastp --sensitive, or NCBI
blastp --outfmt 6) results into a tables-loadable CSV.

Input rows are the 18-column format shared by DIAMOND (-f 6) and blastp
(-outfmt "6 ..."):

    qseqid sseqid pident positive nident length mismatch gapopen
    qstart qend sstart send evalue bitscore qcovhsp qlen slen stitle

The functional-transfer 80-80 flag (Rost 1999) is derived rather than
pre-filtered so the raw evidence is always retained:
    pident >= 80 AND query_cov >= 0.8 AND hit_cov >= 0.8
where query_cov = length/qlen and hit_cov = length/slen.

Rich annotation (protein name, GO, EC, InterPro, Pfam) lives in
swissprot_annot.parquet, keyed by accession; join on swissprot_acc.
"""

import csv
import gzip
import sys
import argparse

# qseqid sseqid pident positive nident length mismatch gapopen qstart qend
# sstart send evalue bitscore qcovhsp qlen slen stitle
COL_QSEQ, COL_SSEQ, COL_PIDENT, COL_POSITIVE, COL_NIDENT, COL_LENGTH = range(6)
COL_MISMATCH, COL_GAPOPEN = 6, 7
COL_QSTART, COL_QEND, COL_SSTART, COL_SEND = 8, 9, 10, 11
COL_EVALUE, COL_BITSCORE = 12, 13
COL_QCOVHSP, COL_QLEN, COL_SLEN = 14, 15, 16
COL_STITLE = 17

HEADER = [
    'species_prefix', 'protein_id', 'swissprot_acc', 'swissprot_entry',
    'swissprot_name', 'pident', 'positive', 'nident', 'aln_length',
    'q_len', 's_len', 'qcovhsp', 'query_cov', 'hit_cov', 'q_start',
    'q_end', 's_start', 's_end', 'gapopen', 'mismatch', 'evalue',
    'bitscore', 'func_transfer_80_80',
]


def swissprot_name_from_stitle(stitle):
    # "sp|P40231|CSK2A_SCHPO Casein kinase II subunit alpha OS=... GN=..."
    desc = stitle.split(' OS=')[0]
    parts = desc.split(' ')
    # Drop the leading "sp|ACC|ENTRY" token.
    return ' '.join(parts[1:]) if len(parts) > 1 else ''


def main():
    parser = argparse.ArgumentParser(description="Merge SwissProt blasttab results")
    parser.add_argument("blasttabs", nargs="+", help="*.blasttab.gz files")
    parser.add_argument("-o", "--outfile", default="swissprot.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(HEADER)
        for f in sorted(args.blasttabs):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row or len(row) < COL_STITLE + 1:
                        continue
                    qseqid = row[COL_QSEQ]
                    sseqid = row[COL_SSEQ]
                    sacc = sseqid.split('|')[1] if '|' in sseqid else sseqid
                    sbits = sseqid.split('|')
                    sentry = sbits[2] if len(sbits) > 2 else ''
                    aln_len = float(row[COL_LENGTH])
                    q_len = float(row[COL_QLEN])
                    s_len = float(row[COL_SLEN])
                    query_cov = aln_len / q_len if q_len else 0.0
                    hit_cov = aln_len / s_len if s_len else 0.0
                    pident = float(row[COL_PIDENT])
                    func80 = 1 if (pident >= 80.0 and query_cov >= 0.8
                                   and hit_cov >= 0.8) else 0
                    w.writerow([
                        qseqid.split('_')[0], qseqid, sacc, sentry,
                        swissprot_name_from_stitle(row[COL_STITLE]),
                        row[COL_PIDENT], row[COL_POSITIVE], row[COL_NIDENT],
                        row[COL_LENGTH], row[COL_QLEN], row[COL_SLEN],
                        row[COL_QCOVHSP], f"{query_cov:.4f}", f"{hit_cov:.4f}",
                        row[COL_QSTART], row[COL_QEND], row[COL_SSTART],
                        row[COL_SEND], row[COL_GAPOPEN], row[COL_MISMATCH],
                        row[COL_EVALUE], row[COL_BITSCORE], func80,
                    ])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
