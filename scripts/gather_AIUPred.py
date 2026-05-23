#!/usr/bin/env python3

# Process AIupred results to call IDPs
# Definition of IDP region will be a stretch 
# A region that have at least 30 consecutive residues with iupred score over 
# 0.5 would be determined as putative disordered region (IDR).
# Any protein that have at least 1 IDR would be classified as IDP.
import os
import csv
import sys
import re
import time
import gzip
import argparse
from itertools import groupby
def average(lst):
    """
    Calculate the average of a list of numbers.
    """
    if len(lst) == 0:
        return 0
    return sum(lst) / len(lst)

def scores_to_idpstatus(iupred_scores, min_length=30):
    disordered_residues = 0
    all_residues = 0
    if max([sum([1 for _ in y]) if x == 1 else 0 for x, y in groupby([1 if x>=0.5 else 0 for x in iupred_scores])]) > 30:
        disordered_protein = 1
        disordered_residues += sum([1 for x in iupred_scores if x>=0.5])
    all_residues += len(iupred_scores)
    return disordered_residues, all_residues

def scores_to_idp_regions(iupred_scores, min_length=30):
    """
    Convert IUPred scores to IDP regions.
    """
    idp_regions = []
    current_start = None
    for i, score in enumerate(iupred_scores):
        if score >= 0.5:
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                s = current_start
                e = i - 1
                mean_score = average(iupred_scores[s:e+1])
                # 1 based not zero based
                idp_regions.append((s+1, e+1, e - s + 1, mean_score))
                current_start = None
    if current_start is not None:
        s = current_start
        e = len(iupred_scores) - 1
        mean_score = average(iupred_scores[s:e+1])
        # 1 based not 0 based
        idp_regions.append((s+1, e+1, e - s + 1, mean_score))
    return idp_regions

def parse_iupred_file(iupred_file):
    """
    Parse IUPred file and return a list of scores.
    """
    iupred_scores = []
    
    with gzip.open(iupred_file, 'rt') as f:
        seqname = None
        iupred_scores = []
        score_set = {}
        for line in f:
            if line.startswith("#>"):
                if seqname is not None:
                    score_set[seqname] = [
                        scores_to_idpstatus(iupred_scores),
                        scores_to_idp_regions(iupred_scores)
                        ]
                    iupred_scores = []
                seqname = line[2:].split()[1]
            elif line.startswith("#"):
                # skip header/comment
                continue
            else:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                try:
                    score = float(parts[2])
                    iupred_scores.append(score)
                except ValueError:
                    continue
        if seqname is not None:
            score_set[seqname] = [
                        scores_to_idpstatus(iupred_scores),
                        scores_to_idp_regions(iupred_scores)
                        ]
    return score_set

def main():
    if len(sys.argv) < 2:
        print("Usage: python gather_AIUPred.py iupred.ouput.txt.gz or -d <dir>")
        sys.exit(1)
        
    parser = argparse.ArgumentParser(description="Collect precomputed AIUPred results into a table",
                                    epilog='Example: gather_AIUPred.py')
    parser.add_argument("iupred_file", nargs='*', help="Input IUPred/AIUpred result file(s)")
    parser.add_argument("-d", "--dir", help="Input dir")
    parser.add_argument("-ext", "--ext", default="iupred.txt.gz",help="file extension when reading folder")
    parser.add_argument("--outfile", default='bigquery/idp.csv', 
                        help="Output IDP region table file")
    parser.add_argument("--outfilesum", default='bigquery/idp_summary.csv', 
                        help="Output IDP summary table file")
    parser.add_argument("--idp_length", default=30, type=int,
                        help="Minimum length of IDP region")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    with open (args.outfile, "w",newline="") as outfh, open (args.outfilesum, "w",newline="") as outfhsum:
        outwriter = csv.writer(outfh)
        outwriter.writerow(["species_prefix","protein_id","IDP_start","IDP_end", "IDP_length", "mean_score"])
        outwritersum = csv.writer(outfhsum)
        outwritersum.writerow(["species_prefix","protein_id","IDP_residues","IDP_fraction","length"])

        if args.iupred_file:
            for file in args.iupred_file:
                if args.debug:
                    print(f"Processing {file}")
                timestart = time.time() 
                iupred_regions = parse_iupred_file(file)
                timeend = time.time()
                if args.debug:
                    print(f"Reading Time elapsed: {timeend-timestart}")
                timestart = time.time()
                for seqid, idpdata in iupred_regions.items():
                    idp_counts, idp_regions = idpdata
                    idp_residue_count, all_residue_count = idp_counts
                    if idp_residue_count == 0:
                        if args.debug:
                            print(f"Skipping {seqid} with no IDP regions")
                        continue
                    prefix = seqid.split("_")[0]                    
                    outwritersum.writerow([prefix,seqid,idp_residue_count, f"{idp_residue_count/all_residue_count:.4f}", all_residue_count ])
                    for idp in idp_regions:
                        start, end, peplength, mean_score = idp
                        if peplength >= args.idp_length:
                            outwriter.writerow([prefix, seqid, start, end, peplength, f"{mean_score:.4f}"])
                timeend = time.time()
                if args.debug:
                    print(f"Writing Time elapsed: {timeend-timestart}")
        else:
            print("not ready to do folder processing, usually this is slow")
if __name__ == "__main__":
    main()