#!/usr/bin/env python3

import os
import csv
import argparse
import sys
from Bio import SeqIO

def load_samples(fh, line_start=None, set_size=1):
    samples = []
    reader = csv.DictReader(fh)
    n = 0
    for row in reader:
        if line_start is None or n >= line_start:
            samples.append(row)
            if line_start is not None and n+1 >= (line_start + set_size):
                break
        n += 1
    return samples

def main():
    parser = argparse.ArgumentParser(description="Collect per chrom stats into a table to support gene density computation",
                                    epilog='Example: collect_chrom_info.py')
    parser.add_argument("-d","--genomedir", default="genomes", help="Directory with genomes *.scaffolds.fa files")
    parser.add_argument("--samples", default="samples.csv", help="samples.csv file for fungi5k")
    parser.add_argument("-o","--outfile", default="bigquery/chrom_info.csv", 
                        help="output file [bigquery/chrom_info.csv]")
    parser.add_argument("--run_with", type=int, required=False, default=None,
                        metavar="LINE_NUMBER",
                        help="run from this starting line from samples file")
    parser.add_argument("--run_set", type=int, required=False, default=1,
                        metavar="LINE_SET",
                        help="run only this many from samples file")

    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    with open(args.samples,"r") as fh, open(args.outfile, "w",newline="") as outfh:
        species = load_samples(fh,args.run_with, args.run_set)
        csvout = csv.writer(outfh)        
        # potentially add RIP index stat here
        header = [ 'LOCUSTAG', 'chrom_name', 'length', 'GC_percent', 'GC_count', 'left_50', 'right_50', 'lower_masked', 'N_masked' ]
        csvout.writerow(header)
            
        for sp in species:
            species_string = sp['SPECIES']
            if len(sp['STRAIN']):
                species_string += "_" + sp['STRAIN']
            species_string = species_string.replace(' ', '_')
            stemname = f"{species_string}.scaffolds.fa"
            genomefile = os.path.join(args.genomedir,stemname)
            
            # remove strain for the edge case of A. niger 
            if not os.path.exists(genomefile):
                species_string = sp['SPECIES'].replace(' ', '_')
                stemname = f"{species_string}.scaffolds.fa"
                genomefile = os.path.join(args.genomedir,stemname)
                
                if not os.path.exists(genomefile):
                    species_string = sp['SPECIESIN']
                    if len(sp['STRAIN']):
                        species_string += "_" + sp['STRAIN']
                    species_string = species_string.replace(' ', '_')                
                    stemname = f"{species_string}.scaffolds.fa"
                    genomefile = os.path.join(args.genomedir,stemname)
                    
                    if not os.path.exists(genomefile):
                        species_string = sp['SPECIESIN'].replace(' ', '_')                
                        stemname = f"{species_string}.scaffolds.fa"
                        genomefile = os.path.join(args.genomedir,stemname)
                        if not os.path.exists(genomefile):
                            print(f"Missing the genomefile (SPECIESIN) {genomefile}")
                            continue            
            if not os.path.exists(genomefile):
                if args.debug:
                    print("Missing the genomefile: ",genomefile)
                continue

            if args.debug:
                print(genomefile,"exists")
                continue

            for record in SeqIO.parse(genomefile, "fasta"):
                chrom_length = len(record.seq)
                seqstr_orig = str(record.seq)
                lowercase = seqstr_orig.count("a") + seqstr_orig.count("g") + seqstr_orig.count("c") + seqstr_orig.count("t")
                seqstr = seqstr_orig.upper()
                GC = seqstr.count("G") + seqstr.count("C")
                row = [ sp['LOCUSTAG'], 
                        record.id, 
                        chrom_length, 
                        round(GC / chrom_length,4),
                        GC,
                        seqstr[:50],
                        seqstr[-50:],
                        lowercase,
                        seqstr.count("N") 
                        ]
                csvout.writerow(row)


if __name__ == "__main__":
    main()