#!/usr/bin/env python3

from Bio import SeqIO
from collections import Counter
import sys
import argparse
import csv
import os
import time


def calculate_aa_frequencies(fasta_file):
    aa_counter = Counter()
    total_aa = 0
    species = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if not species:
            species = record.id.split("_")[0]
        aa_counter.update(record.seq)
        total_aa += len(record.seq)

    aa_frequencies = {aa: count / total_aa for aa, count in aa_counter.items()}
    return (species,aa_frequencies)

def main():
    if len(sys.argv) < 2:
        print("Usage: python calculate_AA_freq.py <fasta_file> or -d <dir>")
        sys.exit(1)
    
    parser = argparse.ArgumentParser(description="Calculate amino acid frequencies from a fasta file",
                                    epilog='Example: calculate_AA_freq.py <fasta_file> or -d <dir>')
    parser.add_argument("fasta_file", nargs='*', help="Input fasta file(s)")
    parser.add_argument("-d", "--dir", help="Input dir")
    parser.add_argument("-ext", "--ext", default="proteins.fa",help="file extension when reading folder")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')

    parser.add_argument("-o","--outfile", default='bigquery/aa_freq.csv', 
                        help="Output file")
    
    
    args = parser.parse_args()
    with open (args.outfile, "w",newline="") as outfh:
        outwriter = csv.writer(outfh)
        outwriter.writerow(["species_prefix","amino_acid","frequency"])
        n = 0
        if args.fasta_file:
            for file in args.fasta_file:
                if args.debug:
                    print(f"Processing {file}")
                timestart = time.time() 
                (prefix,frequencies) = calculate_aa_frequencies(file)
                timeend = time.time()
                print(f"Time elapsed: {timeend-timestart}")
                for aa, freq in sorted(frequencies.items()):
                    outwriter.writerow([prefix,aa,f"{freq:.4f}"])
                n += 1
        else:
            for file in os.listdir(args.dir):
                if file.endswith(args.ext):
                    if args.debug:
                        print(f"Processing {file}")
                    (prefix,frequencies) = calculate_aa_frequencies(os.path.join(args.dir,file))
                    for aa, freq in sorted(frequencies.items()):
                        outwriter.writerow([prefix,aa,f"{freq:.4f}"])
                    n += 1
                if n % 100 == 0:
                    print(f"Processed {n} files")
        print(f"Processed {n} protein files")
                        
if __name__ == "__main__":
    main()