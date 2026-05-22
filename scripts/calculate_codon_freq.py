#!/usr/bin/env python3

from Bio import SeqIO
from collections import Counter
import sys
import argparse
import csv
import os
import time


def validate(seq, alphabet="ACGT"):
    """Checks that a sequence only contains values from an alphabet"""
    alphabet = set(alphabet)
    leftover = set(seq.upper()) - alphabet
    return not leftover

def calculate_codon_frequencies(fasta_file):
    codon_counter = Counter()
    total_codons = 0
    species = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if not species:
            species = record.id.split("_")[0]
        seqstr = str(record.seq).upper()
        codon_count = int(len(seqstr) / 3) # round down to the nearest codon
        for i in range(0, 3 * codon_count,3):
            codon = seqstr[i:i+3]
            if validate(codon):
                codon_counter.update({codon: 1})
        total_codons += codon_count

    codon_frequencies = {codon: count / total_codons for codon, count in codon_counter.items()}
    return (species,codon_frequencies)

def main():
    if len(sys.argv) < 2:
        print("Usage: python calculate_codon_freq.py <fasta_file> or -d <dir>")
        sys.exit(1)
    
    parser = argparse.ArgumentParser(description="Calculate codon frequencies from a fasta file",
                                    epilog='Example: calculate_codon_freq.py <fasta_file> or -d <dir>')
    parser.add_argument("fasta_file", nargs='*', help="Input fasta file(s)")
    parser.add_argument("-d", "--dir", help="Input dir")
    parser.add_argument("-ext", "--ext", default="cds-transcripts.fa",help="file extension when reading folder")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')

    parser.add_argument("-o","--outfile", default='bigquery/codon_freq.csv', 
                        help="Output file")
    
    
    args = parser.parse_args()
    with open (args.outfile, "w",newline="") as outfh:
        outwriter = csv.writer(outfh)
        outwriter.writerow(["species_prefix","codon","frequency"])
        n = 0
        if args.fasta_file:
            for file in args.fasta_file:
                if args.debug:
                    print(f"Processing {file}")
                timestart = time.time() 
                (prefix,frequencies) = calculate_codon_frequencies(file)
                timeend = time.time()
                print(f"Time elapsed: {timeend-timestart}")
                for codon, freq in sorted(frequencies.items()):
                    outwriter.writerow([prefix,codon,f"{freq:.4f}"])
                n += 1
        else:
            for file in os.listdir(args.dir):
                if file.endswith(args.ext):
                    if args.debug:
                        print(f"Processing {file}")
                    (prefix,frequencies) = calculate_codon_frequencies(os.path.join(args.dir,file))
                    for codon, freq in sorted(frequencies.items()):
                        outwriter.writerow([prefix,codon,f"{freq:.4f}"])
                    n += 1
                if n % 100 == 0:
                    print(f"Processed {n} files")
        print(f"Processed {n} coding files")
                        
if __name__ == "__main__":
    main()
