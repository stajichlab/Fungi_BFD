#!/usr/bin/env python3

import csv
import argparse
import os
import gzip

def process_hmmsearch_domtbl_file(fh, csvout):
    i = 0
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.strip().split()
        if len(fields) < 2:
            continue
        # Extract the relevant fields
        protein_id = fields[0]
        # protein_acc = fields[1]
        protein_len = int(fields[2])
        hmm_id = fields[3]
        hmm_acc = fields[4]
        hmm_len = int(fields[5])
        
        full_seq_e_value = fields[6]
        full_seq_score = fields[7]
        full_seq_bias = fields[8]
        domain_num = int(fields[9])
        domain_num_of = int(fields[10])
        domain_c_evalue = fields[11]
        domain_i_evalue = fields[12]
        domain_score = fields[13]
        domain_bias = fields[14]
        hmm_from = int(fields[15])
        hmm_to = int(fields[16])
        ali_from = int(fields[17])
        ali_to = int(fields[18])
        env_from = int(fields[19])
        env_to = int(fields[20])
        outrow = [ protein_id, hmm_id, hmm_acc, hmm_len,
                full_seq_e_value, full_seq_score, full_seq_bias,
                domain_num, domain_num_of, domain_c_evalue, domain_i_evalue,
                domain_score, domain_bias, hmm_from, hmm_to,
                ali_from, ali_to, env_from, env_to ]
        csvout.writerow(outrow)
        i += 1
    return i

def main():    
    parser = argparse.ArgumentParser(description="Convert Pfam TSV from Toronto to long form one line per domain",
                                    epilog='Example: pfamtsv_to_long.py')
    parser.add_argument("-i","--indir", help="Directory with pfam results pre-computed TSV files")
    parser.add_argument("tsv", nargs="*", help="Input TSV file(s)")
    parser.add_argument("-ext","--extension", default="domtblout", help="extension of pfam file [tsv]")
    parser.add_argument("-o","--outfile", default="bigquery/pfam.csv", 
                        help="output file [bigquery/pfam.csv]")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')
    
    args = parser.parse_args()
    
    with open(args.outfile, "w",newline="") as outfh:
        csvout = csv.writer(outfh)
        header = ['protein_id', 'pfam_id', 'pfam_acc', 'pfam_len',
                    'full_seq_e_value', 'full_seq_score', 'full_seq_bias',
                    'domain_num', 'domain_num_of', 'domain_c_evalue', 'domain_i_evalue',
                    'domain_score', 'domain_bias', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to',
                    'env_from','env_to' ]
        csvout.writerow(header)
        n = 0
        if args.indir:
            for file in os.listdir(args.indir):
                if args.debug:
                    print(f"Processing {file}")
                fh = None
                if file.endswith("." + args.extension):
                    filepath = os.path.join(args.indir, file)
                    fh = open(filepath, "r")
                elif file.endswith("." + args.extension + ".gz"):
                    filepath = os.path.join(args.indir, file)
                    fh = gzip.open(filepath, "rt")
                else:
                    continue
                linecount = process_hmmsearch_domtbl_file(fh, csvout)
                if linecount == 0:
                    print(f"Warning: {file} has no valid lines")
                    continue
                n += 1
                if args.debug and n > 1 and n % 1000 == 0:
                    print(f'Processed {n} files')
                fh.close()
        elif args.tsv:
            for file in args.tsv:
                fh = None
                if file.endswith("." + args.extension):
                    fh = open(file, "r")
                elif file.endswith("." + args.extension + ".gz"):
                    fh = gzip.open(file, "rt")
                else:
                    continue
                linecount = process_hmmsearch_domtbl_file(fh, csvout)
                if linecount == 0:
                    print(f"Warning: {file} has no valid lines")
                    continue
                fh.close()


if __name__ == "__main__":
    main()