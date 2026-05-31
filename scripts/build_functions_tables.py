#!/usr/bin/env python3

import argparse
import os
import csv
import gzip
import re
# This script prepares data for loading into BigQuery.
# It reads various function prediction results from a specified directory,
# processes them, and writes the results into CSV files in a specified output directory.
# The output files are named based on the input directory names and are stored in a directory named 'bigquery'.
# The script includes functions for processing MEROPS, CAZY, SignalP, TMHMM, WolfPSORT, and TargetP data.
# The script checks if the output files already exist and skips processing if they do, unless forced to overwrite.
# The output files are compressed with gzip if they already exist.

outdir = 'tables'
def merops(indir="results/function/merops",force=False):
    """Collect MEROPS BLAST tabular results from indir into a gzipped CSV."""
    # load MEROPS data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        writer.writerow(['species_prefix','protein_id','merops_id','percent_identity','aln_length','mismatches','gap_openings','q_start','q_end',
                        's_start','s_end', 'evalue', 'bitscore'])
        for file in os.listdir(indir):
            if file.endswith(".blasttab.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    for row in reader:
                        prefix = row[0].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def cazy_overview(indir="results/function/cazy",force=False):
    """Collect run_dbcan overview TSV files from indir into a gzipped CSV."""
    # load CAZY data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".overview.csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        # Gene_ID	EC	cazyme_fam	sub_fam	diamond_fam	Substrate	#ofTools
        writer.writerow(['species_prefix','protein_id','EC','cazyme_fam','sub_fam','diamond_fam','substrate','toolcount'])
        for spdir in os.listdir(indir):
            infile = os.path.join(indir, spdir, spdir + '.overview.tsv.gz')
            if os.path.exists(infile):
                with gzip.open(infile, "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    next(reader)
                    for row in reader:
                        prefix = row[0].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def cazy_hmm(indir="results/function/cazy",force=False):
    """Collect run_dbcan HMM cazymes TSV files from indir into a gzipped CSV."""
    # load CAZY data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".cazymes_hmm.csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        # HMM_Profile	Profile_Length	Gene_ID	Gene_Length	Evalue	Profile_Start	Profile_End	Gene_Start	Gene_End	Coverage
        writer.writerow(['species_prefix','HMM_id','profile_length','protein_id','protein_length','evalue',
                        'q_start','q_end','s_start','s_end', 'coverage'])
        for spdir in os.listdir(indir):
            infile = os.path.join(indir, spdir, spdir + '.cazymes.tsv.gz')
            if os.path.exists(infile):
                with gzip.open(infile, "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    next(reader)
                    for row in reader:
                        prefix = row[2].split('_')[0]
                        newrow = [prefix]
                        newrow.extend(row)
                        writer.writerow(newrow)

def signalp(indir="results/function/signalp",force=False):
    """Collect SignalP GFF3 output files from indir into a gzipped CSV of signal peptide predictions."""
    # load signalp data
    outfile = os.path.join(outdir,os.path.basename(indir) + ".signal_peptide.csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        # HMM_Profile	Profile_Length	Gene_ID	Gene_Length	Evalue	Profile_Start	Profile_End	Gene_Start	Gene_End	Coverage
        writer.writerow(['species_prefix','protein_id','peptide_start','peptide_end','probability'])
        for file in os.listdir(indir):
            if file.endswith(".signalp.gff3.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    for row in reader:
                        if row[0].startswith('#'):
                            continue
                        id = row[0].split(' ')[0]
                        prefix = id.split('_')[0]                        
                        newrow = [prefix, id, row[3], row[4], row[5]]
                        writer.writerow(newrow)

def tmhmm(indir="results/function/tmhmm",force=False):
    """Collect TMHMM short-format result files from indir into a gzipped CSV; skips proteins with 0 predicted helices."""
    outfile = os.path.join(outdir,os.path.basename(indir) + ".csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        writer.writerow(['species_prefix','protein_id','len','ExpAA','First60','PredHel','Topology'])
        for file in os.listdir(indir):
            if file.endswith(".tmhmm_short.tsv.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:
                    reader = csv.reader(infh, delimiter='\t')
                    for row in reader:
                        if row[0].startswith('#'):
                            continue
                        id = row[0]
                        prefix = id.split('_')[0]
                        newrow = [prefix,id]
                        for n in row[1:]:
                            (key,value) = n.split('=')
                            newrow.append(value)
                        if newrow[-2] == '0':
                            continue

                        writer.writerow(newrow)

def wolfpsort(indir="results/function/wolfpsort",force=False,onlybest=True):
    """Collect WoLF PSORT result files from indir into a gzipped CSV; when onlybest=True, only the top-ranked localization per protein is kept."""
    outfile = os.path.join(outdir,os.path.basename(indir) + ".csv.gz")
    if os.path.exists(outfile) and not force:
        return
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        writer.writerow(['species_prefix','protein_id','localization','score'])
        for file in os.listdir(indir):
            if file.endswith(".wolfpsort.results.txt.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:                    
                    for line in infh:
                        if line.startswith('#'):
                            continue
                        (id,resultstr) = line.strip().split(' ',1)
                        results = resultstr.split(', ')
                        prefix = id.split('_')[0]
                        for scoring in results:
                            (code,score) = scoring.split(' ')                        
                            newrow = [prefix,id,code,score]
                            writer.writerow(newrow)
                            if onlybest:
                                break


# redicts the presence of N-terminal presequences: signal peptide (SP),
# mitochondrial transit peptide (mTP), chloroplast transit peptide (cTP)
# or thylakoid luminal transit peptide (lTP). For the sequences predicted to 
# contain an N-terminal presequence a potential cleavage site is also predicted.

# The type can be

# "SP" for signal peptide,
# "MT" for mitochondrial transit peptide (mTP),
# "CH" for chloroplast transit peptide (cTP),
# "TH" for thylakoidal lumen composite transit peptide (lTP),
# "Other" for no targeting peptide (in this case, the length is given as 0).

# s the position where the sorting signal is cleaved.
# This is encoded as a zero vector of length 200 with 1 in the cleavage site position.

def targetp(indir="results/function/targetP",force=False):
    """Collect TargetP2 summary files from indir into a gzipped CSV; skips proteins predicted as noTP."""
    outfile = os.path.join(outdir,os.path.basename(indir) + ".csv.gz")
    if os.path.exists(outfile) and not force:
        return
    CSmatch = re.compile(r'CS pos:\s+(\d+)-(\d+)\.\s+(\S+)\.\s+Pr:\s+(\S+)')
    with gzip.open(outfile, "wt", newline='') as of:
        writer = csv.writer(of)
        writer.writerow(['species_prefix','protein_id','prediction','probability',
                        'cleavage_position_start', 'cleavage_position_end',
                        'cleavage_probability', 'motif'])
        for file in os.listdir(indir):
            if file.endswith("_summary.targetp2.gz"):
                with gzip.open(os.path.join(indir,file), "rt") as infh:                    
                    for line in infh:
                        if line.startswith('#'):
                            continue
                        (id,prediction,results) = line.strip().split('\t',2)
                        if prediction == "noTP":
                            continue
                        (noTP,SP,mTP,CS) = results.split('\t')
                        prefix = id.split('_')[0]
                        probability = None
                        if prediction == "SP":
                            probability = SP
                        elif prediction == "mTP":
                            probability = mTP

                        m = CSmatch.match(CS)
                        cleavage_position_start = '0'
                        cleavage_position_end = '0'
                        motif = ''
                        cleavage_probability = '0.0'
                        if m:
                            (cleavage_position_start,cleavage_position_end,
                            motif,cleavage_probability) = m.groups()
                        
                        datarow = [prefix,id,prediction, probability,
                                cleavage_position_start, cleavage_position_end,
                                cleavage_probability, motif]
                        
                        writer.writerow(datarow)

def kegg(indir="results/function/kegg",force=False):
    """Placeholder — KEGG annotation collection not yet implemented."""
    print('not running yet')

def busco(indir="results/stats/busco",force=False):
    """Placeholder — BUSCO result collection not yet implemented."""
    print('not running yet')


# no pfam as we run this from Toronto dataset instead of local generation

parser = argparse.ArgumentParser(description="Build function summary tables for BigQuery loading.")
parser.add_argument('--force', action='store_true', help='Overwrite existing output files')
args = parser.parse_args()

merops(force=args.force)
cazy_overview(force=args.force)
cazy_hmm(force=args.force)
signalp(force=args.force)
#kegg()
tmhmm(force=args.force)
wolfpsort(force=args.force)
#busco()
targetp(force=args.force)
