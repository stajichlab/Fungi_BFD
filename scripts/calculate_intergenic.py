#!/usr/bin/env python3

import argparse
import csv
import hashlib
import os
import sys
import time
from contextlib import ExitStack

def parse_gff_get_distances(gff, debug=False):
    """
    Process a GFF file to extract gene distances
    Args:
        gff (str): Path to the GFF file to be parsed.
        debug (bool, optional): If True, enables debug mode which prints additional information and limits the number of processed genes. Defaults to False.
    Returns:
        array: An array of intergenic distances
    Raises:
        ValueError: If the GFF file contains invalid or unexpected data.
    Notes:
        - The function processes only 'gene' features
        - The function prints progress information if debug mode is enabled.
    """
    if debug:
        print(f"DEBUG: parse_gff_get_distances - Reading {gff}")
    genedata = {}
    
    with open(gff, "r") as gff_fh:
        transcript2gene = {}
        time0 = time.time()
        for line in gff_fh:
            timestart = time.time()
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                if debug:
                    print(f"DEBUG: Skipping line with {len(fields)} fields in {gff}")
                continue

            group_data = {}
            for f in fields[8].split(";"):
                if not f or "=" not in f:
                    continue
                (tag, value) = f.split("=")
                group_data[tag] = value
            (fstart, fend) = sorted([int(fields[3]), int(fields[4])])
            fstrand = -1 if fields[6] == "-" else 1

            ftype = fields[2]
            
            if ftype == "gene":
                if "ID" not in group_data:
                    print(
                        f"WARNING: Cannot parse groups {group_data} does not start with ID in {gff}\n{line}"
                    )
                    continue
                gene_id = group_data["ID"]
                if gene_id in genedata:
                    print(f"WARNING: Duplicate gene ID {gene_id} in {gff}")
                    continue
                genedata[gene_id] = {
                    "chrom": fields[0],
                    "start": fstart,
                    "end": fend,
                    "strand": fstrand,
                    "type": "NULL",
                }
            else:
                continue

    if debug and genedata:
        time1 = time.time()
        print(
            f"DEBUG: Processed {len(genedata)} genes in {gff} in {time1-time0} seconds"
            )
        
    last_gene = []
    last_gene_name = ""
    distances = []
    for gene_name in sorted(genedata.keys(), key=lambda x: (genedata[x]["chrom"],genedata[x]["start"])):
        gene = genedata[gene_name]
        if last_gene and last_gene["chrom"] == gene["chrom"]:
            distance = gene["start"] - last_gene["end"]
            distances.append([last_gene_name, gene_name, distance])
        last_gene = gene
        last_gene_name = gene_name

    return distances


def main():
    if len(sys.argv) < 2:
        print("Usage: python calculate_intergenic.py gff or -d <dir>")
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description="Calculate intergenic distances from GFF file(s)",
        epilog="Example: calculate_intergenic.py gff_file.gff",
    )
    parser.add_argument("gff_file", nargs="*", help="Input GFF file(s)")
    parser.add_argument("-g", "--gff_dir", help="GFF files dir")
    parser.add_argument(
        "-gffext",
        "--gffext",
        default="gff3",
        help="file extension when reading gff folder",
    )

    parser.add_argument("-v", "--debug", help="Debugging output", action="store_true")
    parser.add_argument(
        "-o",
        "--outdir",
        default="bigquery",
        help="Output folder for distances",
    )

    args = parser.parse_args()
    if args.debug:
        print(args)

    if args.debug:
        print(
            f"Reading GFF files from {args.gff_dir}"
        )
    output_files = [
        "gene_pairwise_distances.csv",
    ]
    with ExitStack() as stack:
        files = [
            stack.enter_context(open(f"{args.outdir}/{filename}", "w", newline=""))
            for filename in output_files
        ]
        genefile = files[0]
        genecsv = csv.writer(genefile)
        genecsv.writerow(
            [
                "species_prefix",
                "left_gene",
                "right_gene",
                "distance",
            ]
        )

        IDs_to_process = set()
        filenames = []
        if args.gff_dir and os.path.isdir(args.gff_dir):
            for gff_file in os.listdir(args.gff_dir):
                if gff_file.endswith(args.gffext):
                    filename_stem = gff_file.replace(args.gffext, "")
                    if filename_stem.endswith("."):
                        filename_stem = filename_stem[:-1]
                    filenames.append(
                        [filename_stem, os.path.join(args.gff_dir, gff_file)]
                    )
        elif args.gff_file:
            for gff_file in args.gff_file:
                filename_stem = os.path.basename(gff_file).replace(args.gffext, "")
                if filename_stem.endswith("."):
                    filename_stem = filename_stem[:-1]
                filenames.append([filename_stem, gff_file])
        for gff_tuple in filenames:
            (fstem, gff_file_l) = gff_tuple
            genedata = parse_gff_get_distances(gff=gff_file_l, debug=args.debug)
            if genedata:
                IDs_to_process.add(fstem)
            species = None
            for d in genedata:
                if not species:
                    (species) = d[0].split("_")[0]
                row = [species] + d
                genecsv.writerow(row)

if __name__ == "__main__":
    main()