#!/usr/bin/env python3

import argparse
import csv
import hashlib
import os
import sys
import re
import time
import tempfile
from contextlib import ExitStack

from Bio.Seq import Seq
from pyfaidx import Fasta
from pybedtools import BedTool

STOP_CODONS = ("TAA", "TAG", "TGA")

trna_regexp = re.compile(r'Predicted (\S+)')

def getGC(seq):
    """
    Calculate the GC content of a given DNA sequence.

    The GC content is the percentage of nucleotides in the sequence that are either
    guanine (G) or cytosine (C), including both uppercase and lowercase characters.

    Args:
        seq (str): A string representing the DNA sequence.

    Returns:
        float: The GC content as a fraction of the total sequence length.
        Returns 0 if the sequence length is 0.
    """
    gc = sum(str(seq).count(x) for x in "CGScgs")
    length = len(seq)
    if length == 0:
        return 0
    return 100 * (gc / length)


def parse_gff(gff, dna="", codon_table=1, debug=False):
    """
    Process a GFF file to extract gene statistics and optionally include exon/intron DNA sequences.
    Args:
        gff (str): Path to the GFF file to be parsed.
        dna (str, optional): Path to the DNA FASTA file. If provided, DNA sequences will be indexed and used to calculate GC content. Defaults to an empty string.
        debug (bool, optional): If True, enables debug mode which prints additional information and limits the number of processed genes. Defaults to False.
    Returns:
        dict: A dictionary containing gene data, where keys are gene IDs and values are dictionaries with gene statistics and transcript information.
    Raises:
        ValueError: If the GFF file contains invalid or unexpected data.
    Notes:
        - The function processes different feature types such as 'gene', 'mRNA', 'tRNA', 'exon', and 'CDS'.
        - For 'exon' and 'CDS' features, the GC content is calculated if the DNA file is provided.
        - The function prints progress information if debug mode is enabled.
    """
    if debug:
        print(f"DEBUG: parse_gff - Reading {gff}, dna is {dna}")
    genedata = {}
    dnadb = None
    if dna:
        # consider index_db and compressed bgz fasta for speed/space?
        # dnadb = SeqIO.index_db(dna + ".idx",dna,format='fasta')
        # dnadb = SeqIO.index(dna,format='fasta')
        dnadb = Fasta(dna)
    tRNA_gff = os.path.join(os.path.dirname(os.path.realpath(gff)), 
                            "../predict_misc/trnascan.no-overlaps.gff3")
    temp_tRNA = None
    if os.path.exists(tRNA_gff):
        tRNA_gff = os.path.realpath(tRNA_gff)        
        temp_tRNA = tempfile.NamedTemporaryFile(delete=False,suffix=".bed")
    else:
        tRNA_gff = None

    with open(gff, "r") as gff_fh:
        if debug:
            print(f"DEBUG: tRNA_gff is {tRNA_gff}")
            print(f"DEBUG: tRNA_temp is {temp_tRNA.name}")

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
            if debug and genedata and len(genedata) % 1000 == 0 and fields[2] == "gene":
                time1 = time.time()
                print(
                    f"DEBUG: Processed {len(genedata)} genes in {gff} in {time1-time0} seconds"
                )
                time0 = time1

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
                    "transcripts": {},
                }
            elif ftype == "mRNA":
                if not ("ID" in group_data and "Parent" in group_data):
                    print(
                        f"WARNING: Cannot parse groups {group_data} not ID or Parent in {gff}\n{line}"
                    )
                    continue
                mrna_id = group_data["ID"]
                gene_id = group_data["Parent"]
                transcript2gene[mrna_id] = gene_id
                if gene_id not in genedata:
                    print(f"WARNING: mRNA {mrna_id} has no gene in {gff}")
                    continue
                genedata[gene_id]["type"] = "protein_coding"
                genedata[gene_id]["transcripts"][mrna_id] = {
                    "chrom": fields[0],
                    "start": fstart,
                    "end": fend,
                    "strand": fstrand,
                    "is_partial": "NULL",
                    "has_start_codon": "NULL",
                    "has_stop_codon": "NULL",
                    "exon": [],
                    "CDS": [],
                    "intron": [],
                    "CDS_seq": {
                        "id": mrna_id,
                        "length": None,
                        "md5checksum": None,
                        "cdsseq": None,
                    },
                    "protein": {
                        "id": f"{mrna_id}",
                        "parent": mrna_id,
                        "length": None,
                        "md5checksum": None,
                    },
                }
            elif ftype == "tRNA":
                if not ("ID" in group_data and "Parent" in group_data):
                    print(
                        f"WARNING: Cannot parse groups {group_data} no ID or Parent in {gff}\n{line}"
                    )
                    continue
                trna_id = group_data["ID"]
                gene_id = group_data["Parent"]
                amino_acid = group_data.get("product", "NULL")
                amino_acid = amino_acid.replace("tRNA-", "")
                amino_acid = amino_acid.replace(";", "")
                if tRNA_gff:
                    row = "\t".join([fields[0], str(fstart), str(fend), gene_id])+"\n"
                    temp_tRNA.write(b'%s' % row.encode())
                transcript2gene[trna_id] = gene_id
                if gene_id not in genedata:
                    print(f"WARNING: tRNA {trna_id} has no gene in {gff}")
                    continue
                genedata[gene_id]["type"] = "tRNA_gene"                
                genedata[gene_id]["tRNA_amino_acid"] = amino_acid
                genedata[gene_id]["codon"] = ""
                genedata[gene_id]["anticodon"] = ""
                genedata[gene_id]["transcripts"][trna_id] = {
                    "chrom": fields[0],
                    "start": fstart,
                    "end": fend,
                    "strand": fstrand,
                    "is_partial": "FALSE",
                    "has_start_codon": "NULL",
                    "has_stop_codon": "NULL",
                    "exon": [],
                    "intron": [],
                    "CDS_seq": {
                            'length': None,
                            'cdsseq': None,
                            'md5checksum': None,
                        }
                    }
            elif ftype in ("exon", "CDS"):
                if "Parent" not in group_data:
                    print(
                        f"WARNING: Group data {group_data} no Parent in {gff}\n{line}"
                    )
                    continue
                parent_id = group_data["Parent"]
                gene_id = None
                if parent_id not in transcript2gene:
                    print(
                        f"WARNING: Exon from transcript {parent_id} cannot map gene id in {gff}\n{line}"
                    )
                    continue
                else:
                    gene_id = transcript2gene[parent_id]

                if (
                    gene_id not in genedata
                    or parent_id not in genedata[gene_id]["transcripts"]
                ):
                    print(
                        f"WARNING: Exon of {parent_id} has no gene or mRNA in {gff}\n{line}"
                    )
                    continue
                n = len(genedata[gene_id]["transcripts"][parent_id][ftype]) + 1
                exon_id = f"{parent_id}.{ftype}{n}"
                # if "ID" in group_data:  # override with existing value if provided
                #   exon_id = group_data["ID"]
                # zero base indexing
                exonseq_GC = getGC(dnadb[fields[0]][fstart - 1 : fend])
                genedata[gene_id]["transcripts"][parent_id][ftype].append(
                    {
                        "id": exon_id,
                        "chrom": fields[0],
                        "start": fstart,
                        "end": fend,
                        "strand": fstrand,
                        "GC_content": f"{exonseq_GC:0.2f}",
                        "order": None,
                        "frame": fields[7],
                    }
                )
    for gene_name, gene in genedata.items():
        chrom_segment = dnadb[gene["chrom"]]
        for transcript_name, transcript in gene["transcripts"].items():
            if debug:
                print(f"DEBUG: Processing {transcript_name}")
            exonlist = sorted(
                transcript["exon"], key=lambda x: x["strand"] * x["start"]
            )
            e = 0
            lastexon = {}
            n = 0
            for exon in exonlist:
                # give an order for the exons based on 5' to 3' direction
                exon["order"] = e
                if debug:
                    print(f"DEBUG: Exon {exon}")
                e += 1
                # extract introns from the exons
                if lastexon:
                    # reverse complement introns require diff
                    # start/end compare
                    if exon["strand"] == -1:
                        intronstart = exon["end"] + 1
                        intronend = lastexon["start"] - 1
                    else:
                        intronstart = lastexon["end"] + 1
                        intronend = exon["start"] - 1
                    # zero based indexing
                    if intronstart > intronend:
                        print(
                            f"WARNING: improper start/end for exon {exon} lastexon: {lastexon}"
                        )
                        return
                    intron = chrom_segment[intronstart - 1 : intronend]
                    if exon["strand"] == -1:
                        intron = -intron
                    transcript["intron"].append(
                        {
                            "id": f"{transcript_name}.intron{n}",
                            "parent_id": transcript_name,
                            "intron_number": n,
                            "chrom": exon["chrom"],
                            "start": intronstart,
                            "end": intronend,
                            "strand": exon["strand"],
                            "GC_content": f"{getGC(intron):0.2f}",
                            "seq": intron,
                            "splice_5": intron[:2],
                            "splice_3": intron[-2:],
                            "codon_position": None,
                            "codon_frame": None,
                        }
                    )
                    n += 1
                lastexon = exon

            # give the CDS features proper order
            c = 0
            if "CDS" in transcript:
                transcript["CDS"] = sorted(
                    transcript["CDS"], key=lambda x: x["strand"] * x["start"]
                )
                CDS_sequence = ""
                intron_index = None
                lastcds = None
                for cds in transcript["CDS"]:
                    cds["order"] = c
                    frame = int(cds["frame"])
                    if debug:
                        print(f"DEBUG: CDS {cds}")
                    (cds_start, cds_end) = (cds["start"] - 1, cds["end"])
                    if c == 0 and frame > 0:
                        if cds["strand"] == -1:
                            cds_end -= frame
                        else:
                            cds_start += frame
                    CDS_exon_seq = chrom_segment[cds_start:cds_end]
                    if cds["strand"] == -1:
                        CDS_exon_seq = -CDS_exon_seq
                    CDS_sequence += CDS_exon_seq.seq

                    if c == 0:
                        if CDS_exon_seq[:3] == "ATG":
                            transcript["has_start_codon"] = "TRUE"
                        else:
                            if debug:
                                print(
                                    f"DEBUG: start codon for {transcript_name} is {CDS_exon_seq[:3]} strand={cds['strand']} {cds_start}:{cds_end}"
                                )
                            transcript["has_start_codon"] = "FALSE"
                            transcript["is_partial"] = "TRUE"
                    else:
                        if debug:
                            print(
                                f"WARNING: lastcds is {lastcds} intron_index {intron_index}"
                            )
                        # need to sync the intron index with the CDS
                        if intron_index is None:
                            for i, intron in enumerate(transcript["intron"]):
                                # need to check if this intron falls between two CDS
                                if cds["strand"] == -1:
                                    if (
                                        intron["start"] < lastcds["start"]
                                        and intron["start"] > cds["end"]
                                    ):
                                        intron_index = i
                                        break
                                else:
                                    if (
                                        intron["start"] > lastcds["end"]
                                        and intron["start"] < cds["start"]
                                    ):
                                        intron_index = i
                                        break
                            if intron_index is None:
                                print(
                                    f"WARNING: Could not find intron between: lastcds is {lastcds} and cds is {cds}"
                                )
                                print(f'WARNING: introns are {transcript["intron"]}')
                        else:
                            intron_index += 1

                        intronobj = transcript["intron"][intron_index]
                        intronobj["codon_position"] = int(len(CDS_sequence) / 3)
                        intronobj["codon_frame"] = len(CDS_sequence) % 3
                    lastcds = cds
                    c += 1

                if CDS_sequence[-3:] in STOP_CODONS:
                    transcript["has_stop_codon"] = "TRUE"
                else:
                    if debug:
                        print(
                            f"DEBUG: no stop codon for {transcript_name} is {CDS_sequence[-3:]}"
                        )
                        print(f"DEBUG: CDS_seq is {CDS_sequence}")
                    transcript["has_stop_codon"] = "FALSE"
                    transcript["is_partial"] = "TRUE"
                if (
                    transcript["has_stop_codon"] == "TRUE"
                    and transcript["has_start_codon"] == "TRUE"
                ):
                    transcript["is_partial"] = "FALSE"
                # if we have codon table lookup we can use that

                proteinseq = Seq.translate(Seq(str(CDS_sequence)), table=codon_table)
                if proteinseq[-1:] == "*":
                    proteinseq = proteinseq[:-1]  # strip trailing stop codon.
                transcript["protein"]["length"] = len(str(proteinseq))
                transcript["protein"]["pepseq"] = str(proteinseq)
                transcript["protein"]["md5checksum"] = hashlib.md5(
                    str(proteinseq).encode()
                ).hexdigest()
                transcript["CDS_seq"]["cdsseq"] = CDS_sequence
                transcript["CDS_seq"]["length"] = len(CDS_sequence)
                transcript["CDS_seq"]["md5checksum"] = hashlib.md5(
                    str(CDS_sequence).encode()
                ).hexdigest()

    if tRNA_gff:
        temp_tRNA.close()
        tRNA_bed = BedTool(temp_tRNA.name)
        tRNA_gff_bedtools = BedTool(tRNA_gff)
        tRNAs = tRNA_bed.intersect(tRNA_gff_bedtools,wo=True)
        for tRNA in tRNAs:
            chrom = tRNA.chrom
            start = tRNA.start
            end = tRNA.end
            trna_id = tRNA.name            
            if trna_id not in genedata:
                print(f"WARNING: tRNA {trna_id} has no gene in {gff}")
                continue
            gtype = tRNA[6]
            if gtype != "tRNA":
                continue
            for f in tRNA[12].split(";"):
                (tag, note) = f.split("=")
                if tag == "note":
                    m = trna_regexp.search(note)
                    if m:
                        anticodon = m.group(1)
                        codon = Seq(anticodon).reverse_complement()
                        genedata[trna_id]["tRNA_codon"] = str(Seq(anticodon).reverse_complement())
                        genedata[trna_id]["tRNA_anticodon"] = anticodon
                        if debug:
                            print(f"DEBUG: found and adding tRNA_codon for {trna_id} with {anticodon} -> {codon}")
                    else:
                        print(f"WARNING: No anticodon found in {note} for {trna_id}")            
        if os.path.exists(temp_tRNA.name):
            os.unlink(temp_tRNA.name)
    return genedata


def main():
    if len(sys.argv) < 2:
        print("Usage: python gff_build_big_query.py gff or -d <dir>")
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description="Calculate gene stats from GFF file(s)",
        epilog="Example: gff_build_big_query.py gff_file.gff -d dna_dir -p pep_dir",
    )
    parser.add_argument("gff_file", nargs="*", help="Input GFF file(s)")
    parser.add_argument("-g", "--gff_dir", help="GFF files dir")
    parser.add_argument("-d", "--dna_dir", default="genomes", help="DNA files dir")
    parser.add_argument("-p", "--pep_dir", default="input", help="Protein files dir")
    parser.add_argument(
        "-gffext",
        "--gffext",
        default="gff3",
        help="file extension when reading gff folder",
    )
    parser.add_argument(
        "-pepext",
        "--pepext",
        default="proteins.fa",
        help="file extension when reading proteins folder",
    )
    parser.add_argument(
        "-dnaext",
        "--dnaext",
        default="scaffolds.fa",
        help="file extension when reading proteins folder",
    )

    parser.add_argument("-v", "--debug", help="Debugging output", action="store_true")
    parser.add_argument(
        "-o",
        "--outdir",
        default="bigquery",
        help="Output folder for gene info, exons, introns, transcripts, proteins",
    )
    args = parser.parse_args()
    if args.debug:
        print(args)
    if args.debug:
        print(
            f"Reading GFF files from {args.gff_dir} and DNA files from '{args.dna_dir}'"
        )
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    output_files = [
        "gene_info.csv",
        "gene_exons.csv",
        "gene_CDS.csv",
        "gene_introns.csv",
        "gene_transcripts.csv",
        "gene_trnas.csv",
        "gene_proteins.csv",
    ]
    with ExitStack() as stack:
        files = [
            stack.enter_context(open(f"{args.outdir}/{filename}", "w", newline=""))
            for filename in output_files
        ]
        genefile, exonfile, CDSfile, intronsfile, mrnafile, trnafile, pepfile = files
        genecsv = csv.writer(genefile)
        exoncsv = csv.writer(exonfile)
        CDScsv = csv.writer(CDSfile)
        introncsv = csv.writer(intronsfile)
        mrnacsv = csv.writer(mrnafile)
        trnacsv = csv.writer(trnafile)
        pepcsv = csv.writer(pepfile)
        genecsv.writerow(
            [
                "gene_id",
                "LOCUSTAG",
                "chrom",
                "start",
                "end",
                "strand",
                "gene_type",
            ]
        )
        mrnacsv.writerow(
            [
                "gene_id",
                "transcript_id",
                "chrom",
                "start",
                "end",
                "strand",
                "is_partial",
                "has_start_codon",
                "has_stop_codon",
                "CDS_sequence",
                "CDS_length",
                "md5checksum"
            ]
        )
        trnacsv.writerow(
            [
                "gene_id",
                "amino_acid",
                "codon"
            ]
        )
        exoncsv.writerow(
            [
                "exon_id",
                "transcript_id",
                "order",
                "chrom",
                "start",
                "end",
                "strand",
                "GC_content",
            ]
        )
        CDScsv.writerow(
            [
                "cds_id",
                "transcript_id",
                "order",
                "chrom",
                "start",
                "end",
                "strand",
                "GC_content",
            ]
        )
        introncsv.writerow(
            [
                "intron_id",
                "transcript_id",
                "intron_number",
                "chrom",
                "start",
                "end",
                "strand",
                "splice_5",
                "splice_3",
                "GC_content",
                "codon_position",
                "codon_frame",
                "seq",
            ]
        )
        pepcsv.writerow(
            [
                "gene_id",
                "protein_id",
                "transcript_id",
                "length",
                "peptide",
                "md5checksum",
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
            dnafile = ""
            if args.dna_dir and os.path.isdir(args.dna_dir):
                if args.debug:
                    print(f"DEBUG: Adding for {fstem} in {args.dna_dir}")
                dnafile = os.path.join(args.dna_dir, f"{fstem}.{args.dnaext}")
            genedata = parse_gff(gff=gff_file_l, dna=dnafile, debug=args.debug)
            if genedata:
                IDs_to_process.add(fstem)
            species = None
            for genename, gene in genedata.items():
                if not species:
                    (species) = genename.split("_")[0]
                genecsv.writerow(
                    [
                        genename,
                        species,
                        gene["chrom"],
                        gene["start"],
                        gene["end"],
                        gene["strand"],
                        gene["type"]
                    ]
                )
                if gene["type"] == "tRNA_gene":
                    trnacsv.writerow(
                        [
                            genename,
                            gene["tRNA_amino_acid"],
                            gene["tRNA_codon"]
                        ]
                    )
                # consider saving space by only encoding strand on the gene level
                for transcriptname, transcript in gene["transcripts"].items():
                    mrnacsv.writerow(
                        [
                            genename,
                            transcriptname,
                            transcript["chrom"],
                            transcript["start"],
                            transcript["end"],
                            transcript["strand"],
                            transcript["is_partial"],
                            transcript["has_start_codon"],
                            transcript["has_stop_codon"],
                            transcript["CDS_seq"]["cdsseq"],
                            transcript["CDS_seq"]["length"],
                            transcript["CDS_seq"]["md5checksum"]
                        ]
                    )
                    if "exon" in transcript:
                        for exon in transcript["exon"]:
                            exoncsv.writerow(
                                [
                                    exon["id"],
                                    transcriptname,
                                    exon["order"],
                                    exon["chrom"],
                                    exon["start"],
                                    exon["end"],
                                    exon["strand"],
                                    exon["GC_content"],
                                ]
                            )

                    if "CDS" in transcript:
                        for cds in transcript["CDS"]:
                            CDScsv.writerow(
                                [
                                    cds["id"],
                                    transcriptname,
                                    cds["order"],
                                    cds["chrom"],
                                    cds["start"],
                                    cds["end"],
                                    cds["strand"],
                                    cds["GC_content"],
                                ]
                            )
                    if "intron" in transcript:
                        for intron in transcript["intron"]:
                            introncsv.writerow(
                                [
                                    intron["id"],
                                    intron["parent_id"],
                                    intron["intron_number"],
                                    intron["chrom"],
                                    intron["start"],
                                    intron["end"],
                                    intron["strand"],
                                    intron["splice_5"],
                                    intron["splice_3"],
                                    intron["GC_content"],
                                    intron["codon_position"],
                                    intron["codon_frame"],
                                    intron["seq"],
                                ]
                            )
                    if "protein" in transcript:
                        pepcsv.writerow(
                            [
                                genename,
                                transcript["protein"]["id"],
                                transcriptname,
                                transcript["protein"]["length"],
                                transcript["protein"]["pepseq"],
                                transcript["protein"]["md5checksum"],
                            ]
                        )


if __name__ == "__main__":
    main()
