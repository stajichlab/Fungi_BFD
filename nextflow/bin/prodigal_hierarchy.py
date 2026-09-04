#!/usr/bin/env python3
"""Re-emit a Prodigal CDS-only GFF3 as full gene/mRNA/CDS blocks.

EVM assembles gene models from gene/mRNA blocks; CDS-only features can only
act as support, never define structure (validated on Ordospora colligata
OC4: gene Sn 0.584 without this transform vs 0.849 for Prodigal alone --
see Microsporidia_predict PLAN.md 9.15). Each single-exon Prodigal CDS
becomes its own gene/mRNA/CDS block, blocks separated by a blank line
(EVM's lib.readBlocks in funannotate-runEVM.py expects every gene block
immediately preceded by one -- no ##gff-version header, which would join
the first block and crash gene_blocks_to_interlap).

Usage: prodigal_hierarchy.py RAW_PRODIGAL.gff3 OUT_HIER.gff3
"""
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    n = 0
    with open(src) as fh, open(dst, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue
            seqid, _, _, start, end, score, strand, phase, _ = cols
            n += 1
            gid, mid = f"prodigal_g{n}", f"prodigal_m{n}"
            out.write(f"{seqid}\tprodigal\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n")
            out.write(f"{seqid}\tprodigal\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={mid};Parent={gid}\n")
            out.write(f"{seqid}\tprodigal\tCDS\t{start}\t{end}\t{score}\t{strand}\t{phase}\tID={mid}.cds;Parent={mid}\n")
            out.write("\n")
    print(f"[INFO] prodigal_hierarchy: wrote {n} gene/mRNA/CDS blocks to {dst}")


if __name__ == "__main__":
    main()
