#!/usr/bin/env python
"""Review finding 1 (DECISIONS D41): rank order inside getBestModel loci.

Transcript-span, same-strand loci (the current tx_strand rule), with three
rankings of the models in a locus:
  complete_first -- (complete ORF, CDS exons, CDS length, TPM, id)   [current code]
  struct_first   -- (CDS exons, CDS length, complete ORF, TPM, id)
  guarded        -- complete first only among models whose CDS length is
                    >= 80% of the longest CDS in the locus; otherwise structure
Scored against RefSeq with benchmark.score + diversity, for every PASA source
that exists: results/ (rc.1 baseline), r1_fix/ (R1), r1_r2/ (R1+R2).
Output: rank_benchmark.tsv next to this script (or argv[1]).
"""
import csv, gzip, os, shutil, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import benchmark as B
import diversity
from funannotate import library as lib

X = f"{B.R}/do_pasa_rust_vs_perl"
SOURCES = {"rc1": f"{X}/results", "R1": f"{X}/r1_fix", "R1R2": f"{X}/r1_r2"}


def expression(tpm):
    e = {}
    with open(tpm) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t, g, loc, v = line.rstrip("\n").split("\t")
            e[g] = max(float(v), e.get(g, 0.0))
    return e


def best(genes, expr, order):
    def feats(gid):
        g = genes[gid]
        bestf = (False, 0, 0)
        for i, c in enumerate(g.get("CDS") or []):
            f = (lib.is_complete_model(g, i), len(c), sum(e[1] - e[0] + 1 for e in c))
            bestf = max(bestf, f)
        return bestf  # (complete, nCDS, cdsLen)

    models = [(k, v["contig"], v["strand"], v["location"][0], v["location"][1]) for k, v in genes.items()]
    keep = set()
    for c in lib.cluster_overlapping(models, B.OVERLAP / 100.0, strand_aware=True):
        F = {k: feats(k) for k in c}
        if order == "complete_first":
            key = lambda k: (F[k][0], F[k][1], F[k][2], expr.get(k, 0.0), k)
        elif order == "struct_first":
            key = lambda k: (F[k][1], F[k][2], F[k][0], expr.get(k, 0.0), k)
        else:  # guarded
            longest = max(F[k][2] for k in c)
            key = lambda k: (F[k][0] and F[k][2] >= 0.8 * longest, F[k][1], F[k][2], expr.get(k, 0.0), k)
        keep.add(max(c, key=key))
    return {k: genes[k] for k in keep}


def main(out_tsv):
    rows = []
    for src, root in SOURCES.items():
        for name, (asm, acc) in B.CASES.items():
            d = f"{root}/{name}.rust"
            if not (os.path.exists(f"{d}/pasa.step1.gff3") and os.path.exists(f"{d}/kallisto.tsv")):
                continue
            ref = B.load_refseq(f"{X}/refseq/{acc}.gff")
            tmp = tempfile.mkdtemp(dir=os.environ.get("SCRATCH", "/tmp"))
            g = f"{tmp}/g.fa"
            with gzip.open(f"{B.R}/input_clean_genomes/{asm}.masked.fasta.gz", "rb") as i, open(g, "wb") as o:
                shutil.copyfileobj(i, o)
            _, genes = lib.gff2interlap(f"{d}/pasa.step1.gff3", g)
            expr = expression(f"{d}/kallisto.tsv")
            for order in ("complete_first", "struct_first", "guarded"):
                sel = best(genes, expr, order)
                r = {"source": src, "genome": name, "rank": order}
                r.update(B.score(sel, ref))
                r.update(diversity.diversity(sel, ref))
                rows.append(r)
                print("\t".join(str(v) for v in list(r.values())[:12]), flush=True)
            shutil.rmtree(tmp, ignore_errors=True)
    with open(out_tsv, "w") as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "rank_benchmark.tsv"))
