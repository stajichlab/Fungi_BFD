#!/usr/bin/env python
"""Benchmark PASA model selection against RefSeq CDS structures.

For each genome with an rc.1 PASA run (pasa.step1.gff3 + kallisto.tsv), build:
  getBestModel output (PASA evidence passed to EVM) with four locus rules
    old        -- committed funannotate (git HEAD, see old_code_commit.txt)
    tx_strand  -- current working tree: transcript span, same strand, >=30% of shorter
    cds_strand -- CDS span, same strand, >=30% of shorter
    cds_blind  -- CDS span, either strand, >=30% of shorter
  training set (selectTrainingModels) from each, with the NEW selection
  (R3 complete-ORF filter + transitive overlap removal), plus the fully old
  pipeline (old getBestModel -> old selectTrainingModels) as the baseline.
selectTrainingModels gets an empty filterGeneMark keeper GTF (no predict hints
exist at this stage), so keeperCheck is False for every variant.

Each model set is scored against RefSeq protein-coding mRNA CDS chains.
"""
import collections, csv, gzip, logging, os, shutil, sys, tempfile

B = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, B)
logging.basicConfig(level=logging.WARNING)
from funannotate import library as lib
from funannotate import train
import library_old
import train_old
import diversity
from funannotate.interlap import InterLap

lib.log = logging.getLogger("new")
library_old.log = logging.getLogger("old")

R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
# BENCH_RES overrides the PASA results root (e.g. r1_fix/ for the selection x PASA-fix cross)
RES = os.environ.get("BENCH_RES", f"{R}/do_pasa_rust_vs_perl/results")
CASES = {
    "Neurospora_crassa_OR74A": ("GCF_000182925.2_NC12", "GCF_000182925.2"),
    "Aspergillus_nidulans_FGSC_A4": ("GCF_000011425.1_ASM1142v1", "GCF_000011425.1"),
    "Botrytis_cinerea_B05.10": ("GCF_000143535.2_ASM14353v4", "GCF_000143535.2"),
    # low-keeper production genomes (<200 filterGeneMark keepers, reads map >=98%), DECISIONS D29
    "Cryptococcus_neoformans_H99": ("GCF_000149245.1_CNA3", "GCF_000149245.1"),
    "Schizophyllum_commune_H4-8": ("GCF_000143185.2_Schco3", "GCF_000143185.2"),
}
# BENCH_CASES (comma-separated genome names) restricts the run to those genomes
if os.environ.get("BENCH_CASES"):
    CASES = {k: v for k, v in CASES.items() if k in os.environ["BENCH_CASES"].split(",")}
OVERLAP = 30  # funannotate train default --pasa_alignment_overlap


def load_refseq(gff):
    mrna_gene, cds = {}, collections.defaultdict(list)
    coding_genes = set()
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
            if c[2] == "gene" and a.get("gene_biotype") == "protein_coding":
                coding_genes.add(a["ID"])
            elif c[2] == "mRNA":
                mrna_gene[a["ID"]] = a.get("Parent")
            elif c[2] == "CDS" and a.get("Parent"):
                cds[a["Parent"]].append((c[0], c[6], int(c[3]), int(c[4])))
    chains, gene_of_chain = {}, {}
    inter = collections.defaultdict(InterLap)
    for mid, parts in cds.items():
        gid = mrna_gene.get(mid)
        if gid not in coding_genes:
            continue
        parts.sort(key=lambda x: x[2])
        key = (parts[0][0], parts[0][1], tuple((s, e) for _, _, s, e in parts))
        chains[key] = gid
        for ctg, strand, s, e in parts:
            inter[(ctg, strand)].add((s, e))
    introns = {}
    for (ctg, strand, ex), gid in chains.items():
        if len(ex) > 1:
            introns[(ctg, strand, tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))] = gid
    return chains, introns, inter, coding_genes


def model_chain(v):
    ex = sorted(v["CDS"][0])
    return (v["contig"], v["strand"], tuple(ex))


def score(genes, ref):
    chains, introns, inter, coding = ref
    n = exact = ichain = ovl = none = anti = complete = single = 0
    hit_genes, hit_genes_ichain = set(), set()
    for v in genes.values():
        n += 1
        complete += lib.is_complete_model(v)
        ctg, strand, ex = model_chain(v)
        single += len(ex) == 1
        key = (ctg, strand, ex)
        if key in chains:
            exact += 1
            hit_genes.add(chains[key])
        if len(ex) > 1:
            ik = (ctg, strand, tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))
            if ik in introns:
                ichain += 1
                hit_genes_ichain.add(introns[ik])
        same = any(list(inter[(ctg, strand)].find((s, e))) for s, e in ex)
        other = any(list(inter[(ctg, "-" if strand == "+" else "+")].find((s, e))) for s, e in ex)
        if same:
            ovl += 1
        elif other:
            anti += 1
        else:
            none += 1
    pct = lambda a, b: round(100.0 * a / b, 1) if b else 0.0
    return {
        "models": n, "complete_pct": pct(complete, n), "single_exon_pct": pct(single, n),
        "exact_cds": exact, "exact_pct": pct(exact, n),
        "intron_chain_exact": ichain,
        "same_strand_overlap_pct": pct(ovl, n), "antisense_only": anti, "no_refseq_overlap": none,
        "no_overlap_pct": pct(none, n),
        "refseq_genes_exact": len(hit_genes), "recall_exact_pct": pct(len(hit_genes), len(coding)),
        "refseq_genes_ichain": len(hit_genes_ichain),
    }


def best_by_span(step1, fasta, tpm, out, use_cds, strand_aware):
    """Current getBestModel logic with a selectable locus span."""
    Expression = {}
    with open(tpm) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t, g, loc, v = line.rstrip("\n").split("\t")
            Expression[g] = max(float(v), Expression.get(g, 0.0))
    _, Genes = lib.gff2interlap(step1, fasta)

    def rank(gid):
        g = Genes[gid]
        best = (False, 0, 0)
        for i, c in enumerate(g.get("CDS") or []):
            best = max(best, (lib.is_complete_model(g, i), len(c), sum(e[1] - e[0] + 1 for e in c)))
        return best + (Expression.get(gid, 0.0), gid)

    def span(v):
        if not use_cds:
            return v["location"][0], v["location"][1]
        xs = [x for e in v["CDS"][0] for x in e]
        return min(xs), max(xs)

    models = [(k, v["contig"], v["strand"]) + span(v) for k, v in Genes.items()]
    keep = {max(c, key=rank) for c in lib.cluster_overlapping(models, OVERLAP / 100.0, strand_aware)}
    lib.dict2gff3({k: Genes[k] for k in Genes if k in keep}, out)


def main(outdir):
    rows = []
    for name, (asm, acc) in CASES.items():
        d = f"{RES}/{name}.rust"
        if not os.path.exists(f"{d}/pasa.step1.gff3"):
            print(f"skip {name}: no PASA run yet", flush=True)
            continue
        ref = load_refseq(f"{R}/do_pasa_rust_vs_perl/refseq/{acc}.gff")
        tmp = tempfile.mkdtemp(dir=os.environ.get("SCRATCH", "/tmp"))
        g = f"{tmp}/g.fa"
        with gzip.open(f"{R}/input_clean_genomes/{asm}.masked.fasta.gz", "rb") as i, open(g, "wb") as o:
            shutil.copyfileobj(i, o)
        gtf = f"{tmp}/empty.gtf"
        open(gtf, "w").close()
        step1, tpm = f"{d}/pasa.step1.gff3", f"{d}/kallisto.tsv"
        best = {
            "old": f"{outdir}/{name}.best.old.gff3",
            "tx_strand": f"{outdir}/{name}.best.tx_strand.gff3",
            "cds_strand": f"{outdir}/{name}.best.cds_strand.gff3",
            "cds_blind": f"{outdir}/{name}.best.cds_blind.gff3",
        }
        train_old.getBestModel(step1, g, tpm, best["old"], pasa_alignment_overlap=OVERLAP)
        train.getBestModel(step1, g, tpm, best["tx_strand"], pasa_alignment_overlap=OVERLAP)
        best_by_span(step1, g, tpm, best["cds_strand"], use_cds=True, strand_aware=True)
        best_by_span(step1, g, tpm, best["cds_blind"], use_cds=True, strand_aware=False)
        sets = {("pasa_step1_all", "-"): step1}
        for k, p in best.items():
            sets[("getBestModel", k)] = p
            t = f"{outdir}/{name}.train.{k}.newselect.gff3"
            lib.selectTrainingModels(p, g, gtf, t, tmp, min_models=200)
            sets[("training_set", k + "+new_select")] = t
        t = f"{outdir}/{name}.train.old+old_select.gff3"
        library_old.selectTrainingModels(best["old"], g, gtf, t, tmp, min_models=200)
        sets[("training_set", "old+old_select")] = t
        for (stage, variant), p in sets.items():
            G = {}
            G = lib.gff2dict(p, g, G)
            r = {"genome": name, "stage": stage, "variant": variant}
            r.update(score(G, ref))
            r.update(diversity.diversity(G, ref))
            rows.append(r)
            print("\t".join(str(x) for x in r.values()), flush=True)
        shutil.rmtree(tmp, ignore_errors=True)
    with open(f"{outdir}/benchmark.tsv", "w") as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1])
