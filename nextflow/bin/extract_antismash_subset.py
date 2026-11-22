#!/usr/bin/env python3
"""Reduce one antiSMASH (v7/v8) JSON result to the small tables an MS-feature analysis needs.

Keeps: regions (class, coordinates, contig edge, protocluster core), CDS inside regions
(locus tag, protein id, coordinates, gene kind, function text, protein length), polyketide and
nonribosomal domains and modules, biosynthetic-gene protein sequences, and the top
ClusterBlast hits per region.
Drops: genome sequence, all CDS outside regions, all PFAM_domain/gene/mRNA features, the
ClusterBlast protein alignments, HTML, per-region GenBank copies, the zip.

Usage: extract_antismash_subset.py INPUT.json[.gz] OUTDIR [GENOME_ID] [--minimal]
--minimal writes only regions, cds (genes in regions) and known (MIBiG top hits): the BFD representation.
Writes OUTDIR/<GENOME_ID>.{regions,cds,domains,modules,blast_top}.tsv.gz and
OUTDIR/<GENOME_ID>.biosyn_proteins.faa.gz. Reads plain or gzip JSON.
"""
import gzip, json, os, re, sys, pathlib

TOPN = 3


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def span(loc):
    nums = [int(x) for x in re.findall(r"\d+", loc)]
    return (min(nums), max(nums)) if nums else (None, None)


def strand(loc):
    return "-" if loc.rstrip().endswith("(-)") or "(-)" in loc.split(",")[0] else "+"


def q1(f, key, default=""):
    v = f["qualifiers"].get(key)
    if v is None and key == "locus_tag":  # GFF3 input: CDS carry "gene" instead of locus_tag
        v = f["qualifiers"].get("gene")
    if isinstance(v, list):
        return v[0] if v else default
    return default if v is None else v


def main(path, outdir, genome=None, minimal=False):
    path = pathlib.Path(path); outdir = pathlib.Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    genome = genome or path.name.split(".")[0]
    with opener(path) as fh:
        data = json.load(fh)
    keys = ("regions", "cds", "known") if minimal else ("regions", "cds", "domains", "modules", "blast_top", "adomain", "known")
    out = {k: gzip.open(outdir / f"{genome}.{k}.tsv.gz", "wt") for k in keys}
    faa = gzip.open(outdir / f"{genome}.biosyn_proteins.faa.gz", "wt") if not minimal else open(os.devnull, "w")
    NULL = open(os.devnull, "w")
    Z = lambda k: out.get(k, NULL)
    out["regions"].write("genome\trecord\tregion\tstart\tend\tlength\tproducts\tcategories\tcontig_edge\tcore_start\tcore_end\tcandidate_kinds\tn_cds\tn_biosynthetic\tn_additional\n")
    out["cds"].write("genome\trecord\tregion\tlocus_tag\tprotein_id\tstart\tend\tstrand\tgene_kind\tprotein_len\tfunctions\tproduct\n")
    Z("domains").write("genome\trecord\tlocus_tag\tdomain\tsubtypes\tprotein_start\tprotein_end\tevalue\tscore\n")
    Z("modules").write("genome\trecord\tmodule_type\tlocus_tags\tdomains\tcomplete\titerative\n")
    Z("blast_top").write("genome\trecord\tregion\trank\thit_accession\thit_label\thit_cluster_type\thit_description\tsimilarity\tcore_gene_hits\tgene_hits\tblast_score\ttotal_hits\n")
    Z("adomain").write("genome\trecord\tdomain_id\ttool\taa10\tstachelhaus_top_substrates\taa10_score\taa34_score\tphysicochemical_class\tlarge_cluster\tsmall_cluster\n")
    out["known"].write("genome\trecord\tregion\trank\tmibig_accession\thit_description\thit_cluster_type\tsimilarity\tcore_gene_hits\tgene_hits\tblast_score\ttotal_hits\n")
    n = {"known": 0, "adomain": 0, "regions": 0, "cds": 0, "domains": 0, "modules": 0, "blast_top": 0, "faa": 0}
    for rec in data["records"]:
        rid = rec["id"]
        regions = []
        for f in rec["features"]:
            if f["type"] == "region":
                s, e = span(f["location"])
                regions.append((int(q1(f, "region_number", 0)), s, e, f))
        protos = [f for f in rec["features"] if f["type"] == "protocluster"]
        cands = [f for f in rec["features"] if f["type"] == "cand_cluster"]
        cds_in = {}  # region number -> list
        for f in rec["features"]:
            if f["type"] != "CDS":
                continue
            s, e = span(f["location"])
            for num, rs, re_, _ in regions:
                if s is not None and s >= rs and e <= re_:
                    cds_in.setdefault(num, []).append((f, s, e))
                    break
        in_region_tags = set()
        for num, lst in cds_in.items():
            for f, s, e in lst:
                in_region_tags.add(q1(f, "locus_tag"))
        for num, rs, re_, rf in sorted(regions, key=lambda x: x[0]):
            ps = [p for p in protos if rs <= span(p["location"])[0] and span(p["location"])[1] <= re_]
            cats = sorted({q1(p, "category") for p in ps})
            cores = [span(q1(p, "core_location")) for p in ps if q1(p, "core_location")]
            cs = min((c[0] for c in cores), default="")
            ce = max((c[1] for c in cores), default="")
            kinds = sorted({q1(c, "kind") for c in cands if rs <= span(c["location"])[0] and span(c["location"])[1] <= re_})
            lst = cds_in.get(num, [])
            nb = sum(1 for f, _, _ in lst if q1(f, "gene_kind") == "biosynthetic")
            na = sum(1 for f, _, _ in lst if q1(f, "gene_kind") == "biosynthetic-additional")
            out["regions"].write("\t".join(map(str, [genome, rid, num, rs, re_, re_ - rs, ";".join(rf["qualifiers"].get("product", [])),
                                                       ";".join(cats), q1(rf, "contig_edge"), cs, ce, ";".join(kinds), len(lst), nb, na])) + "\n"); n["regions"] += 1
            for f, s, e in lst:
                tr = q1(f, "translation")
                out["cds"].write("\t".join(map(str, [genome, rid, num, q1(f, "locus_tag"), q1(f, "protein_id"), s, e, strand(f["location"]),
                                                      q1(f, "gene_kind"), len(tr), "|".join(f["qualifiers"].get("gene_functions", [])).replace("\t", " "),
                                                      q1(f, "product").replace("\t", " ")])) + "\n"); n["cds"] += 1
                if q1(f, "gene_kind") == "biosynthetic" and tr:
                    faa.write(f">{genome}|{q1(f, 'locus_tag')}|{rid}|region{num}\n{tr}\n"); n["faa"] += 1
        for f in rec["features"]:
            if f["type"] == "aSDomain" and q1(f, "locus_tag") in in_region_tags:
                Z("domains").write("\t".join(map(str, [genome, rid, q1(f, "locus_tag"), q1(f, "aSDomain"), ";".join(f["qualifiers"].get("domain_subtypes", [])),
                                                          q1(f, "protein_start"), q1(f, "protein_end"), q1(f, "evalue"), q1(f, "score")])) + "\n"); n["domains"] += 1
            elif f["type"] == "aSModule":
                q = f["qualifiers"]
                Z("modules").write("\t".join(map(str, [genome, rid, ";".join(q.get("type", [])), ";".join(q.get("locus_tags", [])),
                                                          ";".join(x.split("_")[-1] if False else x for x in q.get("domains", [])), q.get("complete"), q.get("iterative")])) + "\n"); n["modules"] += 1
        cb = rec.get("modules", {}).get("antismash.modules.clusterblast", {}).get("general", {}).get("results", [])
        for res in cb:
            for rank, item in enumerate(res.get("ranking", [])[:TOPN], 1):
                cl, sc = item[0], item[1]
                Z("blast_top").write("\t".join(map(str, [genome, rid, res.get("region_number"), rank, cl.get("accession"), cl.get("cluster_label"),
                                                            cl.get("cluster_type"), str(cl.get("description", "")).replace("\t", " "), sc.get("similarity"),
                                                            sc.get("core_gene_hits"), sc.get("hits"), sc.get("blast_score"), res.get("total_hits")])) + "\n"); n["blast_top"] += 1
        kcl = rec.get("modules", {}).get("antismash.modules.clusterblast", {}).get("knowncluster", {}).get("results", [])
        for res in kcl:
            for rank, item in enumerate(res.get("ranking", [])[:TOPN], 1):
                cl, sc = item[0], item[1]
                out["known"].write("\t".join(map(str, [genome, rid, res.get("region_number"), rank, cl.get("accession"), str(cl.get("description", "")).replace("\t", " "),
                                                         cl.get("cluster_type"), sc.get("similarity"), sc.get("core_gene_hits"), sc.get("hits"), sc.get("blast_score"), res.get("total_hits")])) + "\n"); n["known"] += 1
        preds = rec.get("modules", {}).get("antismash.modules.nrps_pks", {}).get("domain_predictions", {})
        for did, tools in preds.items():
            for tool, v in tools.items():
                if tool != "nrpys" or not isinstance(v, dict):
                    continue
                st = (v.get("stachelhaus_matches") or [{}])[0]
                subs = ",".join(x.get("short", "") for x in st.get("substrates", []))
                Z("adomain").write("\t".join(map(str, [genome, rid, did, tool, v.get("aa10", ""), subs, st.get("aa10_score", ""), st.get("aa34_score", ""),
                                                          (v.get("physiochemical_class") or {}).get("name", ""), (v.get("large_cluster") or {}).get("name", ""),
                                                          (v.get("small_cluster") or {}).get("name", "")])) + "\n"); n["adomain"] += 1
    for h in out.values():
        h.close()
    faa.close()
    print(genome, n)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--minimal"]
    if len(args) < 2:
        sys.exit(__doc__)
    main(*args[:3], minimal="--minimal" in sys.argv)
