#!/usr/bin/env python3
"""Compare funannotate predictions in genome_annotation.pre_rc6 (old) and genome_annotation (rc.6, new).

Three steps:
  pairs      Find genomes that have a GFF3 in both trees. Optionally pick a pilot sample
             stratified by phylum.
  compare    Per genome: gene, mRNA, exon, CDS and intron counts for old and new, the
             same counts restricted to scaffolds that are identical in both assemblies,
             and gffcompare agreement (old = reference, new = query).
  summarize  Per-phylum summary and the distribution of the new/old gene-count ratio.
  reference  Score the old and new predictions against a reference annotation (for example
             NCBI RefSeq) with gffcompare, in two modes: "exon" (the exons as annotated,
             UTRs included) and "cds" (CDS segments as exons, which is the fair mode when the
             prediction has no UTRs). Reference = gffcompare -r.

Why restrict to shared scaffolds: the old and new predictions can come from different
versions of the assembly (the cleaning step drops scaffolds). Genes on a scaffold that is
absent from one assembly cannot agree, so they are reported but not scored.
A scaffold is "shared" when its name and length are the same in both FASTA files.
If a FASTA is missing for either side, shared scaffolds are the seqids that have features
in both GFFs. That is a weaker rule (a scaffold with no genes in one set is dropped); the
column `scaffold_source` says which rule was used.

Only protein-coding genes (genes with an mRNA child) go to gffcompare. tRNA genes are
counted separately.

Reads .gff3, .gff3.gz and .gff3.zst, and the same for FASTA. Needs `gffcompare` on PATH
(module load gffcompare/0.12.10). Run python3.12.

make_sample_tag() copies the function of the same name in
nextflow/bin/pick_representative_strain.py (the Python twin of makeSampleTag in
nextflow/modules/common/utils.nf). Keep them in sync.
"""
import argparse
import csv
import gzip
import io
import math
import multiprocessing as mp
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

OLD_DIR = "genome_annotation.pre_rc6"
NEW_DIR = "genome_annotation"


# ── helpers ───────────────────────────────────────────────────────────────────
def clean_strain(raw_strain):
    s = (raw_strain or "").strip()
    s = re.sub(r"['\"]", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species, raw_strain):
    sp = re.sub(r"['\"]", "", (raw_species or "").strip())
    st = clean_strain(raw_strain)
    tag = "_".join(p for p in (sp, st) if p)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def open_text(path):
    """Open a plain, .gz or .zst text file for reading."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    if p.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", p], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout)
    return open(p, "rt")


def first_existing(*paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return ""


def gff_path(root, tree, out):
    base = os.path.join(root, tree, out, "predict_results", out + ".gff3")
    return first_existing(base, base + ".gz", base + ".zst")


def fasta_path(gff):
    """Scaffold FASTA that sits next to a GFF in predict_results."""
    base = re.sub(r"\.gff3(\.gz|\.zst)?$", "", gff) + ".scaffolds.fa"
    return first_existing(base, base + ".gz", base + ".zst")


def scaffold_lengths(fa):
    """{name: length} from a .fai when one exists, otherwise by streaming the FASTA."""
    if not fa:
        return None
    for fai in (fa + ".fai", re.sub(r"\.(gz|zst)$", "", fa) + ".fai"):
        if os.path.exists(fai):
            lens = {}
            with open(fai) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    lens[f[0]] = int(f[1])
            return lens
    lens, name, n = {}, None, 0
    with open_text(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def attrs(col9):
    d = {}
    for kv in col9.rstrip("\n").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


# ── GFF parsing ───────────────────────────────────────────────────────────────
class Gff:
    """Gene models of one GFF3. mRNA and exon lines are kept to write a restricted GFF."""

    def __init__(self, path):
        self.genes = {}                       # id -> (seqid, start, end)
        self.mrna = {}                        # id -> (gene_id, seqid, start, end, line)
        self.exons = defaultdict(list)        # mrna id -> [(start, end, line)]
        self.cds = defaultdict(int)           # mrna id -> CDS segment count
        self.cds_segs = defaultdict(list)     # mrna id -> [(start, end)]
        self.trna = []                        # seqid of each tRNA feature
        self.seqids = set()
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9:
                    continue
                seqid, typ, s, e = c[0], c[2], int(c[3]), int(c[4])
                self.seqids.add(seqid)
                a = attrs(c[8])
                if typ == "gene":
                    self.genes[a.get("ID", "")] = (seqid, s, e)
                elif typ == "mRNA":
                    self.mrna[a.get("ID", "")] = (a.get("Parent", "").split(",")[0], seqid, s, e, line)
                elif typ == "exon":
                    self.exons[a.get("Parent", "").split(",")[0]].append((s, e, line))
                elif typ == "CDS":
                    pid = a.get("Parent", "").split(",")[0]
                    self.cds[pid] += 1
                    self.cds_segs[pid].append((s, e))
                elif typ == "tRNA":
                    self.trna.append(seqid)

    def metrics(self, keep=None):
        """Counts over protein-coding genes on the seqids in `keep` (None = all)."""
        genes, n_mrna, n_exon, n_cds, n_intron, single = set(), 0, 0, 0, 0, 0
        for mid, (gid, seqid, s, e, _) in self.mrna.items():
            if keep is not None and seqid not in keep:
                continue
            genes.add(gid)
            n_mrna += 1
            ne = len(self.exons.get(mid, ()))
            n_exon += ne
            n_cds += self.cds.get(mid, 0)
            n_intron += max(ne - 1, 0)
            single += 1 if ne == 1 else 0
        lens = []
        for gid in genes:
            if gid in self.genes:
                lens.append(self.genes[gid][2] - self.genes[gid][1] + 1)
        return {
            "genes": len(genes), "mrna": n_mrna, "exons": n_exon, "cds": n_cds,
            "introns": n_intron,
            "mean_gene_len": round(sum(lens) / len(lens), 1) if lens else "",
            "single_exon_frac": round(single / n_mrna, 4) if n_mrna else "",
            "trna": sum(1 for s in self.trna if keep is None or s in keep),
        }

    def introns(self, keep):
        """One record per intron of every mRNA with >= 2 CDS segments on the seqids in `keep`.
        Introns lie between CDS segments. Numbering and phase follow the transcript (5' to 3').
        phase = coding bases 5' of the intron mod 3 (0 = between codons). The CDS phase column of
        the first segment is not used, so models with a partial 5' end can be off."""
        out = []
        for mid, (gid, seqid, s, e, line) in self.mrna.items():
            if seqid not in keep:
                continue
            segs = sorted(self.cds_segs.get(mid, ()))
            if len(segs) < 2:
                continue
            strand = line.split("\t")[6]
            total = sum(b - a + 1 for a, b in segs)
            n = len(segs) - 1
            order = list(range(len(segs))) if strand != "-" else list(range(len(segs) - 1, -1, -1))
            upstream = 0
            for rank, i in enumerate(order[:-1]):
                a, b = segs[i]
                upstream += b - a + 1
                nxt = segs[order[rank + 1]]
                lo, hi = (b + 1, nxt[0] - 1) if strand != "-" else (nxt[1] + 1, a - 1)
                out.append({"transcript": mid, "gene": gid, "seqid": seqid, "start": lo, "end": hi,
                            "strand": strand, "intron_number": rank + 1, "n_introns": n,
                            "codon_position": upstream // 3, "phase": upstream % 3,
                            "cds_length": total})
        return out

    def write_restricted(self, path, keep, mode="exon"):
        """mRNA and exon lines of protein-coding genes on the seqids in `keep`.
        mode "cds": the CDS segments are written as the exons (UTRs dropped)."""
        n = 0
        with open(path, "w") as out:
            out.write("##gff-version 3\n")
            for mid, (gid, seqid, s, e, line) in self.mrna.items():
                if seqid not in keep:
                    continue
                if mode == "cds":
                    segs = sorted(self.cds_segs.get(mid, ()))
                    if not segs:
                        continue
                    c = line.rstrip("\n").split("\t")
                    strand = c[6]
                    c[3], c[4], c[8] = str(segs[0][0]), str(segs[-1][1]), f"ID={mid};Parent={gid}"
                    out.write("\t".join(c) + "\n")
                    for i, (cs, ce) in enumerate(segs):
                        out.write("\t".join([seqid, "ref", "exon", str(cs), str(ce), ".", strand, ".",
                                             f"ID={mid}.e{i};Parent={mid}"]) + "\n")
                else:
                    out.write(line)
                    for _, _, el in self.exons.get(mid, ()):
                        out.write(el)
                n += 1
        return n


# ── gffcompare ────────────────────────────────────────────────────────────────
LEVELS = [("Base", "base"), ("Exon", "exon"), ("Intron", "intron"),
          ("Intron chain", "ichain"), ("Transcript", "transcript"), ("Locus", "locus")]


def run_gffcompare(ref, qry, workdir):
    """gffcompare -r ref qry. Returns a dict of Sn/Pr and match counts."""
    r = subprocess.run(["gffcompare", "-r", ref, "-o", "cmp", qry], cwd=workdir,
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError("gffcompare failed: " + r.stderr.strip()[-200:])
    res = {}
    with open(os.path.join(workdir, "cmp.stats")) as fh:
        txt = fh.read()
    for label, key in LEVELS:
        m = re.search(r"^\s*" + re.escape(label) + r" level:\s+([\d.]+)\s*\|\s*([\d.]+)", txt, re.M)
        res[f"gc_{key}_sn"] = float(m.group(1)) if m else ""
        res[f"gc_{key}_pr"] = float(m.group(2)) if m else ""
    for label, key in (("Matching intron chains", "match_ichains"),
                       ("Matching transcripts", "match_transcripts")):
        m = re.search(label + r":\s+(\d+)", txt)
        res[key] = int(m.group(1)) if m else ""
    for label, key in (("Missed exons", "missed_exons"), ("Novel exons", "novel_exons"),
                       ("Missed introns", "missed_introns"), ("Novel introns", "novel_introns")):
        m = re.search(label + r":\s+(\d+)/(\d+)", txt)
        res[key] = int(m.group(1)) if m else ""
        res[key + "_of"] = int(m.group(2)) if m else ""
    # Transcript classes. tmap: one row per query transcript. refmap: one row per
    # reference transcript that has a match (class =, c, k, m, n, j, e, o).
    q = defaultdict(int)
    for fn in os.listdir(workdir):
        if fn.endswith(".tmap"):
            with open(os.path.join(workdir, fn)) as fh:
                next(fh, None)
                for line in fh:
                    f = line.split("\t")
                    if len(f) > 2:
                        q[f[2]] += 1
        elif fn.endswith(".refmap"):
            r_ref = defaultdict(int)
            with open(os.path.join(workdir, fn)) as fh:
                next(fh, None)
                for line in fh:
                    f = line.split("\t")
                    if len(f) > 2:
                        r_ref[f[2]] += 1
            res["ref_mrna_matched"] = sum(r_ref.values())
            res["ref_mrna_exact"] = r_ref.get("=", 0)
    res["pred_class_exact"] = q.get("=", 0)
    res["pred_class_no_overlap"] = q.get("u", 0)
    res["pred_class_total"] = sum(q.values())
    return res


# ── compare one genome ────────────────────────────────────────────────────────
def compare_one(job):
    row, tmproot, keep_tmp = job
    out = {"out": row["out"], "locustag": row["locustag"], "asmid": row["asmid"],
           "phylum": row["phylum"], "class": row["class"], "error": ""}
    tmp = tempfile.mkdtemp(prefix="cmp_", dir=tmproot)
    try:
        old, new = Gff(row["old_gff"]), Gff(row["new_gff"])
        fa_old, fa_new = fasta_path(row["old_gff"]), fasta_path(row["new_gff"])
        len_old, len_new = scaffold_lengths(fa_old), scaffold_lengths(fa_new)
        if len_old is not None and len_new is not None:
            shared = {k for k, v in len_old.items() if len_new.get(k) == v}
            out["scaffold_source"] = "fasta"
            out["scaffolds_old"], out["scaffolds_new"] = len(len_old), len(len_new)
            out["scaffolds_shared"] = len(shared)
            out["scaffolds_len_mismatch"] = sum(1 for k, v in len_old.items()
                                                if k in len_new and len_new[k] != v)
        else:
            shared = old.seqids & new.seqids
            out["scaffold_source"] = "gff_seqids"
            out["scaffolds_old"], out["scaffolds_new"] = len(old.seqids), len(new.seqids)
            out["scaffolds_shared"], out["scaffolds_len_mismatch"] = len(shared), ""
        for tag, g in (("old", old), ("new", new)):
            for k, v in g.metrics().items():
                out[f"{tag}_{k}"] = v
            for k, v in g.metrics(shared).items():
                out[f"{tag}_shared_{k}"] = v
        for tag in ("old", "new"):
            tot, sh = out[f"{tag}_genes"], out[f"{tag}_shared_genes"]
            out[f"{tag}_frac_genes_on_shared"] = round(sh / tot, 4) if tot else ""
        o, n = out["old_genes"], out["new_genes"]
        out["gene_ratio_new_old"] = round(n / o, 4) if o else ""
        o, n = out["old_shared_genes"], out["new_shared_genes"]
        out["gene_ratio_shared"] = round(n / o, 4) if o else ""
        n_old = old.write_restricted(os.path.join(tmp, "old.gff3"), shared)
        n_new = new.write_restricted(os.path.join(tmp, "new.gff3"), shared)
        if n_old and n_new:
            out.update(run_gffcompare(os.path.join(tmp, "old.gff3"),
                                      os.path.join(tmp, "new.gff3"), tmp))
    except Exception as e:  # keep going: one bad genome must not stop the run
        out["error"] = f"{type(e).__name__}: {e}"[:300]
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return out


# ── pairs ─────────────────────────────────────────────────────────────────────
PAIR_COLS = ["out", "locustag", "asmid", "phylum", "class", "old_gff", "new_gff",
             "has_old", "has_new", "duplicate_out"]


def cmd_pairs(a):
    rows, seen = [], defaultdict(int)
    with open(a.samples, newline="") as fh:
        for r in csv.DictReader(fh):
            species = (r.get("SPECIES") or "").strip().strip("'\"")
            if not species:
                continue
            out = make_sample_tag(species, r.get("STRAIN") or "")
            seen[out] += 1
            rows.append({"out": out, "locustag": (r.get("LOCUSTAG") or "").strip(),
                         "asmid": (r.get("ASMID") or "").strip(),
                         "phylum": (r.get("PHYLUM") or "").strip() or "unknown",
                         "class": (r.get("CLASS") or "").strip()})

    def look(r):
        r["old_gff"] = gff_path(a.runs_root, OLD_DIR, r["out"])
        r["new_gff"] = gff_path(a.runs_root, NEW_DIR, r["out"])
        r["has_old"], r["has_new"] = int(bool(r["old_gff"])), int(bool(r["new_gff"]))
        r["duplicate_out"] = int(seen[r["out"]] > 1)
        return r

    with ThreadPoolExecutor(a.threads) as ex:
        rows = list(ex.map(look, rows))
    write_tsv(a.out, rows, PAIR_COLS)
    both = [r for r in rows if r["has_old"] and r["has_new"] and not r["duplicate_out"]]
    print(f"samples {len(rows)}; old gff {sum(r['has_old'] for r in rows)}; "
          f"new gff {sum(r['has_new'] for r in rows)}; both (no duplicate out) {len(both)}")
    if a.pilot:
        pilot = stratified(both, a.pilot, a.min_per_phylum, a.seed)
        write_tsv(a.pilot_out, pilot, PAIR_COLS)
        by = defaultdict(int)
        for r in pilot:
            by[r["phylum"]] += 1
        print(f"pilot {len(pilot)} genomes -> {a.pilot_out}")
        for k, v in sorted(by.items(), key=lambda kv: -kv[1]):
            print(f"  {k}: {v}")


def stratified(rows, n, min_per, seed):
    """Sample n rows. Each phylum gets a share proportional to its size, at least
    min_per (or all of it if smaller), by largest remainder."""
    groups = defaultdict(list)
    for r in rows:
        groups[r["phylum"]].append(r)
    total = len(rows)
    alloc = {k: min(len(v), max(min_per, int(n * len(v) / total))) for k, v in groups.items()}
    rem = n - sum(alloc.values())
    order = sorted(groups, key=lambda k: -(n * len(groups[k]) / total - alloc[k]))
    i = 0
    while rem > 0 and any(alloc[k] < len(groups[k]) for k in groups):
        k = order[i % len(order)]
        if alloc[k] < len(groups[k]):
            alloc[k] += 1
            rem -= 1
        i += 1
    while rem < 0:  # min_per pushed the total over n: trim the largest groups
        k = max(alloc, key=lambda k: alloc[k])
        alloc[k] -= 1
        rem += 1
    rng = random.Random(seed)
    picked = []
    for k, v in groups.items():
        picked += rng.sample(v, alloc[k])
    return sorted(picked, key=lambda r: r["out"])


def write_tsv(path, rows, cols=None):
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ── compare ───────────────────────────────────────────────────────────────────
def cmd_compare(a):
    if not shutil.which("gffcompare"):
        sys.exit("ERROR: gffcompare not on PATH. Run: module load gffcompare/0.12.10")
    with open(a.pairs, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["old_gff"] and r["new_gff"] and r.get("duplicate_out") != "1"]
    rows.sort(key=lambda r: r["out"])
    if a.shard:
        i, n = (int(x) for x in a.shard.split("/"))
        if not 0 <= i < n:
            sys.exit("ERROR: --shard I/N needs 0 <= I < N")
        rows = rows[i::n]  # interleaved, so large and small genomes spread over shards
    tmproot = a.tmpdir or os.environ.get("SCRATCH") or tempfile.gettempdir()
    os.makedirs(tmproot, exist_ok=True)
    print(f"{len(rows)} genomes, {a.procs} processes, temp dir {tmproot}", flush=True)
    results = []
    with mp.Pool(a.procs) as pool:
        for i, res in enumerate(pool.imap_unordered(
                compare_one, [(r, tmproot, a.keep_tmp) for r in rows], chunksize=1), 1):
            results.append(res)
            if i % 25 == 0 or i == len(rows):
                print(f"  {i}/{len(rows)}", flush=True)
    results.sort(key=lambda r: r["out"])
    write_tsv(a.out, results)
    bad = [r for r in results if r["error"]]
    print(f"wrote {a.out}; {len(bad)} genomes with errors")
    for r in bad[:10]:
        print(f"  {r['out']}: {r['error']}")


# ── reference ─────────────────────────────────────────────────────────────────
def reference_one(job):
    row, which, tmproot, keep_tmp, intron_dir, save_dir = job
    pred_path = row[f"{which}_gff"]
    base = {"out": row["out"], "pred": which, "ref_gff": os.path.basename(row["ref_gff"]), "error": ""}
    tmp = tempfile.mkdtemp(prefix="ref_", dir=tmproot)
    rows = []
    try:
        ref, pred = Gff(row["ref_gff"]), Gff(pred_path)
        lens = scaffold_lengths(fasta_path(pred_path))
        if lens is not None:
            keep = set(lens) & ref.seqids
            src = "pred_fasta_names"
        else:
            keep = pred.seqids & ref.seqids
            src = "gff_seqids"
        info = dict(base, scaffold_source=src, ref_seqids=len(ref.seqids),
                    pred_scaffolds=len(lens) if lens is not None else len(pred.seqids),
                    seqids_used=len(keep))
        for tag, g in (("ref", ref), ("pred", pred)):
            for k, v in g.metrics(keep).items():
                info[f"{tag}_{k}"] = v
        info["gene_ratio_pred_ref"] = (round(info["pred_genes"] / info["ref_genes"], 4)
                                       if info["ref_genes"] else "")
        if intron_dir:
            write_intron_table(os.path.join(intron_dir, f"{row['out']}.{which}.introns.tsv.gz"),
                               ref.introns(keep), pred.introns(keep))
        for mode in ("exon", "cds"):
            r = dict(info, mode=mode)
            d = os.path.join(tmp, mode)
            os.makedirs(d)
            n_ref = ref.write_restricted(os.path.join(d, "ref.gff3"), keep, mode)
            n_pred = pred.write_restricted(os.path.join(d, "pred.gff3"), keep, mode)
            r["ref_mrna_scored"], r["pred_mrna_scored"] = n_ref, n_pred
            if n_ref and n_pred:
                r.update(run_gffcompare(os.path.join(d, "ref.gff3"), os.path.join(d, "pred.gff3"), d))
                if save_dir:
                    os.makedirs(save_dir, exist_ok=True)
                    for fn in os.listdir(d):
                        ext = fn.rsplit(".", 1)[-1]
                        if ext in ("tmap", "refmap", "stats"):
                            with open(os.path.join(d, fn), "rb") as src, gzip.open(
                                    os.path.join(save_dir, f"{row['out']}.{which}.{mode}.{ext}.gz"), "wb") as dst:
                                shutil.copyfileobj(src, dst)
                if r.get("ref_mrna_matched", "") != "":
                    r["ref_mrna_missed"] = n_ref - r["ref_mrna_matched"]
            rows.append(r)
    except Exception as e:
        rows.append(dict(base, error=f"{type(e).__name__}: {e}"[:300]))
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return rows


INTRON_COLS = ["source", "status", "unique_first", "transcript", "gene", "seqid", "start", "end", "strand",
               "intron_number", "n_introns", "codon_position", "phase", "cds_length"]


def write_intron_table(path, ref_introns, pred_introns):
    """Per-intron outcome table. status: reference introns are matched or missed; predicted introns
    are matched or novel. A match is the same seqid, start, end and strand (gffcompare's rule).
    unique_first = 1 on the first row of each distinct intron within its source, so alternative
    isoforms that share an intron can be counted once."""
    key = lambda r: (r["seqid"], r["start"], r["end"], r["strand"])
    rk, pk = {key(r) for r in ref_introns}, {key(r) for r in pred_introns}
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INTRON_COLS, delimiter="\t")
        w.writeheader()
        for src, rows, other, miss in (("ref", ref_introns, pk, "missed"), ("pred", pred_introns, rk, "novel")):
            seen = set()
            for r in rows:
                k = key(r)
                w.writerow(dict(r, source=src, status="matched" if k in other else miss,
                                unique_first=int(k not in seen)))
                seen.add(k)


def cmd_reference(a):
    if not shutil.which("gffcompare"):
        sys.exit("ERROR: gffcompare not on PATH. Run: module load gffcompare/0.12.10")
    with open(a.table, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k in ("ref_gff", "old_gff", "new_gff"):
            if not r.get(k) or not os.path.exists(r[k]):
                sys.exit(f"ERROR: {r.get('out')}: {k} not found: {r.get(k)}")
    tmproot = a.tmpdir or os.environ.get("SCRATCH") or tempfile.gettempdir()
    os.makedirs(tmproot, exist_ok=True)
    jobs = [(r, w, tmproot, a.keep_tmp, a.introns_dir, a.save_dir) for r in rows
            for w in (a.sets.split(",") if a.sets else ("old", "new"))]
    with mp.Pool(min(a.procs, len(jobs))) as pool:
        results = [x for res in pool.imap(reference_one, jobs) for x in res]
    write_tsv(a.out, results)
    bad = [r for r in results if r.get("error")]
    print(f"wrote {a.out}: {len(results)} rows, {len(bad)} errors")
    for r in bad:
        print(f"  {r['out']} {r.get('pred','')}: {r['error']}")


# ── summarize ─────────────────────────────────────────────────────────────────
def cmd_summarize(a):
    import duckdb
    con = duckdb.connect()
    con.execute(f"CREATE VIEW c AS SELECT * FROM read_csv('{a.compare}', delim='\t', header=true, "
                "sample_size=-1, nullstr='', union_by_name=true) WHERE error IS NULL OR error = ''")
    os.makedirs(a.out_dir, exist_ok=True)
    by_phylum = con.sql("""
        SELECT phylum, count(*) AS n_genomes,
               median(old_genes) AS median_old_genes, median(new_genes) AS median_new_genes,
               round(median(gene_ratio_new_old), 3) AS median_gene_ratio,
               round(median(gene_ratio_shared), 3) AS median_gene_ratio_shared,
               round(median(old_frac_genes_on_shared), 3) AS median_old_frac_on_shared,
               round(median(gc_exon_sn), 1) AS median_exon_sn, round(median(gc_exon_pr), 1) AS median_exon_pr,
               round(median(gc_ichain_sn), 1) AS median_ichain_sn, round(median(gc_ichain_pr), 1) AS median_ichain_pr,
               round(median(gc_locus_sn), 1) AS median_locus_sn, round(median(gc_locus_pr), 1) AS median_locus_pr
        FROM c GROUP BY phylum ORDER BY n_genomes DESC""").df()
    by_phylum.to_csv(os.path.join(a.out_dir, "summary_by_phylum.tsv"), sep="\t", index=False)
    dist = con.sql("""
        SELECT count(*) AS n_genomes,
               round(quantile_cont(gene_ratio_new_old, 0.05), 3) AS p05,
               round(quantile_cont(gene_ratio_new_old, 0.25), 3) AS p25,
               round(quantile_cont(gene_ratio_new_old, 0.50), 3) AS p50,
               round(quantile_cont(gene_ratio_new_old, 0.75), 3) AS p75,
               round(quantile_cont(gene_ratio_new_old, 0.95), 3) AS p95,
               count(*) FILTER (WHERE gene_ratio_new_old < 0.5) AS ratio_below_0_5,
               count(*) FILTER (WHERE gene_ratio_new_old > 2) AS ratio_above_2,
               count(*) FILTER (WHERE scaffold_source = 'gff_seqids') AS shared_by_seqids_only
        FROM c""").df()
    dist.to_csv(os.path.join(a.out_dir, "gene_ratio_distribution.tsv"), sep="\t", index=False)
    print(dist.to_string(index=False))
    print()
    print(by_phylum.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pairs", help="find genomes with a GFF3 in both trees")
    p.add_argument("--samples", required=True)
    p.add_argument("--runs-root", required=True, help="Fungi_BFD_runs directory")
    p.add_argument("--out", required=True, help="all-genomes pair table (TSV)")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--pilot", type=int, default=0, help="also write a pilot sample of N genomes")
    p.add_argument("--pilot-out", default="pilot_pairs.tsv")
    p.add_argument("--min-per-phylum", type=int, default=2)
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(fn=cmd_pairs)
    p = sub.add_parser("compare", help="compare old and new predictions per genome")
    p.add_argument("--pairs", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    p.add_argument("--tmpdir", default=None, help="default: $SCRATCH, then the system temp dir")
    p.add_argument("--keep-tmp", action="store_true")
    p.add_argument("--shard", default="", help="I/N: process every N-th genome starting at I (0-based)")
    p.set_defaults(fn=cmd_compare)
    p = sub.add_parser("reference", help="score old and new predictions against a reference GFF3")
    p.add_argument("--table", required=True, help="TSV with columns out, ref_gff, old_gff, new_gff")
    p.add_argument("--out", required=True)
    p.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    p.add_argument("--tmpdir", default=None)
    p.add_argument("--keep-tmp", action="store_true")
    p.add_argument("--sets", default="", help="old,new (default both) or just one of them")
    p.add_argument("--introns-dir", default="", help="write one per-intron outcome table per genome and set here")
    p.add_argument("--save-dir", default="", help="keep the gffcompare .tmap, .refmap and .stats here (gzip)")
    p.set_defaults(fn=cmd_reference)
    p = sub.add_parser("summarize", help="per-phylum summary")
    p.add_argument("--compare", required=True, help="TSV file or glob (for example shard files)")
    p.add_argument("--out-dir", required=True)
    p.set_defaults(fn=cmd_summarize)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
