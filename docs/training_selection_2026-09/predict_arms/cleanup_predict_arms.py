#!/usr/bin/env python3
"""Reduce predict_arms/ to the records needed to audit or rerun the experiments.

Agreed with REVIEW on 2026-09-26 (DECISIONS D89/D90 and the entry that records this run).

KEEP: run records (arm_info, run_status, score, BUSCO summary/log, capture logs),
funannotate logfiles/, PASA evidence inputs (pasa.gff3, titration pasa.subset.gff3,
pasa_inputs/), trinity BAMs (arm.sh --rna_bam), split lists, GeneMark GTFs,
transcripts GFF3, predict_misc training/evidence records, trained parameters,
final predictions (GFF3, proteins, parameters/stats JSON, validation reports),
scripts and job tables.
REMOVE: stock AUGUSTUS config copied into ab_initio_parameters/augustus
(the trained species dirs, BUSCO_* dirs and snap .hmm stay), genome FASTA copies (prep.sh rebuilds them from the split lists),
per-arm augustus config copies, busco_downloads/, gmes/, EVM partitions,
tbl2asn work dirs, predict_misc intermediates, bulky derivable predict_results
files (.gbk, .scaffolds.fa, *-transcripts.fa, .tbl) and the code_new* snapshots
(diffs against funannotate 41a2fd7, which is on GitHub, are in Fungi_BFD docs).

Kept plain-text files larger than 1 MB are compressed with zstd (original removed
after zstd verifies the frame).

Usage: cleanup_predict_arms.py ROOT OUTDIR [--apply]
Without --apply nothing is changed; OUTDIR gets plan_keep.tsv / plan_remove.tsv.
With --apply the removals and compression run, then manifest_kept.tsv
(path, bytes, sha256) and summary.txt are written.
"""
import hashlib
import os
import shutil
import subprocess
import sys

REMOVE_DIR_NAMES = {"busco_downloads", "gmes", "EVM", "tbl2asn", "__pycache__"}
KEEP_MISC_FILES = {
    "weights.evm.txt", "final_training_models.gff3", "final_training_models.gff3.gz",
    "pasa_predictions.gff3", "pasa.training.tmp.f.good.gtf", "pasa.training.tmp.f.bad.gtf",
    "gene_predictions.gff3", "protein_alignments.gff3",
}
KEEP_MISC_DIRS = {"ab_initio_parameters"}
STOCK_AUGUSTUS_DIRS = {"model", "extrinsic", "profile", "parameters", "config"}
STOCK_SPECIES = {"generic", "anidulans"}
RESULT_REMOVE_SUFFIXES = (".gbk", ".scaffolds.fa", "-transcripts.fa", ".tbl")
GENOME_COPY_NAMES = {"genome.fa", "genome_train.fa", "genome.fa.fai", "genome_train.fa.fai"}
COMPRESSED = (".gz", ".zst", ".bam", ".bai", ".png", ".pdf")
COMPRESS_MIN = 1_000_000


def classify(root):
    """Walk ROOT and return (keep_files, remove_paths) as lists of (path, bytes)."""
    keep, remove = [], []

    def tree_size(p):
        total = 0
        for dp, dn, fn in os.walk(p):
            for f in fn:
                try:
                    total += os.lstat(os.path.join(dp, f)).st_size
                except OSError:
                    pass
        return total

    for dp, dn, fn in os.walk(root):
        base = os.path.basename(dp)
        parts = os.path.relpath(dp, root).split(os.sep)
        if "ab_initio_parameters" in parts:
            # Trained parameters: keep the trained species dirs, BUSCO_* training dirs
            # and snap .hmm. Drop the stock AUGUSTUS config copied into every run
            # (model/extrinsic/profile/parameters and species generic/anidulans are
            # identical to lib/augustus/3.5/config and funannotate_db trained_species).
            sub = parts[parts.index("ab_initio_parameters") + 1:]
            if sub == ["augustus"]:
                drop = [d for d in dn if d in STOCK_AUGUSTUS_DIRS]
            elif sub == ["augustus", "species"]:
                drop = [d for d in dn if d in STOCK_SPECIES]
            else:
                drop = []
            for d in drop:
                q = os.path.join(dp, d)
                remove.append((q, tree_size(q)))
                dn.remove(d)
            keep.extend((os.path.join(dp, f), os.lstat(os.path.join(dp, f)).st_size) for f in fn)
            continue
        parent = os.path.basename(os.path.dirname(dp))
        drop = []
        for d in dn:
            p = os.path.join(dp, d)
            if os.path.islink(p):
                continue
            if d in REMOVE_DIR_NAMES or d.startswith("code_new"):
                drop.append(d)
            elif d == "augustus" and base != "ab_initio_parameters":
                drop.append(d)          # per-arm AUGUSTUS config copy
            elif base == "predict_misc" and d not in KEEP_MISC_DIRS:
                drop.append(d)
        for d in drop:
            p = os.path.join(dp, d)
            remove.append((p, tree_size(p)))
            dn.remove(d)
        for f in fn:
            p = os.path.join(dp, f)
            size = os.lstat(p).st_size
            rm = False
            if base == "predict_misc":
                rm = f not in KEEP_MISC_FILES
            elif base == "predict_results":
                rm = f.endswith(RESULT_REMOVE_SUFFIXES)
            elif f in GENOME_COPY_NAMES:
                rm = True
            (remove if rm else keep).append((p, size))
    return keep, remove


def is_text(p):
    with open(p, "rb") as fh:
        return b"\0" not in fh.read(8192)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    root, outdir = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    apply = "--apply" in sys.argv[3:]
    os.makedirs(outdir, exist_ok=True)
    keep, remove = classify(root)
    rel = lambda p: os.path.relpath(p, root)
    for name, rows in (("plan_keep.tsv", keep), ("plan_remove.tsv", remove)):
        with open(os.path.join(outdir, name), "w") as out:
            out.write("path\tbytes\n")
            for p, s in sorted(rows):
                out.write(f"{rel(p)}\t{s}\n")
    kb, rb = sum(s for _, s in keep), sum(s for _, s in remove)
    print(f"plan: keep {len(keep)} files {kb/1e9:.2f} GB; remove {len(remove)} paths {rb/1e9:.2f} GB")
    if not apply:
        return
    for p, _ in remove:
        if os.path.isdir(p) and not os.path.islink(p):
            shutil.rmtree(p)
        elif os.path.lexists(p):
            os.remove(p)
    final = []
    for p, s in keep:
        if (s > COMPRESS_MIN and not p.endswith(COMPRESSED) and not os.path.islink(p)
                and is_text(p)):
            subprocess.run(["zstd", "-q", "-f", "-T4", "-19", "--rm", p], check=True)  # -f: overwrite a partial .zst from an interrupted run
            p = p + ".zst"
            subprocess.run(["zstd", "-q", "-t", p], check=True)
        final.append(p)
    with open(os.path.join(outdir, "manifest_kept.tsv"), "w") as out:
        out.write("path\tbytes\tsha256\n")
        for p in sorted(final):
            if os.path.islink(p):
                out.write(f"{rel(p)}\t0\tsymlink->{os.readlink(p)}\n")
            else:
                out.write(f"{rel(p)}\t{os.path.getsize(p)}\t{sha256(p)}\n")
    after = sum(os.path.getsize(p) for p in final if not os.path.islink(p))
    with open(os.path.join(outdir, "summary.txt"), "w") as out:
        out.write(f"root\t{root}\nplanned_keep_bytes\t{kb}\nremoved_bytes\t{rb}\n"
                  f"kept_files\t{len(final)}\nkept_bytes_after_compression\t{after}\n")
    print(f"done: kept {len(final)} files, {after/1e9:.2f} GB after compression")


if __name__ == "__main__":
    main()
