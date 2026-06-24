#!/usr/bin/env python3
"""Reclaim disk by stubbing normalized RNA-seq reads whose work is fully done.

For each species_tag that satisfies the SAFE predicate (see audit_rnaseq_reclaim.py):
  1. Trinity FASTA exists and is non-empty in rnaseq_data/, AND
  2. EVERY assembly mapped to that species_tag has a non-empty pasa.gff3.

the script:
  a. records a provenance manifest rnaseq_data/<species_tag>.reads_reclaimed.json
     (sha256 + bytes + mtime of each read file, Trinity sha256, trained strains,
      and the SRA accessions used, if a sra_query CSV is present), THEN
  b. replaces rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz with zero-byte
     stubs (default mode), preserving each file's ORIGINAL mtime so the empty stub
     can never be "newer than the GBK" and trip staleRnaseq()/the train guard.

Why zero-byte stubs and not rm: SRA_FETCH declares `storeDir rnaseq_reads/` and is
skipped only when its declared output files exist. Deleting them re-triggers an
expensive re-download; an empty stub keeps storeDir satisfied and every downstream
branch already gates on `.size() > 0`, so a stub reads as "no reads" without churn.

DRY-RUN BY DEFAULT. Pass --apply to act.

Caveat (future strains): if a NEW strain is later added to samples.csv for an
already-reclaimed species, its FUNANNOTATE_TRAIN would receive an empty read +
valid Trinity and hit funannotate.nf's "PASA only, no reads" branch (empty
--left_norm). To genuinely refresh such a species, restore real reads first
(delete Trinity + training and re-run, which re-fetches). The manifest's
"strains_trained" list lets a pre-flight check detect this.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_rnaseq_reclaim import (  # noqa: E402  reuse the audited helpers
    clean_strain, make_sample_tag, species_tag_of, first_existing, human,
)

READ_SUFFIXES = ("_norm_R1.fastq.gz", "_norm_R2.fastq.gz", "_norm_SE.fastq.gz")


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def iso(ts):
    return dt.datetime.fromtimestamp(ts).astimezone().isoformat(timespec="seconds")


def load_species(samples):
    """species_tag -> list of (out, asmid)."""
    species = {}
    with open(samples, newline="") as fh:
        for row in csv.DictReader(fh):
            sp_in = row.get("SPECIES") or ""
            strain = row.get("STRAIN") or ""
            asmid = (row.get("ASMID") or "").strip()
            if not sp_in.strip() or not asmid:
                continue
            stag = species_tag_of(sp_in)
            out = make_sample_tag(sp_in, strain)
            species.setdefault(stag, []).append((out, asmid))
    return species


def read_accessions(sra_query_csv):
    """Best-effort list of accessions used; returns [] if unavailable."""
    if not sra_query_csv.exists():
        return []
    accs = []
    try:
        with open(sra_query_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            # column name varies; grab the first that looks like a run accession
            acc_col = None
            for c in (reader.fieldnames or []):
                if c and c.lower() in ("run", "run_accession", "accession", "acc"):
                    acc_col = c
                    break
            for row in reader:
                if acc_col:
                    v = (row.get(acc_col) or "").strip()
                    if v:
                        accs.append(v)
    except Exception:
        return []
    return accs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="launchDir (default: cwd)")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--reads-dir", default="rnaseq_reads")
    ap.add_argument("--rnaseq-data", default="rnaseq_data")
    ap.add_argument("--training", default="genome_annotation_training")
    ap.add_argument("--mode", choices=("stub", "delete"), default="stub",
                    help="stub: replace with zero-byte file (default, storeDir-safe); "
                         "delete: rm the file (WILL re-trigger SRA_FETCH download)")
    ap.add_argument("--only", action="append", default=[],
                    help="limit to these species_tag(s); repeatable")
    ap.add_argument("--force", action="store_true",
                    help="re-process species that already have a manifest")
    ap.add_argument("--apply", action="store_true",
                    help="actually write manifests and stub/delete (default: dry-run)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    samples = root / args.samples
    reads_dir = root / args.reads_dir
    data_dir = root / args.rnaseq_data
    train_dir = root / args.training
    query_dir = reads_dir / "sra_query"

    if not samples.exists():
        sys.exit(f"ERROR: samples file not found: {samples}")

    species = load_species(samples)
    only = set(args.only)

    planned, skipped_blocked, reclaimed_bytes = [], 0, 0

    for stag, members in sorted(species.items()):
        if only and stag not in only:
            continue

        trinity = data_dir / f"{stag}.trinity-GG.fasta"
        if not trinity.exists() or trinity.stat().st_size == 0:
            continue  # NO_TRINITY / EMPTY_TRINITY -> nothing to reclaim

        # Re-verify the full predicate at action time (never trust a stale audit).
        untrained = []
        for out, _asmid in members:
            gff = first_existing(
                train_dir / out / "training" / "funannotate_train.pasa.gff3",
                train_dir / out / "training" / "funannotate_train.pasa.gff3.gz",
            )
            if gff is None or gff.stat().st_size == 0:
                untrained.append(out)
        if untrained:
            skipped_blocked += 1
            continue

        read_files = [reads_dir / f"{stag}{suf}" for suf in READ_SUFFIXES]
        present = [p for p in read_files if p.exists() and p.stat().st_size > 0]
        if not present:
            continue  # already reclaimed / empty stubs

        manifest_path = data_dir / f"{stag}.reads_reclaimed.json"
        if manifest_path.exists() and not args.force:
            continue

        species_bytes = sum(p.stat().st_size for p in present)
        reclaimed_bytes += species_bytes
        planned.append((stag, present, read_files, trinity, members, manifest_path,
                        species_bytes))

    # ── Report plan ────────────────────────────────────────────────────────────
    verb = "RECLAIM" if args.apply else "WOULD RECLAIM"
    print(f"\n{verb} ({args.mode} mode) — root={root}")
    print(f"  species eligible now : {len(planned)}")
    print(f"  blocked (untrained)  : {skipped_blocked}")
    print(f"  reclaimable          : {human(reclaimed_bytes)}")
    print("-" * 64)

    if not args.apply:
        for stag, present, _all, _t, _m, _mp, b in sorted(planned, key=lambda x: -x[-1]):
            print(f"  {human(b):>10}  {stag}  ({len(present)} file(s))")
        print("\n(dry-run; re-run with --apply to act)")
        return

    # ── Act ──────────────────────────────────────────────────────────────────
    for stag, present, all_files, trinity, members, manifest_path, b in planned:
        manifest = {
            "species_tag": stag,
            "reclaimed_at": iso(dt.datetime.now().timestamp()),
            "mode": args.mode,
            "trinity_fasta": trinity.name,
            "trinity_sha256": sha256_of(trinity),
            "trinity_bytes": trinity.stat().st_size,
            "trinity_mtime": iso(trinity.stat().st_mtime),
            "strains_trained": sorted(out for out, _ in members),
            "accessions": read_accessions(query_dir / f"{stag}.sra_query.csv"),
            "reads": {},
        }
        for p in all_files:
            if not p.exists():
                continue
            st = p.stat()
            key = p.name.split("_norm_")[-1].split(".")[0]  # R1 / R2 / SE
            entry = {"bytes": st.st_size, "mtime": iso(st.st_mtime)}
            if st.st_size > 0:
                entry["sha256"] = sha256_of(p)
            manifest["reads"][key] = entry

        # Write manifest FIRST so provenance survives even if stubbing is interrupted.
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

        for p in present:
            orig_mtime = p.stat().st_mtime
            if args.mode == "delete":
                p.unlink()
            else:  # stub: truncate to zero, restore original mtime (never newer than GBK)
                with open(p, "wb"):
                    pass
                os.utime(p, (orig_mtime, orig_mtime))
        print(f"  done  {human(b):>10}  {stag}  -> manifest {manifest_path.name}")

    print(f"\nReclaimed {human(reclaimed_bytes)} across {len(planned)} species.")


if __name__ == "__main__":
    main()
