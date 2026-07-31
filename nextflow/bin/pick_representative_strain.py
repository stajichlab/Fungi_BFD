#!/usr/bin/env python3
"""
pick_representative_strain.py — select representative strain per species from ANI TSV
and backfill shared ab-initio parameter store.

Ported from nextflow/bin/species_reuse_clusters.py to read from the merged ANI TSV
produced by ANI_REPRESENTATIVE_SELECT (instead of DuckDB). Writes the same
abinitio_reuse_assignments.csv format so FUNANNOTATE_PREDICTION's loadAbinitioReuseMap
consumes it without modification.

Usage:
    pick_representative_strain.py \
        --ani-tsv all_pairs_merged.tsv \
        --predict-input predict_input_for_ani.tsv \
        --samples samples.csv \
        --busco-dir results/genome_stats/BUSCO_genome \
        --out-dir genome_annotation/_reuse_assignments \
        --ani-threshold 99.0 \
        [--shared-root gene_prediction_shared_abinitio] \
        [--target genome_annotation] \
        [--augustus-config lib/augustus/3.5/config]
"""

import argparse
import csv
import fcntl
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

from backfill_abinitio_params import BackfillOutcome, backfill_species_store

_ASMID_EXT_RE = re.compile(r'(\.fasta|\.fna|\.fa|\.faa|\.fas)(\.gz)?$')


def clean_strain(raw_strain: str) -> str:
    s = (raw_strain or "").strip()
    s = re.sub(r"['\"]", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    sp = re.sub(r"['\"]", "", (raw_species or "").strip())
    st = clean_strain(raw_strain)
    parts = [p for p in (sp, st) if p]
    tag = "_".join(parts)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def load_samples(samples_path: Path) -> dict:
    rows = {}
    with open(samples_path, newline="") as fh:
        for row in csv.DictReader(fh):
            asmid = (row.get("ASMID") or "").strip()
            species = (row.get("SPECIES") or "").strip().strip("'\"")
            strain = (row.get("STRAIN") or "").strip()
            if not asmid or not species:
                continue
            out = make_sample_tag(species, strain)
            species_tag = re.sub(r"\s+", "_", species)
            rows[asmid] = {
                "species": species, "strain": strain, "out": out,
                "species_tag": species_tag,
            }
    return rows


_COMPLETE_RE = re.compile(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)")
_N50_RE = re.compile(r"([\d.]+)\s*(Mbp|kbp|bp)\s*Scaffold N50", re.IGNORECASE)
_UNIT_MULT = {"bp": 1, "kbp": 1_000, "Mbp": 1_000_000}


def parse_busco_genome_summary(path: Path):
    text = path.read_text(errors="replace")
    m = _COMPLETE_RE.search(text)
    if not m:
        return None
    n50_bp = 0
    mn = _N50_RE.search(text)
    if mn:
        n50_bp = int(float(mn.group(1)) * _UNIT_MULT[mn.group(2)])
    return {"complete_pct": float(m[1]), "n50_bp": n50_bp}


def load_busco_by_asmid(busco_dir: Path) -> dict:
    """Return {asmid: {complete_pct, n50_bp}} from results/genome_stats/BUSCO_genome/.

    BUSCO_genome is ASMID-keyed (assembly-level statistic) and hash-bucketed
    (nextflow/modules/BFD/BUSCO_GENOME/main.nf writes
    <bucket>/<asmid>.BUSCO_summary.<lineage>.txt) -- the `*/*` glob covers the
    one bucket-subdirectory level; the ASMID is still the filename's own
    prefix before ".BUSCO_summary.", same extraction as before the reorg, just
    keyed by ASMID now instead of the old Genus_species_strain basename.
    """
    result = {}
    if not busco_dir.is_dir():
        return result
    for p in busco_dir.glob("*/*.BUSCO_summary.*.txt"):
        asmid = p.name.split(".BUSCO_summary.")[0]
        parsed = parse_busco_genome_summary(p)
        if parsed:
            result[asmid] = parsed
    return result


def load_predict_input(path: Path) -> dict:
    """Return {asmid: {out, species, strain, ...}} from the predict_input TSV."""
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asmid = row.get("asmid", "").strip()
            if asmid:
                result[asmid] = dict(row)
    return result


def load_ani_pairs(ani_tsv_path: Path) -> dict:
    """Return {(q_asmid, r_asmid): ani} for every pair in the merged TSV.
    Infers asmid from filename stem (everything before the first '.').
    Handles both 3-column (q,r,ani) and 5-column (q,r,ani,af_r,af_q) formats."""
    pairs = {}
    if not ani_tsv_path.exists() or ani_tsv_path.name == '/dev/null':
        return pairs
    with open(ani_tsv_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return pairs
        for row in reader:
            if len(row) < 3:
                continue
            q_file = row[0].strip()
            r_file = row[1].strip()
            try:
                ani = float(row[2])
            except (ValueError, IndexError):
                continue
            # Strip common sequence file extensions to get the asmid
            q_asmid = _ASMID_EXT_RE.sub('', q_file)
            r_asmid = _ASMID_EXT_RE.sub('', r_file)
            if q_asmid and r_asmid:
                pairs[(q_asmid, r_asmid)] = ani
                pairs.setdefault((r_asmid, q_asmid), ani)
    return pairs


def _sparse_ani_warning(species: str, asmids: list, ani_pairs: dict,
                         rep_asmid: str) -> None:
    """Warn if ANI TSV appears sparse — some intra-species pairs are missing."""
    n = len(asmids)
    expected = n * (n - 1)  # (q,r) + (r,q) stored for every comparison
    asmid_set = set(asmids)
    found = sum(
        1 for (q, r) in ani_pairs
        if q in asmid_set and r in asmid_set and q != r
    )
    missing = expected - found
    if missing > 0:
        print(f"[WARN] {species}: ANI TSV has only {found}/{expected} intra-species "
              f"pairs ({missing} missing). This typically means a sparse ANI method "
              f"(skani/mash/sourmash) missed some comparisons. Strain {rep_asmid} "
              f"was selected as representative on BUSCO+N50 alone.",
              file=sys.stderr)


def pick_representative(asmids: list, busco_by_asmid: dict, samples: dict,
                         ani_pairs: dict) -> str:
    asmid_set = set(asmids)
    ani_covered = [
        a for a in asmids
        if any((a, other) in ani_pairs for other in asmid_set if other != a)
    ]
    candidates = ani_covered if ani_covered else asmids

    def sort_key(asmid):
        out = samples.get(asmid, {}).get("out", "")
        b = busco_by_asmid.get(asmid, {"complete_pct": -1.0, "n50_bp": -1})
        return (b["complete_pct"], b["n50_bp"], tuple(-ord(c) for c in out))

    return max(candidates, key=sort_key)


def compute_assignments(predict_input: dict, samples: dict, ani_pairs: dict,
                         busco_by_asmid: dict, threshold: float) -> list:
    by_species = defaultdict(list)
    for asmid, row in predict_input.items():
        sp = row.get("species", "").strip()
        if sp:
            by_species[sp].append(asmid)

    assignments = []
    for species, asmids in sorted(by_species.items()):
        if len(asmids) < 2:
            continue
        rep_asmid = pick_representative(asmids, busco_by_asmid, samples, ani_pairs)
        rep_out = samples.get(rep_asmid, {}).get("out", "")
        _sparse_ani_warning(species, asmids, ani_pairs, rep_asmid)

        for asmid in sorted(asmids, key=lambda a: samples.get(a, {}).get("out", "")):
            out = samples.get(asmid, {}).get("out", "")
            is_rep = (asmid == rep_asmid)
            if is_rep:
                assignments.append({
                    "species": species, "out": out,
                    "is_representative": True,
                    "representative_out": rep_out,
                    "ani_to_representative": 100.0,
                    "reuse_eligible": False,
                })
                continue
            key = (asmid, rep_asmid)
            ani = ani_pairs.get(key)
            eligible = ani is not None and ani >= threshold
            assignments.append({
                "species": species, "out": out,
                "is_representative": False,
                "representative_out": rep_out,
                "ani_to_representative": ani if ani is not None else "",
                "reuse_eligible": eligible,
            })
    return assignments


def _atomic_write_csv(path: Path, fieldnames: list, rows: list):
    """Write rows to path via a same-directory temp file + os.replace(), so any
    concurrent reader (loadAbinitioReuseMap(), a running funannotate pipeline)
    always sees either the fully-old or fully-new file, never a torn one --
    os.replace() is a single atomic rename on POSIX for a plain file (unlike a
    directory replace, no old-file-must-move-aside dance needed here)."""
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(tmp_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    os.replace(tmp_path, path)


def write_assignments(assignments: list, out_dir: Path):
    """Write per-species abinitio_reuse_assignments.{species}.csv files and
    merge this run's rows into the combined abinitio_reuse_assignments.csv
    that loadAbinitioReuseMap() reads.

    The combined file is a merge, not an overwrite: rows for species this run
    didn't touch (out of scope for the current --taxon/--compare/samples
    selection) are preserved from whatever's already on disk, and only rows
    for species this run *did* recompute are replaced. Without this, running
    compare_ani in separate scoped batches (e.g. one phylum at a time) would
    have each run silently delete every other phylum's rows from the combined
    file, since it never read the existing content before overwriting it.

    A flock() around the read-merge-write critical section serializes
    concurrent writers (e.g. two compare_ani runs for different phyla
    finishing around the same time): os.replace() alone guarantees no reader
    ever sees a torn file, but it does NOT prevent a lost update if two
    writers both read the same stale baseline before either one's merge
    lands -- the lock closes that gap. The lock file itself is never deleted
    (removing it would let a writer that already opened it race a new writer
    that creates a fresh inode of the same name, defeating the lock).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    by_species = defaultdict(list)
    for row in assignments:
        by_species[row["species"]].append(row)

    fieldnames = ["species", "out", "is_representative", "representative_out",
                  "ani_to_representative", "reuse_eligible"]

    for species, rows in by_species.items():
        species_tag = re.sub(r"\s+", "_", species)
        out_path = out_dir / f"abinitio_reuse_assignments.{species_tag}.csv"
        _atomic_write_csv(out_path, fieldnames, rows)

    combined_path = out_dir / "abinitio_reuse_assignments.csv"
    lock_path = out_dir / ".abinitio_reuse_assignments.csv.lock"

    with open(lock_path, "a") as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
        try:
            species_this_run = set(by_species.keys())
            existing_rows = []
            if combined_path.exists():
                with open(combined_path, newline="") as fh:
                    existing_rows = list(csv.DictReader(fh))
            kept = [r for r in existing_rows if r.get("species") not in species_this_run]
            merged = kept + assignments
            merged.sort(key=lambda r: (r.get("species", ""), r.get("out", "")))
            _atomic_write_csv(combined_path, fieldnames, merged)
        finally:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)


def write_repr_assignments(assignments: list, out_path: Path):
    """Write a simple TSV for pipeline channel joins:
    out, species, is_representative, representative_out,
    ani_to_representative, reuse_eligible"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("out\tspecies\tis_representative\t"
                 "representative_out\tani_to_representative\treuse_eligible\n")
        for a in assignments:
            fh.write(f"{a['out']}\t{a['species']}\t"
                     f"{a['is_representative']}\t{a['representative_out']}\t"
                     f"{a['ani_to_representative']}\t{a['reuse_eligible']}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ani-tsv", required=True,
                    help="Merged all-pairs TSV from CONCAT_ANI_TSVS (may be /dev/null)")
    ap.add_argument("--predict-input", required=True,
                    help="predict_input_for_ani.tsv written by ANI_REPRESENTATIVE_SELECT")
    ap.add_argument("--samples", required=True, help="samples.csv")
    ap.add_argument("--busco-dir", required=True,
                    help="results/genome_stats/BUSCO_genome/")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for abinitio_reuse_assignments.csv")
    ap.add_argument("--ani-threshold", type=float, default=99.0)
    ap.add_argument("--shared-root", default="")
    ap.add_argument("--target", default="")
    ap.add_argument("--augustus-config", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--allow-missing-busco", action="store_true",
                    help="Proceed even if --busco-dir is missing or has no parseable "
                         "summaries (representative picking then falls back to "
                         "alphabetical order for every species). Default: fail fast, "
                         "since an empty/missing dir almost always means genome_stats "
                         "(--pipeline BFD --run_busco_genome) hasn't been run yet.")
    args = ap.parse_args()

    predict_input = load_predict_input(Path(args.predict_input))
    samples = load_samples(Path(args.samples))
    busco_dir = Path(args.busco_dir)
    busco_by_asmid = load_busco_by_asmid(busco_dir)
    if not busco_by_asmid and not args.allow_missing_busco:
        sys.exit(
            f"[ERROR] No BUSCO summaries found in {busco_dir} "
            f"(exists: {busco_dir.is_dir()}). Representative selection picks by BUSCO "
            "completeness + N50 parsed from these files, so it must run after "
            "genome_stats (--pipeline BFD --run_busco_genome true) has published its "
            "output here — otherwise every species silently falls back to picking "
            "the alphabetically-last strain. Run genome_stats first, or pass "
            "--allow-missing-busco to proceed with alphabetical-only selection anyway."
        )
    ani_pairs = load_ani_pairs(Path(args.ani_tsv))

    print(f"[INFO] predict_input: {len(predict_input)} strains", file=sys.stderr)
    print(f"[INFO] samples: {len(samples)} entries", file=sys.stderr)
    print(f"[INFO] BUSCO summaries: {len(busco_by_asmid)}", file=sys.stderr)
    print(f"[INFO] ANI pairs: {len(ani_pairs)//2} unique pairs", file=sys.stderr)

    assignments = compute_assignments(predict_input, samples, ani_pairs,
                                       busco_by_asmid, args.ani_threshold)

    out_dir = Path(args.out_dir)
    write_assignments(assignments, out_dir)

    n_species = len({a["species"] for a in assignments})
    n_eligible = sum(1 for a in assignments if a["reuse_eligible"])
    n_rep = sum(1 for a in assignments if a["is_representative"])
    n_no_ani = sum(1 for a in assignments
                   if not a["is_representative"] and a["ani_to_representative"] == "")
    print(f"[INFO] Wrote {len(assignments)} rows ({n_species} species, "
          f"{n_rep} representative, {n_eligible} reuse_eligible, "
          f"{n_no_ani} fail-closed/missing-ANI) to {out_dir}/", file=sys.stderr)

    # Write repr_assignments.tsv for pipeline channel joins.
    repr_path = Path(args.out_dir) / "repr_assignments.tsv"
    write_repr_assignments(assignments, repr_path)

    if args.dry_run:
        return

    if args.shared_root:
        target = Path(args.target) if args.target else Path.cwd()
        shared_root = Path(args.shared_root)
        aug_cfg = Path(args.augustus_config) if args.augustus_config else None
        reps = {a["species"]: a["representative_out"]
                for a in assignments if a["is_representative"]}
        n_backfilled = 0
        for species, rep_out in sorted(reps.items()):
            result = backfill_species_store(species, rep_out, target, shared_root,
                                             args.ani_threshold,
                                             aug_cfg=aug_cfg, dry_run=args.dry_run)
            if result.outcome in (BackfillOutcome.BACKFILLED, BackfillOutcome.UP_TO_DATE):
                n_backfilled += 1
        print(f"[INFO] Backfilled {n_backfilled}/{len(reps)} species stores "
              "(includes already-up-to-date)", file=sys.stderr)


if __name__ == "__main__":
    main()
