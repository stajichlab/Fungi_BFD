#!/usr/bin/env python3
"""
pick_representative_strain.py — select representative strain per species from ANI TSV
and backfill shared ab-initio parameter store.

Ported from nextflow/bin/species_reuse_clusters.py to read from the merged ANI TSV
produced by ANI_REPRESENTATIVE_SELECT (instead of DuckDB). Writes the same
abinitio_reuse_assignments.csv format so FUNANNOTATE_PREDICTION's loadAbinitioReuseMap
consumes it without modification.

BUSCO completeness and N50 (used to rank candidates in pick_representative()) are
read from tables/busco_genome.parquet and tables/asm_stats.parquet -- the same
merged, authoritative tables the rest of the pipeline queries -- rather than
re-parsing the raw BUSCO_GENOME *.BUSCO_summary.*.txt reports. Both tables are
keyed by ASMID. Note this is BUSCO_GENOME (assembly-level completeness) only;
BUSCO_PEP (annotation-level, computed after gene prediction) is a separate
downstream QC metric and is intentionally not used here.

Usage:
    pick_representative_strain.py \
        --ani-tsv all_pairs_merged.tsv \
        --predict-input predict_input_for_ani.tsv \
        --samples samples.csv \
        --tables-dir tables \
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

import duckdb

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


def load_busco_by_asmid(tables_dir: Path) -> dict:
    """Return {asmid: {complete_pct, n50_bp}} from the merged Parquet tables.

    complete_pct comes from tables/busco_genome.parquet (assembly-level BUSCO
    completeness, written by MERGE_BUSCO_GENOME from BUSCO_GENOME's per-genome
    *.BUSCO_summary.*.txt reports). n50_bp comes from tables/asm_stats.parquet
    (N50_bp, written by MERGE_ASM_STATS from AAFTF assess reports) rather than
    being re-parsed out of BUSCO's own free-text "Scaffold N50" line -- asm_stats
    is the authoritative assembly-stats table the rest of the pipeline uses.
    Both tables are keyed by ASMID. Missing/empty tables (e.g. run_busco_genome
    was false) yield an empty result, same as the old missing-directory case.
    """
    busco_path = tables_dir / "busco_genome.parquet"
    asm_stats_path = tables_dir / "asm_stats.parquet"
    if not busco_path.exists() or not asm_stats_path.exists():
        return {}
    con = duckdb.connect()
    rows = con.execute(
        """
        SELECT COALESCE(b.ASMID, a.ASMID) AS ASMID,
               CAST(COALESCE(b.complete_pct, -1.0) AS DOUBLE) AS complete_pct,
               CASE WHEN b.ASMID IS NULL THEN CAST(-1 AS BIGINT)
                    ELSE CAST(COALESCE(a.N50_bp, 0) AS BIGINT) END AS n50_bp,
               CAST(a.N50_bp AS BIGINT) AS asm_n50,
               CAST(a.total_length_bp AS BIGINT) AS total_bp,
               CAST(a.contig_count AS BIGINT) AS contigs,
               CAST(a.max_contig_bp AS BIGINT) AS max_contig_bp
        FROM read_parquet(?) AS b
        FULL OUTER JOIN read_parquet(?) AS a USING (ASMID)
        """,
        [str(busco_path), str(asm_stats_path)],
    ).fetchall()
    con.close()
    # complete_pct/n50_bp keep their old meaning (ranking key; a genome with no
    # BUSCO row ranks at the -1/-1 default, same as when it was absent from this
    # dict). asm_n50/total_bp/contigs/max_contig_bp are the raw assembly stats,
    # present even without a BUSCO row, for passes_preflight().
    return {asmid: {"complete_pct": complete_pct, "n50_bp": n50_bp,
                    "asm_n50": asm_n50, "total_bp": total_bp, "contigs": contigs,
                    "max_contig_bp": max_contig_bp}
            for asmid, complete_pct, n50_bp, asm_n50, total_bp, contigs, max_contig_bp
            in rows}


# Fallback values for running this script by hand. The pipeline
# (PICK_REPRESENTATIVE_STRAIN) always passes the predict_* params from
# nextflow.config, which FUNANNOTATE_PREDICT's pre-flight also reads.
DEFAULT_PREFLIGHT = {"min_bp": 9_000_000, "max_n50": 15_000, "max_contigs": 1000,
                     "min_contig_len": 10_000, "min_training_contigs": 1,
                     "abs_min_bp": 100_000}


def passes_preflight(stats: dict, preflight: dict) -> bool:
    """True if a genome would not be skipped by bin/asm_preflight_stats.py.

    Applies that script's three skip verdicts (too_small, no_training_contigs,
    small_fragmented) to asm_stats.parquet values; a 0 threshold disables its
    gate, as there. asm_stats has only the longest contig, so the
    no_training_contigs gate is exact for min_training_contigs == 1 and only a
    necessary condition above that. Fails open: a genome with no asm_stats row
    is treated as passing, so missing data never removes a candidate.
    preflight=None disables the check.
    """
    if not preflight or stats is None or stats.get("total_bp") is None:
        return True
    total, contigs = stats["total_bp"], stats.get("contigs")
    n50 = stats.get("asm_n50", stats.get("n50_bp"))
    if preflight["abs_min_bp"] > 0 and total < preflight["abs_min_bp"]:
        return False
    max_ctg = stats.get("max_contig_bp")
    if (preflight["min_contig_len"] > 0
            and preflight.get("min_training_contigs", 1) > 0
            and max_ctg is not None and max_ctg < preflight["min_contig_len"]):
        return False
    if preflight["min_bp"] > 0 and total < preflight["min_bp"]:
        fragmented = ((preflight["max_n50"] > 0 and n50 is not None
                       and n50 < preflight["max_n50"])
                      or (preflight["max_contigs"] > 0 and contigs is not None
                          and contigs > preflight["max_contigs"]))
        if fragmented:
            return False
    return True


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
            # Strip any directory prefix (e.g. skani's 'query/' staging dir)
            # before stripping the sequence file extension to get the asmid.
            q_asmid = _ASMID_EXT_RE.sub('', os.path.basename(q_file))
            r_asmid = _ASMID_EXT_RE.sub('', os.path.basename(r_file))
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
                         ani_pairs: dict, preflight: dict = None) -> str:
    asmid_set = set(asmids)
    ani_covered = [
        a for a in asmids
        if any((a, other) in ani_pairs for other in asmid_set if other != a)
    ]
    candidates = ani_covered if ani_covered else asmids
    # A representative that predict's pre-flight skips never produces the shared
    # ab-initio store, which blocks every sibling (allow_independent_fallback=false).
    # Prefer candidates that pass; if none do, keep the old behaviour.
    if preflight:
        predictable = [a for a in candidates
                       if passes_preflight(busco_by_asmid.get(a), preflight)]
        if predictable:
            candidates = predictable

    def sort_key(asmid):
        out = samples.get(asmid, {}).get("out", "")
        b = busco_by_asmid.get(asmid, {"complete_pct": -1.0, "n50_bp": -1})
        return (b["complete_pct"], b["n50_bp"], tuple(-ord(c) for c in out))

    return max(candidates, key=sort_key)


def compute_assignments(predict_input: dict, samples: dict, ani_pairs: dict,
                         busco_by_asmid: dict, threshold: float,
                         preflight: dict = None) -> list:
    by_species = defaultdict(list)
    for asmid, row in predict_input.items():
        sp = row.get("species", "").strip()
        if sp:
            by_species[sp].append(asmid)

    assignments = []
    for species, asmids in sorted(by_species.items()):
        if len(asmids) < 2:
            continue
        rep_asmid = pick_representative(asmids, busco_by_asmid, samples, ani_pairs,
                                        preflight)
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
        # lineterminator="\n": csv.DictWriter defaults to "\r\n" regardless
        # of platform, which -- combined with the newline="" this function
        # already sets on open() to stop that "\r\n" being doubled -- would
        # still leave every row ending in a literal "\r\n" rather than "\n".
        w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
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
    ap.add_argument("--tables-dir", required=True,
                    help="tables/ directory containing busco_genome.parquet and "
                         "asm_stats.parquet (MERGE_BUSCO_GENOME / MERGE_ASM_STATS output)")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for abinitio_reuse_assignments.csv")
    ap.add_argument("--ani-threshold", type=float, default=99.0)
    ap.add_argument("--shared-root", default="")
    ap.add_argument("--target", default="")
    ap.add_argument("--augustus-config", default="")
    ap.add_argument("--dry-run", action="store_true")
    # Pre-flight thresholds: same names and meaning as bin/asm_preflight_stats.py
    # (params.predict_* in nextflow.config). 0 disables a gate.
    ap.add_argument("--min-asm-bp", type=int, default=DEFAULT_PREFLIGHT["min_bp"],
                    help="predict_min_asm_bp (small_fragmented gate)")
    ap.add_argument("--frag-max-n50", type=int, default=DEFAULT_PREFLIGHT["max_n50"],
                    help="predict_frag_max_n50")
    ap.add_argument("--frag-max-contigs", type=int,
                    default=DEFAULT_PREFLIGHT["max_contigs"],
                    help="predict_frag_max_contigs")
    ap.add_argument("--min-contig-len", type=int,
                    default=DEFAULT_PREFLIGHT["min_contig_len"],
                    help="predict_min_training_contig_len")
    ap.add_argument("--min-training-contigs", type=int,
                    default=DEFAULT_PREFLIGHT["min_training_contigs"],
                    help="predict_min_training_contigs")
    ap.add_argument("--abs-min-bp", type=int, default=DEFAULT_PREFLIGHT["abs_min_bp"],
                    help="predict_abs_min_asm_bp (too_small gate)")
    ap.add_argument("--allow-missing-busco", action="store_true",
                    help="Proceed even if busco_genome.parquet/asm_stats.parquet are "
                         "missing or empty (representative picking then falls back to "
                         "alphabetical order for every species). Default: fail fast, "
                         "since missing tables almost always means genome_stats "
                         "(--pipeline BFD --run_busco_genome true) hasn't been merged yet.")
    args = ap.parse_args()

    predict_input = load_predict_input(Path(args.predict_input))
    samples = load_samples(Path(args.samples))
    tables_dir = Path(args.tables_dir)
    busco_by_asmid = load_busco_by_asmid(tables_dir)
    if not busco_by_asmid and not args.allow_missing_busco:
        sys.exit(
            f"[ERROR] No BUSCO/asm_stats data found in {tables_dir} "
            f"(busco_genome.parquet exists: {(tables_dir / 'busco_genome.parquet').exists()}, "
            f"asm_stats.parquet exists: {(tables_dir / 'asm_stats.parquet').exists()}). "
            "Representative selection picks by BUSCO completeness + N50 from these "
            "tables, so it must run after genome_stats (--pipeline BFD "
            "--run_busco_genome true) has been merged into tables/ — otherwise every "
            "species silently falls back to picking the alphabetically-last strain. "
            "Run genome_stats first, or pass --allow-missing-busco to proceed with "
            "alphabetical-only selection anyway."
        )
    ani_pairs = load_ani_pairs(Path(args.ani_tsv))

    print(f"[INFO] predict_input: {len(predict_input)} strains", file=sys.stderr)
    print(f"[INFO] samples: {len(samples)} entries", file=sys.stderr)
    print(f"[INFO] BUSCO summaries: {len(busco_by_asmid)}", file=sys.stderr)
    print(f"[INFO] ANI pairs: {len(ani_pairs)//2} unique pairs", file=sys.stderr)

    preflight = {"min_bp": args.min_asm_bp, "max_n50": args.frag_max_n50,
                 "max_contigs": args.frag_max_contigs,
                 "min_contig_len": args.min_contig_len,
                 "min_training_contigs": args.min_training_contigs,
                 "abs_min_bp": args.abs_min_bp}
    assignments = compute_assignments(predict_input, samples, ani_pairs,
                                       busco_by_asmid, args.ani_threshold, preflight)

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
