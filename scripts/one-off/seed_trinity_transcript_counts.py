#!/usr/bin/env python3
"""
seed_trinity_transcript_counts.py — one-time backfill of
rnaseq_data/counts/<tag>.n_transcripts.txt so COUNT_TRINITY_TRANSCRIPTS's storeDir
cache (nextflow/modules/funannotate/rnaseq/COUNT_TRINITY_TRANSCRIPTS/main.nf) is
already warm the next time any do_annotation_* pipeline runs, instead of submitting
one SLURM job per species for a `grep -c '^>'` that has usually already been computed.

Why this exists
----------------
rnaseq_data/ is shared across all do_annotation_*/ run directories. Before
COUNT_TRINITY_TRANSCRIPTS had a storeDir, every pipeline pass re-ran it for every
species with reads, leaving thousands of already-computed
<tag>.n_transcripts.txt files scattered across each run's work/funannotate/*/*/
directories -- most of that work is redundant with what a fresh grep -c would
produce right now, and reusing it avoids re-reading gigabytes of Trinity FASTA.

But an old work-dir count file is only trustworthy if it is not stale relative to
the current rnaseq_data/<tag>.trinity-GG.fasta -- the recovery scripts
(pick_rnaseq_representative_override.py, fix_low_trinity.py,
one-off/reset_failed_trinity.py) delete/move that fasta and its
rnaseq_data/counts/ entry to force RNASEQ_PREPARE to rebuild it, and an old
work-dir copy from before such a rebuild would silently reintroduce the exact
staleness those scripts exist to fix.

What it does
------------
For every rnaseq_data/*.trinity-GG.fasta (skipping the
*.composite-parents.trinity-GG.fasta hybrid inputs, which COUNT_TRINITY_TRANSCRIPTS
never consumes):
  1. Skip it if rnaseq_data/counts/<tag>.n_transcripts.txt already exists AND is not
     older than the fasta (already seeded/fresh).
  2. Otherwise, look for a <tag>.n_transcripts.txt anywhere under
     <work-search-root>/do_annotation_*/work/funannotate/*/*/ and reuse the newest
     one found IF its mtime is >= the fasta's mtime (i.e. it was written from this
     exact fasta content, not a since-rebuilt one).
  3. Otherwise, compute it directly: grep -c '^>' fasta (same command
     COUNT_TRINITY_TRANSCRIPTS's script runs), falling back to 0 on grep's
     no-match exit code, exactly like the Nextflow process does.

Writes are atomic (write to a temp file, then os.replace) so a killed run can be
safely re-invoked -- it just re-does whatever wasn't finished.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

COMPOSITE_SUFFIX = ".composite-parents.trinity-GG.fasta"
FASTA_SUFFIX = ".trinity-GG.fasta"


def species_tag_from_filename(fname: str) -> str:
    return fname[: -len(FASTA_SUFFIX)]


def find_reusable_counts(work_search_root: Path) -> dict[str, Path]:
    """tag -> newest matching <tag>.n_transcripts.txt found under any
    do_annotation_*/work/funannotate/*/*/ directory."""
    best: dict[str, Path] = {}
    best_mtime: dict[str, float] = {}
    pattern = "do_annotation*/work/funannotate/*/*/*.n_transcripts.txt"
    for f in work_search_root.glob(pattern):
        tag = f.name[: -len(".n_transcripts.txt")]
        mtime = f.stat().st_mtime
        if tag not in best or mtime > best_mtime[tag]:
            best[tag] = f
            best_mtime[tag] = mtime
    return best


def grep_count(fasta_path: Path) -> int:
    """Mirror COUNT_TRINITY_TRANSCRIPTS's script exactly: grep -c '^>' fasta,
    falling back to 0 when grep exits non-zero (no matches, e.g. empty fasta)."""
    proc = subprocess.run(
        ["grep", "-c", "^>", str(fasta_path)],
        capture_output=True, text=True,
    )
    if proc.returncode == 0:
        return int(proc.stdout.strip())
    return 0


def write_count(count_path: Path, n: int, dry_run: bool) -> None:
    if dry_run:
        return
    count_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=count_path.parent, prefix=f".{count_path.name}.")
    try:
        with open(fd, "w") as fh:
            fh.write(f"{n}\n")
        Path(tmp_name).replace(count_path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def copy_count(src: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=dest.parent, prefix=f".{dest.name}.")
    try:
        with open(fd, "w") as fh:
            fh.write(src.read_text())
        Path(tmp_name).replace(dest)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--rnaseq-data", default="rnaseq_data",
                    help="Shared rnaseq_data/ dir (default: rnaseq_data)")
    ap.add_argument("--work-search-root", default=".",
                    help="Directory containing the do_annotation_*/ run dirs to scan "
                         "for reusable work-dir count files (default: current dir)")
    ap.add_argument("--project-dir", default=".")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would happen; write nothing")
    args = ap.parse_args()

    root = Path(args.project_dir).resolve()
    rnaseq_data = root / args.rnaseq_data
    work_search_root = Path(args.work_search_root).resolve()
    counts_dir = rnaseq_data / "counts"

    if not rnaseq_data.is_dir():
        sys.exit(f"[ERROR] rnaseq_data not found: {rnaseq_data}")

    print(f"[INFO] Scanning {work_search_root} for reusable work-dir count files ...",
          file=sys.stderr)
    reusable = find_reusable_counts(work_search_root)
    print(f"[INFO] Found {len(reusable)} reusable species-tagged count files "
          f"(newest per tag) under {work_search_root}", file=sys.stderr)

    fasta_files = sorted(
        f for f in rnaseq_data.glob(f"*{FASTA_SUFFIX}")
        if not f.name.endswith(COMPOSITE_SUFFIX)
    )
    print(f"[INFO] {len(fasta_files)} trinity-GG fasta files in {rnaseq_data}", file=sys.stderr)

    n_already_fresh = 0
    n_reused = 0
    n_computed = 0

    for fasta in fasta_files:
        tag = species_tag_from_filename(fasta.name)
        count_path = counts_dir / f"{tag}.n_transcripts.txt"
        fasta_mtime = fasta.stat().st_mtime

        if count_path.exists() and count_path.stat().st_mtime >= fasta_mtime:
            n_already_fresh += 1
            continue

        candidate = reusable.get(tag)
        if candidate is not None and candidate.stat().st_mtime >= fasta_mtime:
            print(f"  REUSE     {tag}  <- {candidate}")
            copy_count(candidate, count_path, args.dry_run)
            n_reused += 1
            continue

        n = grep_count(fasta)
        print(f"  COMPUTE   {tag}  n_transcripts={n}")
        write_count(count_path, n, args.dry_run)
        n_computed += 1

    mode = "[DRY-RUN] " if args.dry_run else ""
    print(f"\n{mode}Done: {n_already_fresh} already fresh, "
          f"{n_reused} reused from work dirs, {n_computed} freshly computed "
          f"(total {len(fasta_files)}).", file=sys.stderr)


if __name__ == "__main__":
    main()
