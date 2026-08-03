#!/usr/bin/env python3
"""Rewrite absolute cross-project symlinks in funannotate training/ dirs to local relative links.

funannotate train leaves a handful of convenience symlinks in <out>/training/:

    funannotate_train.coordSorted.bam   -> trinity.alignments.bam
    funannotate_train.transcripts.gff3  -> transcript.alignments.gff3
    funannotate_train.trinity-GG.fasta  -> trinity.fasta
    transcript.alignments.bam           -> trinity.alignments.bam

When a training directory was produced under a different project root (or moved
afterwards) those links carry the *old* absolute path, e.g.

    funannotate_train.trinity-GG.fasta -> /bigdata/.../Fungi_BFD/genome_annotation/<out>/training/trinity.fasta

which breaks as soon as the old tree goes away and makes the directory
non-relocatable. The real file almost always sits right next to the link with
the same basename, so this script relinks each one to a plain relative name.

Rule, applied per symlink:
  * relative link with no '/'          -> already local, leave alone
  * absolute (or ../-escaping) link    -> if a file named basename(target) exists
                                          in the same directory and is not the
                                          link itself, relink to that basename
  * no local counterpart               -> report, leave alone (nothing safe to do)

Safety: if the old absolute target is still readable and its size differs from
the local counterpart, the link is left alone and reported unless --force.

Default is a dry run; pass --apply to actually rewrite.

Usage:
    relink_training_symlinks.py genome_annotation_training                  # audit whole tree
    relink_training_symlinks.py --apply genome_annotation_training
    relink_training_symlinks.py --apply --quiet genome_annotation_training/<out>/training
"""

import argparse
import os
import sys

# Links funannotate creates in training/. Used only to order/limit reporting when
# --only-known is given; the generic basename rule handles everything otherwise.
KNOWN_LINKS = (
    "funannotate_train.coordSorted.bam",
    "funannotate_train.transcripts.gff3",
    "funannotate_train.trinity-GG.fasta",
    "transcript.alignments.bam",
)


def find_training_dirs(root):
    """Yield training directories at or below root.

    Accepts a single training/ dir, a single <out>/ dir, or a tree of <out>/ dirs.
    """
    root = os.path.abspath(root)
    if os.path.basename(root) == "training" and os.path.isdir(root):
        yield root
        return
    direct = os.path.join(root, "training")
    if os.path.isdir(direct):
        yield direct
        return
    try:
        entries = sorted(os.scandir(root), key=lambda e: e.name)
    except OSError as exc:
        print(f"[ERROR] cannot read {root}: {exc}", file=sys.stderr)
        return
    for entry in entries:
        if not entry.is_dir(follow_symlinks=False):
            continue
        sub = os.path.join(entry.path, "training")
        if os.path.isdir(sub):
            yield sub


def fix_dir(traindir, apply_changes=False, force=False, only_known=False, quiet=False):
    """Relink absolute symlinks in one training dir. Returns (fixed, skipped, orphan)."""
    fixed = skipped = orphan = 0
    try:
        entries = sorted(os.scandir(traindir), key=lambda e: e.name)
    except OSError as exc:
        print(f"[ERROR] cannot read {traindir}: {exc}", file=sys.stderr)
        return 0, 0, 0

    for entry in entries:
        if not entry.is_symlink():
            continue
        if only_known and entry.name not in KNOWN_LINKS:
            continue
        link = entry.path
        target = os.readlink(link)
        if "/" not in target:
            continue  # already a bare local name
        base = os.path.basename(target)
        if not base:
            continue
        if base == entry.name:
            # Relinking would make the link point at itself.
            print(f"[SKIP] {link}: target basename equals link name ({base})")
            skipped += 1
            continue
        local = os.path.join(traindir, base)
        if not os.path.exists(local):
            print(f"[ORPHAN] {link} -> {target} (no local {base})")
            orphan += 1
            continue
        # Size sanity check against the old target, when it is still reachable.
        if not force:
            try:
                old_size = os.stat(target).st_size
                new_size = os.stat(local).st_size
                if old_size != new_size:
                    print(
                        f"[SKIP] {link}: old target {old_size} bytes != local {base} "
                        f"{new_size} bytes (use --force to relink anyway)"
                    )
                    skipped += 1
                    continue
            except OSError:
                pass  # old target gone/unreadable — local copy is all we have
        if apply_changes:
            tmp = link + ".relink.tmp"
            try:
                if os.path.lexists(tmp):
                    os.remove(tmp)
                os.symlink(base, tmp)
                os.replace(tmp, link)  # atomic swap; never leaves the name missing
            except OSError as exc:
                print(f"[ERROR] {link}: {exc}", file=sys.stderr)
                if os.path.lexists(tmp):
                    os.remove(tmp)
                skipped += 1
                continue
            if not quiet:
                print(f"[FIXED] {link} -> {base}")
        else:
            print(f"[DRY-RUN] {link}: {target} -> {base}")
        fixed += 1
    return fixed, skipped, orphan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+",
                    help="training dir, <out> dir, or a tree containing <out>/training dirs")
    ap.add_argument("--apply", action="store_true",
                    help="actually rewrite the symlinks (default is a dry run)")
    ap.add_argument("--force", action="store_true",
                    help="relink even when the old target's size differs from the local file")
    ap.add_argument("--only-known", action="store_true",
                    help="only touch the four symlinks funannotate train creates")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress per-link FIXED lines (problems are still printed)")
    args = ap.parse_args()

    tot_fixed = tot_skipped = tot_orphan = tot_dirs = 0
    for path in args.paths:
        for traindir in find_training_dirs(path):
            tot_dirs += 1
            f, s, o = fix_dir(traindir, apply_changes=args.apply, force=args.force,
                              only_known=args.only_known, quiet=args.quiet)
            tot_fixed += f
            tot_skipped += s
            tot_orphan += o

    verb = "relinked" if args.apply else "would relink"
    print(f"[SUMMARY] {tot_dirs} training dirs: {verb} {tot_fixed}, "
          f"skipped {tot_skipped}, orphaned {tot_orphan}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
