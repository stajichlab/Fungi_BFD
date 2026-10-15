#!/usr/bin/env python3
"""compress_predict_results.py -- report on and unify gzip state of predict_results/.

Each genome folder <root>/<genome>/predict_results/ holds a mix of plain and
gzip-compressed result files. FUNANNOTATE_PREDICT now compresses these at the
end of a run (same suffix list as gz_compress() in modules/funannotate/utils.nf).
This script brings older folders to the same state.

Default is a DRY RUN: it reads and reports, and changes nothing. Pass --execute
to change files.

Per file (plain = <name>, gz = <name>.gz), the state is one of:
  plain_only        only the plain file exists            -> compress
  gz_only           only the .gz exists                   -> nothing (already done)
  both_plain_newer  plain mtime is newer than the .gz     -> recompress, replace the .gz
  both_gz_current   .gz mtime is not older than the plain -> with --resolve-both and a
                                                             byte-identical check, remove plain
  plain_empty       plain file has size 0                 -> skipped

Compress step: pigz -c to a temp file in the same directory, pigz -t on the temp,
set the temp mtime to the plain file's mtime (ns), rename over <name>.gz, then
remove the plain file. The kept mtime matters: BFD's clearIfStale() compares result
mtimes with the input mtime, so a fresh mtime would force every result to rerun.

Safety:
  * A folder is skipped if any file in predict_results/ changed in the last
    --min-age-minutes (default 60) or an *.inflate.* temp file is present
    (funannotate update/annotate inflate .gz copies next to the .gz while running).
  * The plain file is removed only after the .gz passes pigz -t.
  * .scaffolds.fa is never compressed (funannotate update needs it plain).
  * Do not run --execute --resolve-both while update/annotate runs on those genomes.

Usage:
  compress_predict_results.py --root genome_annotation --report report.tsv        # dry run
  compress_predict_results.py --root genome_annotation --execute --jobs 8
  compress_predict_results.py --root genome_annotation --genomes list.txt --limit 20
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Same set as gz_compress() in modules/funannotate/utils.nf.
SUFFIXES = (".txt", ".mrna-transcripts.fa", ".gbk", ".gff3",
            ".cds-transcripts.fa", ".proteins.fa", ".tbl")


def is_target(name):
    return name.endswith(SUFFIXES)


def pigz_test(path):
    return subprocess.run(["pigz", "-t", path], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0


def identical(plain, gz):
    """True if gunzip(gz) is byte-identical to plain."""
    with open(plain, "rb") as pf:
        p = subprocess.Popen(["pigz", "-dc", gz], stdout=subprocess.PIPE)
        try:
            while True:
                a = p.stdout.read(1 << 20)
                b = pf.read(1 << 20)
                if a != b:
                    return False
                if not a:
                    break
        finally:
            p.stdout.close()
            p.wait()
        return p.returncode == 0


def compress(plain, gz):
    """Compress plain -> gz atomically, keep the mtime, verify. Returns True on success."""
    st = os.stat(plain)
    tmp = f"{gz}.tmp.{os.getpid()}"
    try:
        with open(tmp, "wb") as out:
            r = subprocess.run(["pigz", "-p", "1", "-c", plain], stdout=out)
        if r.returncode != 0 or not pigz_test(tmp):
            return False
        os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
        os.replace(tmp, gz)
        return True
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def process_genome(root, genome, execute, resolve_both, verify_gz, min_age_s):
    """Return a list of (genome, file, state, action, plain_bytes) rows."""
    d = os.path.join(root, genome, "predict_results")
    rows = []
    try:
        entries = list(os.scandir(d))
    except (FileNotFoundError, NotADirectoryError):
        return [(genome, "", "no_predict_results", "none", 0)]
    names = {e.name: e for e in entries}
    now = time.time()
    for e in entries:
        if ".inflate." in e.name:
            return [(genome, e.name, "skipped_inflate_tempfile", "none", 0)]
        try:
            if now - e.stat().st_mtime < min_age_s:
                return [(genome, e.name, "skipped_recent_change", "none", 0)]
        except FileNotFoundError:
            continue
    for name in sorted(names):
        if name.endswith(".gz"):
            base = name[:-3]
            if is_target(base) and base not in names:
                action = "none"
                if verify_gz and not pigz_test(os.path.join(d, name)):
                    rows.append((genome, name, "gz_only_CORRUPT", "none", 0))
                    continue
                rows.append((genome, name, "gz_only", action, 0))
            continue
        if not is_target(name):
            continue
        plain = os.path.join(d, name)
        gz = plain + ".gz"
        st = os.stat(plain)
        if st.st_size == 0:
            rows.append((genome, name, "plain_empty", "none", 0))
            continue
        if name + ".gz" not in names:
            state, action = "plain_only", "compress"
        else:
            gst = os.stat(gz)
            if st.st_mtime_ns > gst.st_mtime_ns:
                state, action = "both_plain_newer", "recompress"
            else:
                state = "both_gz_current"
                action = "remove_plain_if_identical" if resolve_both else "none"
        if execute:
            if action in ("compress", "recompress"):
                action += "_ok" if compress(plain, gz) and _remove(plain) else "_FAILED"
            elif action == "remove_plain_if_identical":
                if pigz_test(gz) and identical(plain, gz):
                    action = "removed_plain" if _remove(plain) else "remove_FAILED"
                else:
                    action = "kept_both_DIFFER_or_corrupt"
        rows.append((genome, name, state, action, st.st_size))
    if not rows:
        rows.append((genome, "", "no_result_files", "none", 0))
    return rows


def _remove(path):
    try:
        os.unlink(path)
        return True
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="genome_annotation directory")
    ap.add_argument("--genomes", help="file listing genome folder names (default: all in --root)")
    ap.add_argument("--limit", type=int, default=0, help="only the first N genomes (0 = all)")
    ap.add_argument("--execute", action="store_true", help="change files (default: dry run)")
    ap.add_argument("--resolve-both", action="store_true",
                    help="with --execute: remove a plain file when a current, identical .gz exists")
    ap.add_argument("--verify-gz", action="store_true",
                    help="run pigz -t on gz_only files (reads every .gz; slow)")
    ap.add_argument("--min-age-minutes", type=float, default=60,
                    help="skip a folder with a file changed more recently than this (default 60)")
    ap.add_argument("--jobs", type=int, default=8, help="parallel folders (default 8)")
    ap.add_argument("--report", help="write the per-file TSV here (.gz or .zst not applied)")
    args = ap.parse_args()

    if args.genomes:
        with open(args.genomes) as fh:
            genomes = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    else:
        genomes = sorted(e.name for e in os.scandir(args.root) if e.is_dir())
    if args.limit:
        genomes = genomes[:args.limit]

    mode = "EXECUTE" if args.execute else "DRY RUN (no files changed)"
    print(f"[info] {mode}: {len(genomes)} genomes under {args.root}", file=sys.stderr)

    all_rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(process_genome, args.root, g, args.execute, args.resolve_both,
                          args.verify_gz, args.min_age_minutes * 60) for g in genomes]
        for f in as_completed(futs):
            all_rows.extend(f.result())
            done += 1
            if done % 1000 == 0:
                print(f"[info] {done}/{len(genomes)} genomes", file=sys.stderr)
    all_rows.sort()

    if args.report:
        with open(args.report, "w") as out:
            out.write("genome\tfile\tstate\taction\tplain_bytes\n")
            for r in all_rows:
                out.write("\t".join(str(x) for x in r) + "\n")

    # Summary: file-state counts, then per-genome class.
    from collections import Counter, defaultdict
    states = Counter(r[2] for r in all_rows)
    actions = Counter(r[3] for r in all_rows)
    plain_bytes = sum(r[4] for r in all_rows if r[3].startswith(("compress", "recompress")))
    per = defaultdict(set)
    for g, _f, s, _a, _b in all_rows:
        per[g].add(s)
    cls = Counter()
    for g in genomes:
        s = per.get(g, {"no_predict_results"})
        if s & {"no_predict_results", "no_result_files"}:
            cls["no_predict_results_or_files"] += 1
        elif s & {"skipped_recent_change", "skipped_inflate_tempfile"}:
            cls["skipped_folder"] += 1
        elif s <= {"gz_only"}:
            cls["all_gz"] += 1
        elif "gz_only" not in s and not (s & {"both_plain_newer", "both_gz_current"}):
            cls["all_plain"] += 1
        else:
            cls["mixed"] += 1
    print("\nfile states:")
    for k, v in sorted(states.items()):
        print(f"  {k:28s}{v:>9d}")
    print("planned/performed actions:")
    for k, v in sorted(actions.items()):
        print(f"  {k:28s}{v:>9d}")
    print(f"plain bytes to compress: {plain_bytes / 1e9:.1f} GB")
    print("genome folders:")
    for k, v in sorted(cls.items()):
        print(f"  {k:28s}{v:>9d}")
    bad = sum(v for k, v in actions.items() if k.endswith("FAILED") or "DIFFER" in k) \
        + states.get("gz_only_CORRUPT", 0)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
