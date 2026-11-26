#!/usr/bin/env python3
"""Compress <store>/<out>/<out>.genemark.gtf to <out>.genemark.gtf.gz, in place.

Default is a dry run. Pass --apply to change files.

For each GTF the script:
  1. skips it when it is empty, when it was modified within --min-age-min minutes
     (a run may still be writing it), or when a .gz that is not older already exists;
  2. writes <gtf>.gz.tmp with pigz;
  3. inflates the temp file and compares byte count and sha256 with the original;
  4. sets the .gz mtime to the original's mtime (storedGenemarkFor() in
     nextflow/modules/funannotate/utils.nf compares mtime with the masked genome);
  5. renames the temp file to .gz, then deletes the original.
A failed check leaves the original untouched and removes the temp file.
With --remove-superseded, a plain GTF is also deleted when a .gz that passes `pigz -t` and is
not older already exists (a pipeline run publishes only the .gz, so an older plain GTF may remain).

The .mod and .other.gff3 files are not touched.

Usage:
  compress_genemark_store.py --store Fungi_BFD_runs/genemark_store            # dry run
  compress_genemark_store.py --store ... --apply --workers 8 --limit 20      # first 20
"""
import argparse
import concurrent.futures as cf
import hashlib
import os
import subprocess
import sys
import time


def sha_and_size(cmd_or_path):
    """sha256 and byte count of a file (path) or of a command's stdout (list)."""
    h = hashlib.sha256()
    n = 0
    if isinstance(cmd_or_path, list):
        p = subprocess.Popen(cmd_or_path, stdout=subprocess.PIPE)
        stream = p.stdout
    else:
        p = None
        stream = open(cmd_or_path, "rb")
    with stream:
        while True:
            b = stream.read(1 << 20)
            if not b:
                break
            h.update(b)
            n += len(b)
    if p is not None and p.wait() != 0:
        raise RuntimeError("command failed: " + " ".join(cmd_or_path))
    return h.hexdigest(), n


def convert(gtf, apply, min_age_s, remove_superseded=False):
    gz = gtf + ".gz"
    st = os.stat(gtf)
    if st.st_size == 0:
        return "skip_empty", 0, 0
    if time.time() - st.st_mtime < min_age_s:
        return "skip_recent", 0, 0
    if os.path.exists(gz) and os.stat(gz).st_mtime >= st.st_mtime and os.path.getsize(gz) > 0:
        # A newer-or-equal .gz exists (e.g. a pipeline run published it next to an older plain
        # GTF). storedGenemarkFor() already prefers the .gz, so the plain file is superseded.
        if remove_superseded:
            if subprocess.run(["pigz", "-t", gz]).returncode != 0:
                print(f"ERROR {gz}: failed pigz -t; keeping {gtf}", file=sys.stderr)
                return "error", 0, 0
            if apply:
                os.remove(gtf)
                return "removed_superseded", st.st_size, 0
            return "would_remove_superseded", st.st_size, 0
        return "skip_has_gz", 0, 0
    if not apply:
        return "would_compress", st.st_size, 0
    tmp = gz + ".tmp"
    try:
        with open(tmp, "wb") as out:
            subprocess.run(["pigz", "-p", "1", "-6", "-n", "-c", gtf], stdout=out, check=True)
        want = sha_and_size(gtf)
        got = sha_and_size(["pigz", "-dc", tmp])
        if want != got:
            raise RuntimeError(f"verify mismatch {want} != {got}")
        if os.stat(gtf).st_mtime != st.st_mtime:
            raise RuntimeError("source changed while compressing")
        os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
        os.replace(tmp, gz)
        os.remove(gtf)
        return "compressed", st.st_size, os.path.getsize(gz)
    except Exception as e:  # noqa: BLE001 - report and keep the original
        if os.path.exists(tmp):
            os.remove(tmp)
        print(f"ERROR {gtf}: {e}", file=sys.stderr)
        return "error", st.st_size, 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", required=True, help="a genemark_store directory")
    ap.add_argument("--apply", action="store_true", help="change files (default: dry run)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="stop after this many GTFs (0 = all)")
    ap.add_argument("--remove-superseded", action="store_true",
                    help="delete a plain GTF when a verified .gz that is not older already exists")
    ap.add_argument("--min-age-min", type=float, default=60, help="skip GTFs newer than this")
    a = ap.parse_args()
    if not os.path.isdir(a.store):
        sys.exit(f"ERROR: {a.store} is not a directory")

    gtfs = []
    for out in sorted(os.listdir(a.store)):
        g = os.path.join(a.store, out, out + ".genemark.gtf")
        if os.path.isfile(g):
            gtfs.append(g)
            if a.limit and len(gtfs) >= a.limit:
                break
    print(f"{len(gtfs)} GTFs found in {a.store} ({'APPLY' if a.apply else 'dry run'})", file=sys.stderr)

    counts, raw, packed = {}, 0, 0
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for status, r, p in ex.map(lambda g: convert(g, a.apply, a.min_age_min * 60, a.remove_superseded), gtfs):
            counts[status] = counts.get(status, 0) + 1
            if status in ("compressed", "would_compress", "removed_superseded", "would_remove_superseded"):
                raw += r
                packed += p
    for k, v in sorted(counts.items()):
        print(f"{k}\t{v}")
    print(f"raw_bytes\t{raw}")
    if a.apply:
        print(f"compressed_bytes\t{packed}")
    if counts.get("error"):
        sys.exit(1)


if __name__ == "__main__":
    main()
