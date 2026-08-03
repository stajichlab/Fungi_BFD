#!/usr/bin/env python3
"""
profile_train_pasa_trinity.py — profile PASA/Trinity step runtimes across
funannotate-train.log files to see whether tool/version upgrades over time
are changing performance.

Scans <root>/*/logfiles/funannotate-train.log, parses the timestamped log
lines, and extracts per-run:
  - run date, funannotate version, Trinity version, PASA version
  - trinity_prep_seconds : poly-A cleaning (seqclean) + minimap2 alignment of
                            the Trinity assembly to the genome, up to the
                            start of the PASA alignment step
  - pasa_align_seconds   : Launch_PASA_pipeline.pl step (transcripts -> loci)
  - pasa_total_seconds   : PASA alignment step through "PASA finished"
  - total_runtime_seconds: first log line to last log line
  - status               : completed / error / incomplete
  - *_cached flags        : True if funannotate reused a prior result instead
                            of re-running that stage (duration is not
                            meaningful in that case)

Writes a CSV with one row per run and prints a summary of median step
durations grouped by funannotate version (and by PASA/Trinity version, since
those track funannotate's bundled toolchain here).
"""

import argparse
import csv
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

TS_RE = re.compile(r"^\[(\d{2}/\d{2}/\d{2}) (\d{2}:\d{2}:\d{2})\]:\s?(.*)$")
FUNANNOTATE_VERSION_RE = re.compile(r"^Running (\d\S*)$")
TRINITY_VERSION_RE = re.compile(r"Trinity version=(\S+)")
PASA_VERSION_RE = re.compile(r"PASA/(\d+\.\d+(?:\.\d+)?)")
PASA_ALIGN_START_RE = re.compile(r"^Running PASA alignment step using ([\d,]+) transcripts")
PASA_ASSIGNED_RE = re.compile(
    r"^PASA assigned ([\d,]+) transcripts to ([\d,]+) loci"
)
PASA_FINISHED_RE = re.compile(r"^PASA finished\.")
POLYA_START_RE = re.compile(r"^Removing poly-A sequences")
TRAIN_COMPLETE_RE = re.compile(r"^Trinity/PASA has completed")
CMD_ERROR_RE = re.compile(r"^(CMD ERROR|Error running|Traceback)")


def parse_timestamp(datestr, timestr):
    return datetime.strptime(f"{datestr} {timestr}", "%m/%d/%y %H:%M:%S")


def parse_log(path):
    """Return a dict of profiled fields for one funannotate-train.log, or None
    if the file has no recognizable timestamped content."""
    raw_text = path.read_text(errors="replace")

    lines = []
    for raw in raw_text.splitlines():
        m = TS_RE.match(raw)
        if m:
            dt = parse_timestamp(m.group(1), m.group(2))
            lines.append((dt, m.group(3)))

    if not lines:
        return None

    result = {
        "run_dir": path.parent.parent.name,
        "log_path": str(path),
        "run_date": lines[0][0].isoformat(sep=" "),
        "funannotate_version": None,
        "trinity_version": None,
        "pasa_version": None,
        "n_transcripts_input": None,
        "n_transcripts_assigned": None,
        "n_loci": None,
        "trinity_prep_seconds": None,
        "trinity_prep_cached": False,
        "pasa_align_seconds": None,
        "pasa_align_cached": False,
        "pasa_total_seconds": None,
        "total_runtime_seconds": (lines[-1][0] - lines[0][0]).total_seconds(),
        "status": "incomplete",
    }

    # PASA version can appear in untimestamped subprocess-echoed lines
    # (e.g. seqclean's psx invocation), so search the raw file, not just
    # the timestamped log lines.
    pasa_ver_m = PASA_VERSION_RE.search(raw_text)
    if pasa_ver_m:
        result["pasa_version"] = pasa_ver_m.group(1)

    polya_start = None
    pasa_align_start = None
    pasa_assigned_time = None
    pasa_finished_time = None

    for dt, msg in lines:
        if result["funannotate_version"] is None:
            fv = FUNANNOTATE_VERSION_RE.match(msg)
            if fv:
                result["funannotate_version"] = fv.group(1)
        if result["trinity_version"] is None:
            tv = TRINITY_VERSION_RE.search(msg)
            if tv:
                result["trinity_version"] = tv.group(1)

        if "Existing SeqClean output found" in msg or "Existing BAM alignments found" in msg:
            result["trinity_prep_cached"] = True
        if POLYA_START_RE.match(msg) and polya_start is None:
            polya_start = dt

        if "Existing PASA assemblies found" in msg:
            result["pasa_align_cached"] = True
        pam = PASA_ALIGN_START_RE.match(msg)
        if pam and pasa_align_start is None:
            pasa_align_start = dt
            result["n_transcripts_input"] = int(pam.group(1).replace(",", ""))

        pasm = PASA_ASSIGNED_RE.match(msg)
        if pasm and pasa_assigned_time is None:
            pasa_assigned_time = dt
            result["n_transcripts_assigned"] = int(pasm.group(1).replace(",", ""))
            result["n_loci"] = int(pasm.group(2).replace(",", ""))

        if PASA_FINISHED_RE.match(msg) and pasa_finished_time is None:
            pasa_finished_time = dt

        if TRAIN_COMPLETE_RE.match(msg):
            result["status"] = "completed"
        if CMD_ERROR_RE.match(msg) and result["status"] != "completed":
            result["status"] = "error"

    # trinity prep duration: poly-A start -> PASA alignment step start
    if polya_start is not None and not result["trinity_prep_cached"]:
        end = pasa_align_start or pasa_assigned_time
        if end is not None:
            result["trinity_prep_seconds"] = (end - polya_start).total_seconds()

    # PASA alignment duration: PASA alignment start -> transcripts assigned to loci
    if pasa_align_start is not None and pasa_assigned_time is not None and not result["pasa_align_cached"]:
        result["pasa_align_seconds"] = (pasa_assigned_time - pasa_align_start).total_seconds()

    # PASA total duration: PASA alignment start -> "PASA finished"
    if pasa_align_start is not None and pasa_finished_time is not None and not result["pasa_align_cached"]:
        result["pasa_total_seconds"] = (pasa_finished_time - pasa_align_start).total_seconds()

    return result


def find_logs(root):
    return sorted(Path(root).glob("*/logfiles/funannotate-train.log"))


def write_csv(rows, out_path):
    fieldnames = [
        "run_dir", "run_date", "status",
        "funannotate_version", "trinity_version", "pasa_version",
        "n_transcripts_input", "n_transcripts_assigned", "n_loci",
        "trinity_prep_seconds", "trinity_prep_cached",
        "pasa_align_seconds", "pasa_align_cached",
        "pasa_total_seconds", "total_runtime_seconds",
        "log_path",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fieldnames})


def summarize(rows, field, group_field="funannotate_version"):
    groups = {}
    for row in rows:
        val = row.get(field)
        if val is None:
            continue
        groups.setdefault(row.get(group_field), []).append(val)

    print(f"\n== median {field} by {group_field} ==")
    for group, vals in sorted(groups.items(), key=lambda kv: (kv[0] is None, kv[0])):
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        med = statistics.median(vals_sorted)
        print(f"  {group!s:30s} n={n:4d}  median={med:9.1f}s  "
              f"min={vals_sorted[0]:9.1f}s  max={vals_sorted[-1]:9.1f}s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="genome_annotation_training",
                     help="Directory containing per-strain run dirs "
                          "(each with logfiles/funannotate-train.log)")
    ap.add_argument("--out", default="analysis/funannotate_train_pasa_trinity_profile.csv",
                     help="Output CSV path (one row per run)")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        sys.exit(f"root directory not found: {root}")

    logs = find_logs(root)
    if not logs:
        sys.exit(f"no funannotate-train.log files found under {root}/*/logfiles/")

    rows = []
    for log_path in logs:
        parsed = parse_log(log_path)
        if parsed is not None:
            rows.append(parsed)

    rows.sort(key=lambda r: r["run_date"])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(rows, out_path)

    n_completed = sum(1 for r in rows if r["status"] == "completed")
    n_error = sum(1 for r in rows if r["status"] == "error")
    n_incomplete = sum(1 for r in rows if r["status"] == "incomplete")
    print(f"Parsed {len(rows)} logs -> {out_path}")
    print(f"  completed={n_completed}  error={n_error}  incomplete={n_incomplete}")

    versions = sorted({r["funannotate_version"] for r in rows if r["funannotate_version"]})
    print(f"\nfunannotate versions seen: {versions}")
    trinity_versions = sorted({r["trinity_version"] for r in rows if r["trinity_version"]})
    print(f"Trinity versions seen: {trinity_versions}")
    pasa_versions = sorted({r["pasa_version"] for r in rows if r["pasa_version"]})
    print(f"PASA versions seen: {pasa_versions}")

    summarize(rows, "trinity_prep_seconds")
    summarize(rows, "pasa_align_seconds")
    summarize(rows, "pasa_total_seconds")
    summarize(rows, "total_runtime_seconds")


if __name__ == "__main__":
    main()
