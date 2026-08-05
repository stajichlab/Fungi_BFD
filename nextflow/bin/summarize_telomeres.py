#!/usr/bin/env python3
"""Aggregate per-genome telomere feature tables into merged summary + detail tables.

Inputs are supplied as a manifest file (--manifest, one TAB-delimited
``<path>\\t<mtime>\\t<size>`` record per line) or as a directory to glob
(--reportdir).  Two outputs are written:

  --summary-out  one row per ASMID: total tract length/count, scaffold counts,
                 5'/3'/strand counts, top monomer, and telomere_scaffolds_both_ends
                 (scaffolds with a *terminal* tract on both ends -- the
                 chromosome-completeness signal).
  --tracts-out   one row per raw tract (scaffold/end_type/strand/monomer/sequence),
                 for sequence- and length-level mining and ad hoc both-ends
                 queries. Carries a surrogate tract_id since a scaffold can, in
                 principle, carry more than one tract call at the same end
                 (e.g. degenerate repeats).

Both outputs carry ASMID and LOCUSTAG (looked up from --samples) so they join
directly to the `species` table without going through asm_stats.
"""

import argparse
import csv
import glob
import gzip
import os
import sys
from collections import Counter


SUMMARY_COLUMNS = [
    "ASMID", "LOCUSTAG", "SPECIES", "STRAIN",
    "telomere_scaffolds",
    "telomere_scaffolds_both_ends",
    "telomere_tracts",
    "telomere_total_length_bp",
    "telomere_total_repeats",
    "telomere_5prime_count",
    "telomere_3prime_count",
    "telomere_plus_count",
    "telomere_minus_count",
    "telomere_terminal_count",
    "telomere_internal_count",
    "telomere_top_monomer",
    "telomere_monomers",
]

TRACT_COLUMNS = [
    "tract_id", "ASMID", "LOCUSTAG",
    "scaffold", "end_type", "strand", "monomer",
    "repeat_count", "tract_length",
    "start", "end_coord", "terminal", "distance_to_end",
    "tract_seq", "flank_seq",
]


def parse_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def parse_bool(value):
    return value.strip().lower() in ("true", "1", "yes", "t")


def load_samples(path):
    """Return a dict mapping ASMID -> {LOCUSTAG, SPECIES, STRAIN}."""
    samples = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = row["ASMID"]
            samples[asmid] = {
                "LOCUSTAG": row.get("LOCUSTAG", ""),
                "SPECIES": row.get("SPECIES_IN", row.get("SPECIES", "")),
                "STRAIN": row.get("STRAIN", ""),
            }
    return samples


def read_manifest(path):
    """Read a toManifest file; return the list of report paths."""
    files = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            files.append(line.split("\t", 1)[0])
    return files


def open_in(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def open_out(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def process_genome(path, asmid, locustag):
    """Read one genome's telomere TSV; return (summary_dict, list-of-tract-rows)."""
    scaffolds = set()
    both_ends_scaffolds = {"5prime": set(), "3prime": set()}
    total_length = 0
    total_repeats = 0
    end_5 = 0
    end_3 = 0
    strand_plus = 0
    strand_minus = 0
    terminal = 0
    internal = 0
    monomer_counts = Counter()
    tract_rows = []

    with open_in(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            scaffold = row["scaffold"]
            end = row["end"]
            is_terminal = parse_bool(row.get("terminal", "True"))

            scaffolds.add(scaffold)
            total_length += parse_int(row["tract_length"])
            total_repeats += parse_int(row["repeat_count"])
            if end == "5prime":
                end_5 += 1
            elif end == "3prime":
                end_3 += 1
            if row["strand"] == "+":
                strand_plus += 1
            elif row["strand"] == "-":
                strand_minus += 1
            if is_terminal:
                terminal += 1
                # Only a *terminal* call counts toward "both ends telomere-capped" --
                # an interstitial telomeric sequence (ITS) at an internal locus is
                # not evidence the scaffold end itself is capped.
                if end in both_ends_scaffolds:
                    both_ends_scaffolds[end].add(scaffold)
            else:
                internal += 1
            monomer_counts[row["monomer"]] += 1

            tract_rows.append({
                "tract_id": f"{asmid}:{scaffold}:{end}:{row['start']}",
                "ASMID": asmid,
                "LOCUSTAG": locustag,
                "scaffold": scaffold,
                "end_type": end,
                "strand": row["strand"],
                "monomer": row["monomer"],
                "repeat_count": parse_int(row["repeat_count"]),
                "tract_length": parse_int(row["tract_length"]),
                "start": parse_int(row["start"]),
                "end_coord": parse_int(row["end_coord"]),
                "terminal": is_terminal,
                "distance_to_end": parse_int(row["distance_to_end"]),
                "tract_seq": row.get("tract_seq", ""),
                "flank_seq": row.get("flank_seq", ""),
            })

    both_ends = both_ends_scaffolds["5prime"] & both_ends_scaffolds["3prime"]
    top_monomer = monomer_counts.most_common(1)[0][0] if monomer_counts else ""
    summary = {
        "telomere_scaffolds": len(scaffolds),
        "telomere_scaffolds_both_ends": len(both_ends),
        "telomere_tracts": sum(monomer_counts.values()),
        "telomere_total_length_bp": total_length,
        "telomere_total_repeats": total_repeats,
        "telomere_5prime_count": end_5,
        "telomere_3prime_count": end_3,
        "telomere_plus_count": strand_plus,
        "telomere_minus_count": strand_minus,
        "telomere_terminal_count": terminal,
        "telomere_internal_count": internal,
        "telomere_top_monomer": top_monomer,
        "telomere_monomers": ";".join(sorted(monomer_counts)),
    }
    return summary, tract_rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--manifest", help="TAB-delimited manifest of report paths")
    src.add_argument("--reportdir", help="Directory to glob for *.telomeres.tsv.gz")
    p.add_argument("--samples", required=True, help="samples CSV with ASMID/LOCUSTAG/SPECIES_IN/STRAIN")
    p.add_argument("--summary-out", required=True, help="output TSV (.gz -> gzip-compressed), one row per genome")
    p.add_argument("--tracts-out", required=True, help="output TSV (.gz -> gzip-compressed), one row per tract")
    args = p.parse_args()

    samples = load_samples(args.samples)

    if args.manifest:
        report_files = read_manifest(args.manifest)
    else:
        report_files = glob.glob(os.path.join(args.reportdir, "**/*.telomeres.tsv.gz"), recursive=True)
    report_files = sorted(report_files)

    if not report_files:
        src_desc = args.manifest or args.reportdir
        print(f"No telomere report files found ({src_desc})", file=sys.stderr)
        sys.exit(1)

    summary_rows = []
    tract_rows = []
    missing_meta = []
    for path in report_files:
        stem = os.path.basename(path).replace(".telomeres.tsv.gz", "").replace(".telomeres.tsv", "")

        meta = samples.get(stem)
        if meta is None:
            missing_meta.append(stem)
            locustag, species, strain = "", "", ""
        else:
            locustag = meta["LOCUSTAG"]
            species = meta["SPECIES"]
            strain = meta["STRAIN"]

        summary, tracts = process_genome(path, stem, locustag)

        row = {"ASMID": stem, "LOCUSTAG": locustag, "SPECIES": species, "STRAIN": strain}
        row.update(summary)
        summary_rows.append(row)
        tract_rows.extend(tracts)

    if missing_meta:
        print(f"WARNING: {len(missing_meta)} ASMIDs not found in {args.samples}:", file=sys.stderr)
        for m in missing_meta[:10]:
            print(f"  {m}", file=sys.stderr)
        if len(missing_meta) > 10:
            print(f"  ... and {len(missing_meta) - 10} more", file=sys.stderr)

    with open_out(args.summary_out) as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary_rows)

    with open_out(args.tracts_out) as fh:
        writer = csv.DictWriter(fh, fieldnames=TRACT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(tract_rows)

    print(f"Wrote {len(summary_rows)} genome summary rows to {args.summary_out}")
    print(f"Wrote {len(tract_rows)} tract rows to {args.tracts_out}")


if __name__ == "__main__":
    main()
