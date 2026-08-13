#!/usr/bin/env python3
"""Compare baseline vs. GeneMark-disabled funannotate predict reruns.

For each genome in outputs/predict_runs/<OUT>__{baseline,nogenemark}/, parses
the raw EVM-input model counts and final gene counts from
logfiles/funannotate-predict.log, then computes genomic-coordinate overlap
between the two variants' final gene sets (bedtools intersect) to find genes
present in one variant but not the other -- the actual GeneMark-attributable
gene set, not just a raw count delta (EVM consensus can substitute a
GeneMark model for what would otherwise be an Augustus/snap model at the
same locus, which a count-only comparison would miss).

This is a deterministic paired rerun (same genome, same training data, same
funannotate binary, only -w genemark:0 differs) -- not a sampled measurement,
so no bootstrap/permutation testing applies; the reportable quantity is the
exact set difference between two runs, not an estimate with a confidence
interval.
"""
import argparse
import ast
import csv
import re
import subprocess
import sys
from pathlib import Path

GENE_SUMMARY_RE = re.compile(r"Summary of gene models: (\{.*\})")
FINAL_COUNT_RE = re.compile(r"Collecting final annotation files for ([\d,]+) total gene models")
EVM_COUNT_RE = re.compile(r"([\d,]+) total gene models from EVM")


def parse_log(log_path):
    """Return (raw_model_summary_dict_or_None, evm_count_or_None, final_count_or_None)."""
    if not log_path.exists():
        return None, None, None
    text = log_path.read_text(errors="replace")
    summary = None
    for m in GENE_SUMMARY_RE.finditer(text):
        # last match wins (some runs log a dry-run summary earlier)
        summary = ast.literal_eval(m.group(1))
    evm_count = None
    for m in EVM_COUNT_RE.finditer(text):
        evm_count = int(m.group(1).replace(",", ""))
    final_count = None
    for m in FINAL_COUNT_RE.finditer(text):
        final_count = int(m.group(1).replace(",", ""))
    return summary, evm_count, final_count


def gene_bed(gff3_path, out_bed):
    """Extract gene-feature intervals (chrom, start0, end, id) to a sorted BED."""
    rows = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            gid_m = re.search(r"ID=([^;]+)", f[8])
            gid = gid_m.group(1) if gid_m else "NA"
            rows.append((chrom, start, end, gid, 0, strand))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_bed, "w") as out:
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    return len(rows)


def unique_to_a(bed_a, bed_b):
    """Count features in bed_a with zero overlap in bed_b (bedtools intersect -v)."""
    p = subprocess.run(
        ["bedtools", "intersect", "-sorted", "-v", "-a", str(bed_a), "-b", str(bed_b)],
        capture_output=True, text=True, check=True,
    )
    return len([l for l in p.stdout.splitlines() if l.strip()])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs-dir", default=Path(__file__).parent.parent / "outputs" / "predict_runs")
    ap.add_argument("--out", default=Path(__file__).parent.parent / "outputs" / "genemark_contribution_summary.csv")
    args = ap.parse_args()
    runs_dir = Path(args.runs_dir)

    # Discover every genome with a completed baseline rerun under this analysis's own
    # outputs/predict_runs/ -- run.sh writes exactly one baseline dir per genome in
    # GENOMES, so this must pick up all of them rather than a hardcoded fixed list.
    outs = sorted({p.name.rsplit("__", 1)[0] for p in runs_dir.glob("*__baseline")})  # ANALYSIS_OK[file-selection]: enumerates all completed runs in this analysis's own dir, not a latest/backup pick; verified against GENOMES list in run.sh
    if not outs:
        print(f"No baseline runs found under {runs_dir}", file=sys.stderr)
        sys.exit(1)

    rows = []
    for out_name in outs:
        base_dir = runs_dir / f"{out_name}__baseline"
        nogm_dir = runs_dir / f"{out_name}__nogenemark"
        base_log = base_dir / "logfiles" / "funannotate-predict.log"
        nogm_log = nogm_dir / "logfiles" / "funannotate-predict.log"

        base_summary, base_evm, base_final = parse_log(base_log)
        nogm_summary, nogm_evm, nogm_final = parse_log(nogm_log)

        row = {
            "genome": out_name,
            "baseline_raw_total": base_summary.get("total") if base_summary else None,
            "baseline_genemark_raw": base_summary.get("GeneMark") if base_summary else None,
            "baseline_evm_count": base_evm,
            "baseline_final_gene_count": base_final,
            "nogenemark_raw_total": nogm_summary.get("total") if nogm_summary else None,
            "nogenemark_genemark_raw": nogm_summary.get("GeneMark") if nogm_summary else None,
            "nogenemark_evm_count": nogm_evm,
            "nogenemark_final_gene_count": nogm_final,
        }

        base_gff = base_dir / "predict_results" / f"{out_name}.gff3"
        nogm_gff = nogm_dir / "predict_results" / f"{out_name}.gff3"
        if base_gff.exists() and nogm_gff.exists():
            base_bed = base_dir / "genes.bed"
            nogm_bed = nogm_dir / "genes.bed"
            n_base = gene_bed(base_gff, base_bed)
            n_nogm = gene_bed(nogm_gff, nogm_bed)
            only_base = unique_to_a(base_bed, nogm_bed)
            only_nogm = unique_to_a(nogm_bed, base_bed)
            row.update({
                "baseline_gff_gene_count": n_base,
                "nogenemark_gff_gene_count": n_nogm,
                "genes_only_in_baseline": only_base,
                "genes_only_in_nogenemark": only_nogm,
                "genes_shared_by_coordinate": n_base - only_base,
            })
        else:
            missing = [str(p) for p in (base_gff, nogm_gff) if not p.exists()]
            print(f"[WARN] {out_name}: missing gff3 ({missing}); skipping coordinate overlap", file=sys.stderr)

        rows.append(row)

    fieldnames = list(rows[0].keys())
    for r in rows[1:]:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {args.out}")
    for r in rows:
        print(f"\n{r['genome']}:")
        for k, v in r.items():
            if k != "genome":
                print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
