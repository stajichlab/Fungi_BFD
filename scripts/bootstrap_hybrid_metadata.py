#!/usr/bin/env python3
"""bootstrap_hybrid_metadata.py — generator/linter for nextflow/assets/hybrid_parentage.csv,
the file consumed by the hybrid-cross RNA-seq handling in FUNANNOTATE_RNASEQ.nf. See
nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md for the full design.

Writes hybrid_parentage.csv: one row per (hybrid_species_tag, parent_species) pair,
narrow/long format so 2-way through N-way crosses need no schema change. Populated
by splitting SPECIES on the "x"/"×" token for the composite-cross form ("Genus
speciesA x Genus speciesB[...]"); the formal nothospecies shorthand ("Genus x
epithet", e.g. "Saccharomyces x bayanus") doesn't decompose from the string alone
and is resolved via NOTHOSPECIES_PARENTS below instead -- literature-sourced
parentage, not inferred.

samples.csv itself is read-only here and never modified -- hybrid status is
detected by scanning its SPECIES column each run (regex below), not by a persisted
column, since interspecific hybrids are <0.5% of rows and hybrid_parentage.csv's
membership already *is* the runtime signal FUNANNOTATE_RNASEQ.nf acts on (see
loadHybridParentage() in nextflow/modules/funannotate/utils.nf).

Re-run whenever samples.csv gains a new hybrid row; existing hybrid_parentage.csv
rows for species already present are left untouched (see --no-clobber below) so a
manually-corrected entry never gets silently overwritten by a re-run.

Usage:
    scripts/bootstrap_hybrid_metadata.py --samples samples.csv \\
        --parentage-csv nextflow/assets/hybrid_parentage.csv

    # Lint only -- report drift, write nothing, nonzero exit if any species is
    # unresolved (no hybrid_parentage.csv rows and not decomposable/looked up):
    scripts/bootstrap_hybrid_metadata.py --samples samples.csv \\
        --parentage-csv nextflow/assets/hybrid_parentage.csv --lint
"""
import argparse
import csv
import re
import sys

# Literature-sourced parentage for nothospecies shorthand names that don't
# decompose from the SPECIES string itself. One entry today -- extend only with
# a citation-backed parentage, never a guess.
#
# Saccharomyces x bayanus (historically S. pastorianus / S. carlsbergensis): the
# lager-yeast hybrid lineage. Three parents, not two -- the CBS 380 type strain
# is documented as an S. uvarum x S. eubayanus hybrid with a minor S. cerevisiae
# contribution, per bioinformatics review 2026-08-28 (original 2-parent entry
# under-represented its actual composition).
NOTHOSPECIES_PARENTS = {
    "Saccharomyces x bayanus": [
        "Saccharomyces cerevisiae",
        "Saccharomyces eubayanus",
        "Saccharomyces uvarum",
    ],
}

HYBRID_TOKEN_RE = re.compile(r"(?:\sx\s|×)", re.IGNORECASE)


def species_tag(species: str) -> str:
    return re.sub(r"\s+", "_", species.strip())


def is_hybrid(species: str) -> bool:
    return bool(HYBRID_TOKEN_RE.search(species))


def resolve_parents(species: str) -> list:
    """Return the list of parent SPECIES binomials for a hybrid SPECIES string,
    or [] if it can't be resolved (falls through to the genus-wide runtime
    fallback in FUNANNOTATE_RNASEQ.nf instead)."""
    if species in NOTHOSPECIES_PARENTS:
        return NOTHOSPECIES_PARENTS[species]
    parts = [p.strip() for p in re.split(r"\sx\s", species) if p.strip()]
    if len(parts) < 2:
        return []
    if all(len(p.split()) >= 2 for p in parts):
        return parts
    return []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="samples.csv to scan (read-only)")
    ap.add_argument("--parentage-csv", required=True, help="hybrid_parentage.csv to write/update")
    ap.add_argument("--no-clobber", action="store_true", default=True,
                     help="keep existing hybrid_parentage.csv rows for species already present (default on)")
    ap.add_argument("--lint", action="store_true",
                     help="report drift only, write nothing; nonzero exit if any hybrid species is unresolved")
    args = ap.parse_args()

    with open(args.samples, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    if "SPECIES" not in fieldnames:
        sys.exit(f"ERROR: {args.samples} has no SPECIES column")

    existing_parentage = {}
    try:
        with open(args.parentage_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                existing_parentage.setdefault(row["hybrid_species_tag"], []).append(row["parent_species"])
    except FileNotFoundError:
        pass

    n_hybrid = 0
    n_resolved = 0
    n_unresolved = []
    new_parentage_rows = []
    seen_tags = set()

    for row in rows:
        species = (row.get("SPECIES") or "").strip()
        if not is_hybrid(species):
            continue
        n_hybrid += 1
        tag = species_tag(species)
        if tag in seen_tags:
            continue
        seen_tags.add(tag)
        if args.no_clobber and tag in existing_parentage:
            continue
        parents = resolve_parents(species)
        if parents:
            n_resolved += 1
            for p in parents:
                new_parentage_rows.append((tag, p))
        else:
            n_unresolved.append(species)

    print(f"[INFO] {len(rows)} samples.csv rows scanned; {n_hybrid} hybrid rows across {len(seen_tags)} distinct species")
    print(f"[INFO] {n_resolved} hybrid species resolved to parents; {len(n_unresolved)} unresolved (will use genus-wide runtime fallback):")
    for s in n_unresolved:
        print(f"       UNRESOLVED: {s}")

    if args.lint:
        if n_unresolved:
            sys.exit(f"[LINT FAIL] {len(n_unresolved)} hybrid species have no hybrid_parentage.csv coverage and don't auto-resolve")
        print("[LINT OK] every hybrid species in samples.csv resolves to hybrid_parentage.csv rows (existing or newly resolvable)")
        return

    all_parentage_rows = [
        (tag, p) for tag, parents in existing_parentage.items() for p in parents
    ] + new_parentage_rows
    all_parentage_rows.sort()
    # lineterminator='\n' -- csv.writer defaults to '\r\n' even on Unix, which
    # would silently CRLF-ify this file against every other LF-only file in
    # the repo (confirmed 2026-08-28: an earlier run of this script did exactly
    # that to both this file and samples.csv, making every unrelated line show
    # as changed in `git diff`).
    with open(args.parentage_csv, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["hybrid_species_tag", "parent_species"])
        writer.writerows(all_parentage_rows)
    print(f"[INFO] Wrote {len(all_parentage_rows)} rows to {args.parentage_csv}")


if __name__ == "__main__":
    main()
