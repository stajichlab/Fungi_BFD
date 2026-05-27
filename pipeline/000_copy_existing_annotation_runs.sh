#!/usr/bin/bash -l
# Copy pre-existing funannotate annotation folders for samples that have no RNAseq
# training data (zero-byte norm files) from the 1KFG common annotation store.
#
# Skips:
#   - samples where TRANSL_TABLE=12 (non-standard codon table; need re-annotation)
#   - destinations that already exist in genome_annotation/
#   - basenames with no matching folder in the source annotation store
#
# Run from the project root (no SLURM needed — this is just rsync):
#   bash pipeline/000_copy_existing_annotation_runs.sh
#   bash pipeline/000_copy_existing_annotation_runs.sh --dry-run

set -euo pipefail

DRY_RUN=0
for arg in "$@"; do
    [[ "$arg" == "--dry-run" ]] && DRY_RUN=1
done
[[ "$DRY_RUN" == "1" ]] && echo "[DRY-RUN] No files will be copied."

RNASEQ_DIR="rnaseq_reads"
SAMPLES_CSV="samples.csv"
SRC_BASE="/bigdata/stajichlab/shared/projects/1KFG/common_annotate/annotate/"
DST_BASE="genome_annotation"

mkdir -p "$DST_BASE"

# Build lookups from samples.csv keyed on SPECIES_underscored:
#   TRANSL12: species -> 1 if any row has TRANSL_TABLE=12
#   STRAINS:  species -> pipe-delimited list of STRAIN values (for folder naming)
declare -A TRANSL12
declare -A STRAINS

while IFS=',' read -r ASMID SPECIES_IN STRAIN BIOPROJECT TAXONID BUSCO PHYLUM SUBPHYLUM \
                       CLASS SUBCLASS ORDER FAMILY GENUS SPECIES TRANSL_TABLE LOCUSTAG; do
    [[ "$ASMID" == "ASMID" ]] && continue   # skip header
    key="${SPECIES// /_}"
    if [[ "$TRANSL_TABLE" == "12" ]]; then
        TRANSL12["$key"]=1
    else
        # Only set to 0 if not already marked as 12
        : "${TRANSL12[$key]:=0}"
    fi
    # Accumulate strains for this species (pipe-delimited)
    strain_norm="${STRAIN// /_}"
    STRAINS["$key"]+="${strain_norm}|"
done < "$SAMPLES_CSV"

copied=0
skipped_transl=0
skipped_exists=0
skipped_nosrc=0

for r1 in "$RNASEQ_DIR"/*_norm_R1.fastq.gz; do
    [[ -f "$r1" ]] || continue

    # Only process zero-byte norm files (no RNAseq data available)
    [[ -s "$r1" ]] && continue

    # RNAseq filenames use species only (no strain)
    species="${r1##*/}"
    species="${species%_norm_R1.fastq.gz}"

    # Check TRANSL_TABLE for this species
    if [[ "${TRANSL12[$species]+set}" == "set" && "${TRANSL12[$species]}" == "1" ]]; then
        echo "[SKIP-TRANSL12]  $species"
        (( skipped_transl++ )) || true
        continue
    fi

    # Annotation folders include strain: Genus_species_STRAIN
    # A species may have multiple strains; iterate over all of them.
    if [[ -z "${STRAINS[$species]+set}" ]]; then
        echo "[NO-SPECIES]     $species (not found in $SAMPLES_CSV)"
        (( skipped_nosrc++ )) || true
        continue
    fi

    IFS='|' read -ra strain_list <<< "${STRAINS[$species]%|}"
    for strain in "${strain_list[@]}"; do
        [[ -z "$strain" ]] && continue

        if [[ "$strain" == "NA" || "$strain" == "na" ]]; then
            folder="$species"
        else
            folder="${species}_${strain}"
        fi

        src="$SRC_BASE/$folder"
        dst="$DST_BASE/$folder"

        # Skip if source annotation folder does not exist
        if [[ ! -d "$src" ]]; then
            echo "[NO-SOURCE]      $folder  (looked in $src)"
            (( skipped_nosrc++ )) || true
            continue
        fi

        # Skip if destination already exists
        if [[ -d "$dst" ]]; then
            echo "[ALREADY-EXISTS] $folder"
            (( skipped_exists++ )) || true
            continue
        fi

        if [[ "$DRY_RUN" == "1" ]]; then
            echo "[DRY-RUN]        rsync -av --progress --exclude annotate_misc $src/ $dst/"
        else
            echo "[COPYING]        $folder"
            rsync -a --exclude annotate_misc "$src/" "$dst/"
        fi
        (( copied++ )) || true
    done
done

echo ""
echo "Done."
if [[ "$DRY_RUN" == "1" ]]; then
    echo "  Would copy:          $copied"
else
    echo "  Copied:              $copied"
fi
echo "  Skipped (TRANSL12):  $skipped_transl"
echo "  Skipped (exists):    $skipped_exists"
echo "  Skipped (no source): $skipped_nosrc"
