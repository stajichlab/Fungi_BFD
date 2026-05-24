#!/usr/bin/bash -l
#SBATCH -p epyc -c 128 --mem 96gb --out logs/fastmask.log

module load funannotate

INFOLDER="${1:-input_clean_genomes}"
JOBS=64 # parallel jobs; tune relative to --cpus 2 and core count

mask_one() {
    genome_fa="$1"
    STEM="${genome_fa%.fa}"
    STEM="${STEM%.fasta}"
    outfile="${STEM}.masked.fasta"
    if [[ -f "$outfile" ]]; then
        echo "SKIP: $outfile already exists"
        return 0
    fi
    echo "Masking: $genome_fa"
    funannotate mask -i "$genome_fa" -o "$outfile" -m tantan --cpus 2
}

export -f mask_one

find "$INFOLDER" -maxdepth 1 -name "*.fa" -o -name "*.fasta" \
  | parallel -j "$JOBS" mask_one {}
