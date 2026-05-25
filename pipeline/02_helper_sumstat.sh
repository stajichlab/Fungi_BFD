#!/usr/bin/bash -l
#SBATCH -c 128 --mem 16gb --out logs/sumstat.log --time 2:00:00 -p short


module load AAFTF

INFOLDER="${1:-input_clean_genomes}"
JOBS=128 # parallel jobs; tune relative to --cpus 2 and core count

sumstats() {
    genome_fa="$1"
    outdir="$2"
    STEM=$(basename ${genome_fa} .masked.fasta)
    outfile="$outdir/${STEM}.stats.txt"

#    echo "${genome_fa} - $STEM - $outfile"
    if [[ -f "$outfile" ]]; then
        echo "SKIP: $outfile already exists"
        return 0
    fi
    echo "Summary stats: $genome_fa"
    AAFTF assess -i "$genome_fa" -r "$outfile" 
}

export -f sumstats
mkdir -p results/asm_reports
find "$INFOLDER" -maxdepth 1 -name "*.masked.fasta" \
  | parallel -j "$JOBS" sumstats {} results/asm_reports
