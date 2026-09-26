#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]   # abu_dhabi nodes lack BMI2: bam2hints/Trinity SIGILL (exit 132)
#SBATCH -p epyc -c 16 --mem 64G -t 12:00:00 -J armPrep --array=0-2
#SBATCH -o logs/%x_%A_%a.out -e logs/%x_%A_%a.err
# Per-genome inputs shared by every predict arm (see REPORT_training_selection.md D3):
#   chromosome split (peer predict_scorer.py split), train-only genome,
#   GeneMark-ES GTF (full genome, braker3 image = pipeline GENEMARK_RUN),
#   Trinity-GG minimap2 splice BAMs (full + train-only), transcript GFF3.
set -euo pipefail
module load apptainer samtools
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
X=$R/do_pasa_rust_vs_perl
A=$X/predict_arms
PIXI=/rhome/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/bin
FSIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.1.sif
GSIF=/bigdata/stajichlab/shared/singularity_cache/braker3-v3.1.1.sif
CASES=(
  "Neurospora_crassa_OR74A GCF_000182925.2_NC12 GCF_000182925.2 Neurospora_crassa"
  "Aspergillus_nidulans_FGSC_A4 GCF_000011425.1_ASM1142v1 GCF_000011425.1 Aspergillus_nidulans"
  "Botrytis_cinerea_B05.10 GCF_000143535.2_ASM14353v4 GCF_000143535.2 Botrytis_cinerea"
  "Cryptococcus_neoformans_H99 GCF_000149245.1_CNA3 GCF_000149245.1 Cryptococcus_neoformans"
  "Schizophyllum_commune_H4-8 GCF_000143185.2_Schco3 GCF_000143185.2 Schizophyllum_commune"
)
read -r NAME ASM ACC SP <<< "${CASES[$SLURM_ARRAY_TASK_ID]}"
G=$A/$NAME; mkdir -p $G; cd $G
export TMPDIR=${SCRATCH:?}

[ -s genome.fa ] || pigz -dc $R/input_clean_genomes/$ASM.masked.fasta.gz > genome.fa
samtools faidx genome.fa
[ -s split.train_chroms.txt ] || $PIXI/python $X/refseq_benchmark/predict_scorer.py split \
    --ref-gff3 $X/refseq/$ACC.gff --out-prefix $G/split
samtools faidx genome.fa $(cat split.train_chroms.txt) > genome_train.fa
samtools faidx genome_train.fa
echo "train: $(tr '\n' ' ' < split.train_chroms.txt)"; echo "holdout: $(tr '\n' ' ' < split.holdout_chroms.txt)"

# Trinity-GG alignments (same aligner/settings for every arm)
TR=$(readlink -f $R/rnaseq_data/$SP.trinity-GG.fasta)
B="--bind $G:$G,$(dirname $TR):$(dirname $TR),$TMPDIR:$TMPDIR"
for ref in genome genome_train; do
  if [ ! -s trinity.$ref.bam ]; then
    apptainer exec $B $FSIF bash -c "minimap2 -ax splice -u b -G 3000 --secondary=no -t 16 $G/$ref.fa $TR | samtools sort -@4 -o $G/trinity.$ref.bam -"
    samtools index trinity.$ref.bam
  fi
done
# transcript alignments GFF3: the rc.1 train output (trinity.alignments.gff3), full + train subset
T=$X/results/$NAME.rust/trinity.alignments.gff3
cp $T transcripts.genome.gff3
awk -F'\t' 'NR==FNR{k[$1]=1;next} /^#/||($1 in k)' split.train_chroms.txt $T > transcripts.genome_train.gff3

# GeneMark-ES, full genome, fungal mode (identical input for every arm)
if [ ! -s genemark.genome.gtf ]; then
  GMK=$(readlink -f $HOME/.gm_key); mkdir -p gmes; cd gmes
  apptainer exec --bind $G:$G,$(dirname $GMK):$(dirname $GMK),$TMPDIR:$TMPDIR $GSIF \
      gmes_petap.pl --ES --fungus --cores 16 --sequence $G/genome.fa > gmes.log 2>&1
  cp genemark.gtf ../genemark.genome.gtf; cd ..
fi
awk -F'\t' 'NR==FNR{k[$1]=1;next} ($1 in k)' split.train_chroms.txt genemark.genome.gtf > genemark.genome_train.gtf
echo "prep done $NAME: genemark genes=$(awk -F'\t' '$3=="CDS"' genemark.genome.gtf | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)"
