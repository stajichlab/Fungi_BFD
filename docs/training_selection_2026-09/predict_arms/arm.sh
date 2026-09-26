#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]   # abu_dhabi nodes lack BMI2: bam2hints/Trinity SIGILL (exit 132)
#SBATCH -p epyc -c 16 --mem 64G -t 16:00:00
#SBATCH -o logs/%x_%j.out -e logs/%x_%j.err
# One predict arm. Usage: sbatch -J NAME arm.sh GENOME ARM STEP EVID
#   GENOME  Neurospora_crassa_OR74A | Aspergillus_nidulans_FGSC_A4 | Botrytis_cinerea_B05.10
#   ARM     old | oldnew | tx | cdsS | cdsB | busco
#   STEP    A (train on train chromosomes) | B (full genome, -p params from A)
#   EVID    step B PASA evidence: fixed (old getBestModel for every arm) | own (arm's own)
# See REPORT_training_selection.md D3 and DECISIONS.md D11.
set -euo pipefail
module load apptainer
GENOME=$1; ARM=$2; STEP=$3; EVID=${4:-fixed}
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
X=$R/do_pasa_rust_vs_perl; A=$X/predict_arms; BM=$X/refseq_benchmark
G=$A/$GENOME
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.1.sif
SITE=/pixi/.pixi/envs/base/lib/python3.8/site-packages/funannotate
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
case $GENOME in
  Neurospora_crassa_OR74A)       SPECIES="Neurospora crassa";    STRAIN="OR74A";   TAG=FC69C3D3; ACC=GCF_000182925.2 ;;
  Aspergillus_nidulans_FGSC_A4)  SPECIES="Aspergillus nidulans"; STRAIN="FGSC A4"; TAG=F197696E; ACC=GCF_000011425.1 ;;
  Botrytis_cinerea_B05.10)       SPECIES="Botrytis cinerea";     STRAIN="B05.10";  TAG=F9A9046A; ACC=GCF_000143535.2 ;;
  Cryptococcus_neoformans_H99)   SPECIES="Cryptococcus neoformans"; STRAIN="H99";  TAG=F416552B; ACC=GCF_000149245.1 ;;
  Schizophyllum_commune_H4-8)    SPECIES="Schizophyllum commune"; STRAIN="H4-8";   TAG=F43E66A0; ACC=GCF_000143185.2 ;;
esac
case $ARM in
  old) VAR=old ;; oldnew) VAR=old ;; tx) VAR=tx_strand ;; cdsS) VAR=cds_strand ;; cdsB) VAR=cds_blind ;; busco) VAR=tx_strand ;;
  tx2) VAR=tx_strand ;; se) VAR=tx_strand ;;   # code_new2 snapshot (review fixes + R6 option b), D45
  txG) VAR=tx_strand; PSRC=G ;;   # code_new3 (guarded getBestModel ranking, D55)
  txRel) VAR=tx_strand; PSRC=REL ;;   # code_new + relaxed PASA validation (xs/mm2fix_blat_relaxed)
  seR1R2) VAR=tx_strand; PSRC=R1R2 ;;     # code_new2 + --training_single_exon + R1R2 PASA (D67)
  buscoR1R2) VAR=tx_strand; PSRC=R1R2 ;;  # code_new, BUSCO-forced training, R1R2 PASA evidence (D67)
  busco4) VAR=tx_strand ;;   # code_new4 (titration code), BUSCO-forced; REP env var names repeat dirs (D96)
  txR1) VAR=tx_strand; PSRC=R1 ;; txR1R2) VAR=tx_strand; PSRC=R1R2 ;; txID90) VAR=tx_strand; PSRC=ID90 ;;   # code_new + fixed-PASA inputs (pasa_inputs/<src>/)
esac
AD=$ARM${REP:+_r$REP}   # repeat runs: REP=1,2,3 -> busco4_r1.A etc.
W=$G/$AD.$STEP; [ $STEP = B ] && W=$G/$AD.B.$EVID
mkdir -p $W; cd $W
export TMPDIR=${SCRATCH:?}
rm -rf $W/augustus_config $W/augustus $W/out
# basename must be 'config': predict.py only sets AUGUSTUS_BASE when it is
mkdir -p $W/augustus; cp -r $R/lib/augustus/3.5/config $W/augustus/config

EXTRA=()
BINDS="$W:$W,$G:$G,$BM:$BM,$A:$A,$R/lib:$R/lib,$FDB:$FDB,$TMPDIR:$TMPDIR"
if [ $STEP = A ]; then
  GEN=$G/genome_train.fa; BAM=$G/trinity.genome_train.bam; TRX=$G/transcripts.genome_train.gff3; GMK=$G/genemark.genome_train.gtf
  SRCBM=$BM; [ -n "${PSRC:-}" ] && SRCBM=$A/pasa_inputs/$PSRC
  awk -F'\t' 'NR==FNR{k[$1]=1;next} /^#/||($1 in k)' $G/split.train_chroms.txt $SRCBM/$GENOME.best.$VAR.gff3 > $W/pasa.gff3
  if [ $ARM != old ]; then
    # new selection code (frozen snapshot) over the image's funannotate package
    CODE=$A/code_new; [[ $ARM == tx2 || $ARM == se || $ARM == seR1R2 ]] && CODE=$A/code_new2; [[ $ARM == txG ]] && CODE=$A/code_new3; [[ $ARM == busco4 ]] && CODE=$A/code_new4
    BINDS="$BINDS,$CODE/funannotate:$SITE"
    if [[ $ARM == busco || $ARM == buscoR1R2 || $ARM == busco4 ]]; then EXTRA+=(--min_pasa_complete_models 1000000000); else EXTRA+=(--min_pasa_complete_models 0); fi
    [[ $ARM == se || $ARM == seR1R2 ]] && EXTRA+=(--training_single_exon)
  fi
else
  GEN=$G/genome.fa; BAM=$G/trinity.genome.bam; TRX=$G/transcripts.genome.gff3; GMK=$G/genemark.genome.gtf
  EV=old; [ $EVID = own ] && EV=$VAR
  SRCBM=$BM; [ -n "${PSRC:-}" ] && [ $EVID = own ] && SRCBM=$A/pasa_inputs/$PSRC
  cp $SRCBM/$GENOME.best.$EV.gff3 $W/pasa.gff3
  P=$(ls $G/$AD.A/out/predict_results/*.parameters.json)
  EXTRA+=(-p "$P")
fi
echo "arm=$ARM step=$STEP evid=$EVID pasa_variant=$VAR genome=$GEN extra=${EXTRA[*]:-none}" | tee $W/arm_info.txt

export APPTAINERENV_FUNANNOTATE_DB=$FDB
export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/augustus/config
start=$(date +%s)
set +e
apptainer exec --bind $BINDS $SIF funannotate predict --name $TAG -i $GEN --strain "$STRAIN" \
    -o $W/out -s "$SPECIES" --cpu 16 --busco_db dikarya \
    --AUGUSTUS_CONFIG_PATH $W/augustus/config -w codingquarry:0 glimmerhmm:0 genemark:1 \
    --min_training_models 30 --tmpdir $TMPDIR --SeqCenter NCBI --keep_no_stops --header_length 24 \
    --protein_evidence $R/lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10 \
    --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf $GMK \
    --rna_bam $BAM --pasa_gff $W/pasa.gff3 --transcript_alignments $TRX "${EXTRA[@]}" \
    > $W/predict.capture.log 2>&1
status=$?
set -e
echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
[ $status -eq 0 ] || { echo "predict failed ($status)"; exit $status; }

if [ $STEP = B ]; then
  PRED=$(ls $W/out/predict_results/*.gff3 | grep -v '\.tbl' | head -1)
  module load gffcompare
  /rhome/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/bin/python $BM/predict_scorer.py score \
      --ref-gff3 $X/refseq/$ACC.gff --pred-gff3 $PRED --holdout $G/split.holdout_chroms.txt \
      --genome $GENOME --variant "$AD.B.$EVID" --out $W/score.tsv
  module load busco/5.8.0
  busco -i $(ls $W/out/predict_results/*.proteins.fa) -m proteins -l /srv/projects/db/BUSCO/v10/lineages/fungi_odb10 \
      -o busco -c 16 --offline -f --out_path $TMPDIR > $W/busco.log 2>&1 || true
  cp $TMPDIR/busco/short_summary*.txt $W/busco_short_summary.txt 2>/dev/null || true
fi
echo done
