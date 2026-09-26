#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]   # abu_dhabi nodes lack BMI2/AVX2 (DECISIONS D49)
#SBATCH -p epyc -c 16 --mem 64G -t 6:00:00 -J titr
#SBATCH -o logs/%x_%A_%a.out -e logs/%x_%A_%a.err
# Experiment A (DECISIONS D81): one titration task = one (genome, N, draw) from tasks.tsv.
#   step A: train Augustus/SNAP on the TRAIN chromosomes from N randomly drawn complete
#           PASA models (code_new4: R3/R5, review fixes, R6 (b) on by default; gate off)
#   step B: full-genome predict with -p params and the SAME fixed evidence as the
#           predict_arms fixed arms (PASA = rc1 old getBestModel), image-native code
#   score : predict_scorer.py on holdout chromosomes + proteome BUSCO -> rows/<task>.tsv
set -uo pipefail
module load apptainer
T=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/predict_arms/titration
A=$(dirname $T); X=$(dirname $A); R=$(dirname $X); BM=$X/refseq_benchmark
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.1.sif
SITE=/pixi/.pixi/envs/base/lib/python3.8/site-packages/funannotate
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
PIXI=/rhome/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/bin
read -r G N D < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $T/tasks.tsv)
case $G in
  Neurospora_crassa_OR74A)       SPECIES="Neurospora crassa";       STRAIN="OR74A";   TAG=FC69C3D3; ACC=GCF_000182925.2 ;;
  Aspergillus_nidulans_FGSC_A4)  SPECIES="Aspergillus nidulans";    STRAIN="FGSC A4"; TAG=F197696E; ACC=GCF_000011425.1 ;;
  Botrytis_cinerea_B05.10)       SPECIES="Botrytis cinerea";        STRAIN="B05.10";  TAG=F9A9046A; ACC=GCF_000143535.2 ;;
  Cryptococcus_neoformans_H99)   SPECIES="Cryptococcus neoformans"; STRAIN="H99";     TAG=F416552B; ACC=GCF_000149245.1 ;;
esac
GD=$A/$G; W=$T/$G/N${N}_d${D}; rm -rf $W; mkdir -p $W; cd $W
export TMPDIR=${SCRATCH:?}
start=$(date +%s)
BINDS="$W:$W,$GD:$GD,$BM:$BM,$X/refseq_benchmark_lowkeep:$X/refseq_benchmark_lowkeep,$T:$T,$R/lib:$R/lib,$FDB:$FDB,$TMPDIR:$TMPDIR"
COMMON=(--name $TAG --strain "$STRAIN" -s "$SPECIES" --cpu 16 --busco_db dikarya
  -w codingquarry:0 glimmerhmm:0 genemark:1 --min_training_models 30 --tmpdir $TMPDIR
  --SeqCenter NCBI --keep_no_stops --header_length 24 --protein_evidence $R/lib/swissprot_fungi.faa
  --max_intronlen 3000 --min_intronlen 10 --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark)

# subset of N complete models (train chromosomes), fixed seed
$PIXI/python $T/subset.py $G $N $D $W/pasa.subset.gff3 > $W/n_subset.txt || { echo "subset failed"; exit 1; }

# step A: training
mkdir -p $W/A/augustus; cp -r $R/lib/augustus/3.5/config $W/A/augustus/config
export APPTAINERENV_FUNANNOTATE_DB=$FDB APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/A/augustus/config
apptainer exec --bind $BINDS,$A/code_new4/funannotate:$SITE $SIF funannotate predict -i $GD/genome_train.fa -o $W/A/out \
  "${COMMON[@]}" --AUGUSTUS_CONFIG_PATH $W/A/augustus/config --genemark_gtf $GD/genemark.genome_train.gtf \
  --rna_bam $GD/trinity.genome_train.bam --pasa_gff $W/pasa.subset.gff3 --transcript_alignments $GD/transcripts.genome_train.gff3 \
  --min_pasa_complete_models 0 > $W/A.capture.log 2>&1
sa=$?
PLOG=$W/A/out/logfiles/funannotate-predict.log
keepers=$(grep -m1 -oE "[0-9,]+ from filterGeneMark" $PLOG 2>/dev/null | tr -d , | cut -d' ' -f1)
ntrain=$(grep -oE "[0-9,]+ of [0-9,]+ models pass training parameters" $PLOG 2>/dev/null | tail -1 | cut -d' ' -f1 | tr -d ,)
single=$(grep -oE "Single-exon training models: [0-9,]+ admitted" $PLOG 2>/dev/null | grep -oE "[0-9,]+" | head -1 | tr -d ,)
P=$(ls $W/A/out/predict_results/*.parameters.json 2>/dev/null | head -1)

# step B: full genome, fixed evidence
sb=1
if [ $sa -eq 0 ] && [ -n "$P" ]; then
  mkdir -p $W/B/augustus; cp -r $R/lib/augustus/3.5/config $W/B/augustus/config
  export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/B/augustus/config
  apptainer exec --bind $BINDS $SIF funannotate predict -i $GD/genome.fa -o $W/B/out \
    "${COMMON[@]}" --AUGUSTUS_CONFIG_PATH $W/B/augustus/config --genemark_gtf $GD/genemark.genome.gtf \
    --rna_bam $GD/trinity.genome.bam --pasa_gff $BM/$G.best.old.gff3 --transcript_alignments $GD/transcripts.genome.gff3 \
    -p "$P" > $W/B.capture.log 2>&1
  sb=$?
fi

# score
row_sc=""
if [ $sb -eq 0 ]; then
  PRED=$(ls $W/B/out/predict_results/*.gff3 | head -1)
  module load gffcompare
  $PIXI/python $X/refseq_benchmark/predict_scorer.py score --ref-gff3 $X/refseq/$ACC.gff --pred-gff3 $PRED \
    --holdout $GD/split.holdout_chroms.txt --genome $G --variant "N${N}_d${D}" --out $W/score.tsv > /dev/null 2>&1
  module load busco/5.8.0
  busco -i $(ls $W/B/out/predict_results/*.proteins.fa) -m proteins -l /srv/projects/db/BUSCO/v10/lineages/fungi_odb10 \
    -o busco -c 16 --offline -f --out_path $TMPDIR > $W/busco.log 2>&1 || true
  cp $TMPDIR/busco/short_summary*.txt $W/busco_short_summary.txt 2>/dev/null || true
fi
$PIXI/python - "$W" "$G" "$N" "$D" "$(cat $W/n_subset.txt 2>/dev/null)" "${ntrain:-}" "${keepers:-}" "${single:-}" "$sa" "$sb" "$(( $(date +%s) - start ))" "$(hostname -s)" <<'EOF'
import csv, os, re, sys
W, g, n, d, nsub, ntrain, keepers, single, sa, sb, wall, node = sys.argv[1:13]
row = {"genome": g, "N": n, "draw": d, "n_subset": nsub, "n_complete_used": ntrain, "n_keepers": keepers,
       "n_single_admitted": single}
sc = f"{W}/score.tsv"
cols = ["locus_sn", "locus_pr", "exon_sn", "exon_pr", "intron_chain_sn", "intron_chain_pr", "pred_genes"]
if os.path.exists(sc):
    rec = list(csv.DictReader(open(sc), delimiter="\t"))[-1]
    row.update({c: rec.get(c, "") for c in cols})
else:
    row.update({c: "" for c in cols})
b = f"{W}/busco_short_summary.txt"
m = re.search(r"C:([\d.]+)%", open(b).read()) if os.path.exists(b) else None
row["busco_c"] = m.group(1) if m else ""
row.update({"exit_stepA": sa, "exit_stepB": sb, "runtime_s": wall, "node": node})
with open(f"{os.path.dirname(os.path.dirname(W))}/rows/{g}_N{n}_d{d}.tsv", "w") as f:
    f.write("\t".join(row) + "\n" + "\t".join(str(v) for v in row.values()) + "\n")
EOF
[ $sa -eq 0 ] && [ $sb -eq 0 ] || exit 1
