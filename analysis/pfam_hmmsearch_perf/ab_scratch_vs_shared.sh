#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 1 -c 4 --mem 16gb --time 4:00:00 --job-name pfam_ab_scratch
#SBATCH --out logs/pfam_ab_scratch.%j.log

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load hmmer/3.4
module load db-pfam
module load workspace/scratch

OUTDIR="analysis/pfam_hmmsearch_perf"
PROT="input/pep/Malassezia_brasiliensis_CBS_14135.proteins.fa"
CPUS=4

echo "=== Host: $(hostname) ==="
echo "=== PFAM_DB (shared): $PFAM_DB ==="
echo "=== Protein file: $PROT ($(grep -c '^>' "$PROT") proteins) ==="

echo ""
echo "=== Copying Pfam-A DB to scratch ==="
mkdir -p "$SCRATCH/pfam_db"
time cp "$PFAM_DB"/Pfam-A.hmm* "$SCRATCH/pfam_db/"
du -sh "$SCRATCH/pfam_db"

echo ""
echo "=== TEST A: hmmsearch against SHARED storage DB ==="
/usr/bin/time -v hmmsearch --cut_ga --noali --cpu $CPUS \
    --domtbl "$OUTDIR/testA_shared.domtblout" \
    --tblout "$OUTDIR/testA_shared.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testA_shared.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testA_shared.timing.txt"

echo ""
echo "=== TEST B: hmmsearch against SCRATCH-copied DB ==="
/usr/bin/time -v hmmsearch --cut_ga --noali --cpu $CPUS \
    --domtbl "$OUTDIR/testB_scratch.domtblout" \
    --tblout "$OUTDIR/testB_scratch.tblout" \
    "$SCRATCH/pfam_db/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testB_scratch.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testB_scratch.timing.txt"

echo ""
echo "=== TEST A (repeat, warm cache check) ==="
/usr/bin/time -v hmmsearch --cut_ga --noali --cpu $CPUS \
    --domtbl "$OUTDIR/testA2_shared.domtblout" \
    --tblout "$OUTDIR/testA2_shared.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testA2_shared.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testA2_shared.timing.txt"

echo ""
echo "=== Sanity: outputs equivalent (domain counts) ==="
zcat -f "$OUTDIR/testA_shared.tblout" 2>/dev/null | grep -vc '^#' || grep -vc '^#' "$OUTDIR/testA_shared.tblout"
zcat -f "$OUTDIR/testB_scratch.tblout" 2>/dev/null | grep -vc '^#' || grep -vc '^#' "$OUTDIR/testB_scratch.tblout"

echo ""
echo "=== DONE ==="
