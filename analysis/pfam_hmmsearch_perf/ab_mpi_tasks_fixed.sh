#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 4 -c 1 --mem 16gb --time 4:00:00 --job-name pfam_ab_mpi_fixed
#SBATCH --out logs/pfam_ab_mpi_fixed.%j.log

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load hmmer/3.4-mpi
module load db-pfam

OUTDIR="analysis/pfam_hmmsearch_perf"
PROT="input/pep/Malassezia_brasiliensis_CBS_14135.proteins.fa"

echo "=== TEST C (fixed): hmmsearch --mpi, 4 MPI tasks, --cpu-bind=none ==="
/usr/bin/time -v srun --cpu-bind=none -N 1 -n 4 hmmsearch --mpi --cut_ga --noali \
    --domtbl "$OUTDIR/testC_mpi4_fixed.domtblout" \
    --tblout "$OUTDIR/testC_mpi4_fixed.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testC_mpi4_fixed.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testC_mpi4_fixed.timing.txt"

echo ""
echo "=== TEST D (fixed): hmmsearch --mpi, 2 MPI tasks, --cpu-bind=none ==="
/usr/bin/time -v srun --cpu-bind=none -N 1 -n 2 hmmsearch --mpi --cut_ga --noali \
    --domtbl "$OUTDIR/testD_mpi2_fixed.domtblout" \
    --tblout "$OUTDIR/testD_mpi2_fixed.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testD_mpi2_fixed.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testD_mpi2_fixed.timing.txt"

echo ""
echo "=== Sanity: domain counts (should roughly match Test A baseline) ==="
grep -vc '^#' "$OUTDIR/testC_mpi4_fixed.tblout" || true
grep -vc '^#' "$OUTDIR/testD_mpi2_fixed.tblout" || true
