#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 16 -c 1 --mem 16gb --time 4:00:00 --job-name pfam_ab_mpi_scale
#SBATCH --out logs/pfam_ab_mpi_scale.%j.log

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load hmmer/3.4-mpi
module load db-pfam

OUTDIR="analysis/pfam_hmmsearch_perf"
PROT="input/pep/Malassezia_brasiliensis_CBS_14135.proteins.fa"

for n in 8 16; do
    echo "=== hmmsearch --mpi, $n MPI tasks, --cpu-bind=none ==="
    /usr/bin/time -v srun --cpu-bind=none -N 1 -n $n hmmsearch --mpi --cut_ga --noali \
        --domtbl "$OUTDIR/testE_mpi${n}.domtblout" \
        --tblout "$OUTDIR/testE_mpi${n}.tblout" \
        "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testE_mpi${n}.timing.txt"
    grep -E "Elapsed|Maximum resident" "$OUTDIR/testE_mpi${n}.timing.txt"
    echo ""
done

echo "=== Sanity: domain counts ==="
for n in 8 16; do
    grep -vc '^#' "$OUTDIR/testE_mpi${n}.tblout" || true
done
