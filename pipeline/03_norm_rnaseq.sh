#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 64gb --time 4:00:00
#SBATCH --job-name=norm_rnaseq
#SBATCH --output=logs/norm_rnaseq.%A_%a.log
#SBATCH --array=0-48

# Trim and normalize raw RNA-seq reads already in rnaseq_reads/.
# Input:  rnaseq_reads/<TAG>_R1.fastq.gz  +  <TAG>_R2.fastq.gz  (no "norm" in name)
# Output: rnaseq_reads/<TAG>_norm_R1.fastq.gz + <TAG>_norm_R2.fastq.gz
#
# Submit from the project root:
#   sbatch --array=0-$(( $(ls rnaseq_reads/*_R1.fastq.gz | grep -v norm | wc -l) - 1 ))%20 \
#          pipeline/03_norm_rnaseq.sh
#
# After verifying norm outputs, the following originals can be deleted:
#   rm rnaseq_reads/*_R1.fastq.gz rnaseq_reads/*_R2.fastq.gz
#   (glob matches only files without "norm"; norm files use the _norm_R[12] pattern)
# Or more selectively:
#   find rnaseq_reads -maxdepth 1 -name '*_R[12].fastq.gz' ! -name '*_norm_*' -delete

set -euo pipefail

module load fastp
module load BBTools

READDIR="${1:-rnaseq_reads}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB=$(( (SLURM_MEM_PER_NODE / 1024) - 4 ))   # leave 4 GB headroom for the OS
TMPDIR="${SCRATCH:-/tmp}"

# ── Build ordered list of unnormalized R1 files ───────────────────────────────
mapfile -t R1_FILES < <(
    ls "${READDIR}"/*_R1.fastq.gz 2>/dev/null \
    | grep -v '_norm_R1\.fastq\.gz$' \
    | sort
)

NTASKS="${#R1_FILES[@]}"
if [[ "${NTASKS}" -eq 0 ]]; then
    echo "[INFO] No unnormalized R1 files found in ${READDIR}; nothing to do."
    exit 0
fi

IDX="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "${IDX}" -ge "${NTASKS}" ]]; then
    echo "[INFO] Array index ${IDX} >= ${NTASKS} files; no work for this task."
    exit 0
fi

R1="${R1_FILES[${IDX}]}"
R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

if [[ ! -f "${R2}" ]]; then
    echo "[ERROR] Expected R2 not found: ${R2}" >&2
    exit 1
fi

# Derive the species tag (everything before _R1.fastq.gz, basename only).
BASENAME=$(basename "${R1}" _R1.fastq.gz)
TEMP_NORM_R1="${TMPDIR}/${BASENAME}_norm_R1.fastq.gz"
TEMP_NORM_R2="${TMPDIR}/${BASENAME}_norm_R2.fastq.gz"
NORM_R1="${READDIR}/${BASENAME}_norm_R1.fastq.gz"
NORM_R2="${READDIR}/${BASENAME}_norm_R2.fastq.gz"

if [[ -s "${NORM_R1}" && -s "${NORM_R2}" ]]; then
    echo "[INFO] Normalized reads already exist for ${BASENAME}; skipping."
    exit 0
fi

echo "[INFO] Processing: ${BASENAME}"
echo "[INFO]   R1: ${R1}"
echo "[INFO]   R2: ${R2}"

# ── Step 1: bbnorm coverage normalization ────────────────────────────────────
echo "[INFO] Running bbnorm for ${BASENAME}..."
bbnorm.sh \
    in="${R1}"  in2="${R2}" \
    out="${TEMP_NORM_R1}" out2="${TEMP_NORM_R2}" \
    target=30 \
    threads="${CPUS}" -Xmx${MEM_GB}g
FASTP_HTML="${TMPDIR}/${BASENAME}_fastp.html"
FASTP_JSON="${TMPDIR}/${BASENAME}_fastp.json"

# ── Step 1: fastp adapter trimming and quality filtering ─────────────────────
echo "[INFO] Running fastp for ${BASENAME}..."
fastp \
    --in1  "${TEMP_NORM_R1}"      --in2  "${TEMP_NORM_R2}" \
    --out1 "${NORM_R1}" --out2 "${NORM_R2}" \
    --thread "${CPUS}" \
    --detect_adapter_for_pe \
    --cut_front  --cut_front_window_size  1 --cut_front_mean_quality  5 \
    --cut_tail   --cut_tail_window_size   1 --cut_tail_mean_quality   5 \
    --cut_right  --cut_right_window_size  4 --cut_right_mean_quality  5 \
    --length_required 25 \
    --html "${FASTP_HTML}" --json "${FASTP_JSON}"

# ── Cleanup trimmed intermediates ────────────────────────────────────────────
rm -f "${TEMP_NORM_R1}" "${TEMP_NORM_R2}" "${FASTP_HTML}" "${FASTP_JSON}"

NPAIRS=$(zcat "${NORM_R1}" 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0)
echo "[INFO] Done: ${NORM_R1} / ${NORM_R2}  (${NPAIRS} read pairs)"

# ── Cleanup note ─────────────────────────────────────────────────────────────
# Once all array tasks complete and norm outputs are verified, the original
# unnormalized files can be removed:
#
#   find rnaseq_reads -maxdepth 1 -name '*_R[12].fastq.gz' ! -name '*_norm_*' -delete
#
# Verify first with a dry-run:
#   find rnaseq_reads -maxdepth 1 -name '*_R[12].fastq.gz' ! -name '*_norm_*'
