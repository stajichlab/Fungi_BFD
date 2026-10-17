// DRAFT. One chunk of contigs predicted with a fixed GeneMark model (scatter step).
// Same gmes_petap.pl call as GENEMARK_TRAIN's fast-reuse branch (--predict_with with the same
// --max_intron / --soft_mask 2000 / --fungus), so a chunk gives the same genes the monolithic
// run would for those contigs. --soft_mask must stay numeric: 'auto' derives the masking
// threshold from input file size, which differs between a chunk and the whole genome.
//
// The production shell is `/bin/bash -l` (nextflow.config process.shell), which has NO `set -e`,
// so this script sets its own strict mode. A chunk passes only when all of these hold:
//   - every contig piece gmes made (output/data/dna.fa_N) has a result file (output/gmhmm/dna.fa_N.out)
//   - no "warning on: ... gmhmme3" line in gmes.log
//   - every contig name in the GTF is a contig of this chunk
// A chunk with no genes at all is accepted only when those checks pass (an empty GTF is then correct).
process GENEMARK_PREDICT_CHUNK {
    tag "${out}:${idx}"

    cpus   8
    memory '8 GB'
    time   '2h'

    input:
    tuple val(out), val(idx), val(n_chunks), path(chunk_fa), path(model), val(transl_table)

    output:
    // idx is the zero-padded chunk number; GENEMARK_MERGE checks it, never relies on arrival order
    tuple val(out), val(idx), val(n_chunks), path("${out}.chunk_${idx}.gtf"), emit: gtf

    script:
    // gmes_petap.pl (CheckGcode) accepts only --gcode 1, 6 and 29; table 26 gives
    // "error, input genetic code 26 is not supported" (exit 1). Table 26 (CUG=Ala)
    // keeps the standard stop codons, so GeneMark runs with its default code 1.
    def gm_table = (transl_table.toString() == '26') ? '1' : transl_table
    """
    set -euo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer
    # Node-local scratch is required: a chunk makes ~20k small contig files. SCRATCH is set inside tasks.
    WORK=\$(mktemp -d "\${SCRATCH:?SCRATCH not set}/gmchunk.XXXXXX")
    trap 'rm -rf "\$WORK"' EXIT
    echo "[INFO] GENEMARK_PREDICT_CHUNK ${out}:${idx} work dir \$WORK (SCRATCH='\${SCRATCH}' host=\$(hostname))" >&2
    export TMPDIR="\$SCRATCH"

    GM_KEY_REAL=\$(readlink -f "\$HOME/.gm_key" 2>/dev/null || echo "\$HOME/.gm_key")
    GM_KEY_DIR=\$(dirname "\$GM_KEY_REAL")
    SING="apptainer exec --bind \$WORK:\$WORK,\$PWD:\$PWD,${workflow.projectDir}:${workflow.projectDir},\$GM_KEY_DIR:\$GM_KEY_DIR,\$TMPDIR:\$TMPDIR ${params.genemark_sif}"

    pigz -dc "${chunk_fa}" > "\$WORK/genome.fa"
    [ -s "\$WORK/genome.fa" ] || { echo "ERROR: chunk ${chunk_fa} inflated to an empty file" >&2; exit 1; }
    [ -s "${model}" ] || { echo "ERROR: model ${model} is missing or empty" >&2; exit 1; }
    cp "${model}" "\$WORK/genemark-shared.mod"

    # Same --gcode decision as GENEMARK_TRAIN. Capture the help text first: piping gmes into
    # `grep -q` under pipefail can fail on SIGPIPE. An empty help text means the probe itself failed.
    GCODE_ARGS=()
    if [ "${gm_table}" = "6" ] || [ "${gm_table}" = "26" ]; then
        GMES_HELP=\$(\$SING gmes_petap.pl 2>&1 || true)
        [ -n "\$GMES_HELP" ] || { echo "ERROR: could not run gmes_petap.pl to probe --gcode support (transl_table ${transl_table})" >&2; exit 1; }
        if printf '%s' "\$GMES_HELP" | grep -q -- '--gcode'; then GCODE_ARGS=(--gcode "${gm_table}"); fi
    fi

    CHUNK_CONTIGS=\$(grep -c '^>' "\$WORK/genome.fa")
    CHUNK_BP=\$(grep -v '^>' "\$WORK/genome.fa" | tr -d '\n' | wc -c)
    T_START=\$SECONDS
    cd "\$WORK"
    # gmes_petap.pl exits 1 on its own errors; its exit status does not tell us about per-contig
    # gmhmme3 failures, so check artifacts below instead.
    \$SING gmes_petap.pl --predict_with genemark-shared.mod \\
        --sequence genome.fa --max_intron ${params.max_intronlen} --soft_mask 2000 \\
        --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" > gmes_stdout.log 2>&1 || true
    cd "\$OLDPWD"
    echo "[INFO] GENEMARK_PREDICT_CHUNK ${out}:${idx} contigs=\$CHUNK_CONTIGS bp=\$CHUNK_BP gmes_wall_s=\$((SECONDS - T_START)) cpus=${task.cpus}" >&2
    cp "\$WORK/gmes_stdout.log" "${out}.chunk_${idx}.gmes_stdout.log"
    tail -n 20 "\$WORK/gmes_stdout.log" >&2

    # 1. completeness: one result file per contig piece
    PIECES=\$(ls "\$WORK"/output/data/dna.fa_* 2>/dev/null | wc -l)
    RESULTS=\$(ls "\$WORK"/output/gmhmm/dna.fa_*.out 2>/dev/null | wc -l)
    if [ "\$PIECES" -eq 0 ] || [ "\$PIECES" -ne "\$RESULTS" ]; then
        echo "ERROR: GENEMARK_PREDICT_CHUNK ${out}:${idx}: \$PIECES contig pieces but \$RESULTS result files" >&2
        exit 1
    fi
    # 2. no gmhmme3 failures
    FAILS=\$(grep -c "warning on: .*gmhmme3" "\$WORK/gmes.log" 2>/dev/null || true)
    if [ "\${FAILS:-0}" -gt 0 ]; then
        echo "ERROR: GENEMARK_PREDICT_CHUNK ${out}:${idx}: gmhmme3 failed \$FAILS time(s); likely out of memory at ${task.memory}" >&2
        exit 1
    fi
    # 3. the final GTF must exist (empty is fine when checks 1 and 2 passed: a chunk with no genes)
    [ -e "\$WORK/genemark.gtf" ] || { echo "ERROR: GENEMARK_PREDICT_CHUNK ${out}:${idx}: gmes wrote no genemark.gtf" >&2; exit 1; }
    if [ ! -s "\$WORK/genemark.gtf" ]; then
        echo "[INFO] GENEMARK_PREDICT_CHUNK ${out}:${idx}: no genes predicted in \$PIECES contig pieces (all results present, no warnings)" >&2
    fi
    # 4. contig names must belong to this chunk (empty contig list = failure, not a pass)
    awk '/^>/ { print substr(\$1, 2) }' "\$WORK/genome.fa" > "\$WORK/contig_names.txt"
    [ -s "\$WORK/contig_names.txt" ] || { echo "ERROR: no contig names in chunk FASTA" >&2; exit 1; }
    if ! awk 'NR == FNR { ids[\$1] = 1; next }
              !/^#/ && NF >= 9 && !(\$1 in ids) { print "unknown sequence " \$1 > "/dev/stderr"; bad = 1; exit 1 }
              END { exit bad }' "\$WORK/contig_names.txt" "\$WORK/genemark.gtf"; then
        echo "ERROR: GENEMARK_PREDICT_CHUNK ${out}:${idx}: GTF names sequences that are not in the chunk" >&2
        exit 1
    fi
    cp "\$WORK/genemark.gtf" "${out}.chunk_${idx}.gtf"
    """

    stub:
    """
    touch "${out}.chunk_${idx}.gtf"
    """
}
