include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

process CALC_ASM_STATS {
    label    'asmstats'
    tag      "${meta.asmid}"
    storeDir { "${params.genome_stats_outdir}/asm_stats/${hashBucketForType('asm_stats', meta.asmid)}" }

    input:
    tuple val(meta), path(genome)

    output:
    path "${meta.asmid}.stats.txt", emit: stats

    script:
    // The genome may be <id>.scaffolds.fa.gz (bgzip, FUNANNOTATE_PREDICT since 2026-11-07).
    // Inflate it into the task dir: this tool may not read gzip.
    def genome_plain = plainName(genome)
    def genome_rm = genome.name.endsWith('.gz') ? "rm -f '${genome_plain}'" : ''
    """
    ${gzInflate(genome)}
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    # ${genome_plain} may be a symlink Nextflow staged from outside \$PWD (e.g. a
    # shared genome_dir elsewhere under /bigdata); bind its real parent dir
    # too, or apptainer won't see the symlink target inside the container.
    GENOME_REAL_DIR=\$(dirname \$(readlink -f ${genome_plain}))
    SING_BINDS="--bind \${PWD}:\${PWD},\${GENOME_REAL_DIR}:\${GENOME_REAL_DIR},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.aaftf_sif} \\
        AAFTF assess -i ${genome_plain} -r ${meta.asmid}.stats.txt
    # Read by bin/summarize_asm_stats.py; a stats file without this line is annotated.
    echo "ANNOTATION STATUS = ${meta.annotated == false ? 'unannotated' : 'annotated'}" >> ${meta.asmid}.stats.txt
    ${genome_rm}
    """

    stub:
    """
    cat > ${meta.asmid}.stats.txt <<'EOF'
Assembly statistics for: ${meta.id}.scaffolds.fa
   CONTIG COUNT  =  10
   TOTAL LENGTH  =  1000000
            MIN  =  500
            MAX  =  200000
         MEDIAN  =  50000
           MEAN  =  100000.0
            L50  =  3
            N50  =  150000
            L90  =  8
            N90  =  20000
            GC%  =  50.00
    N GAP COUNT  =  0
  TOTAL N BASES  =  0
   BASES MASKED  =  10000
 PERCENT MASKED  =  1.00
  T2T SCAFFOLDS  =  0
   TELOMERE FWD  =  0
   TELOMERE REV  =  0
EOF
    """
}
