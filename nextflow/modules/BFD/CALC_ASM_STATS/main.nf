include { hashBucketForType } from '../../common/utils.nf'

process CALC_ASM_STATS {
    label    'asmstats'
    tag      "${meta.asmid}"
    storeDir { "${params.genome_stats_outdir}/asm_stats/${hashBucketForType('asm_stats', meta.asmid)}" }

    input:
    tuple val(meta), path(genome)

    output:
    path "${meta.asmid}.stats.txt", emit: stats

    script:
    """
    module load singularity
    singularity exec ${params.aaftf_sif} \\
        AAFTF assess -i ${genome} -r ${meta.asmid}.stats.txt
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
