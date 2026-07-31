include { hashBucketForType } from '../../common/utils.nf'

process BUSCO_GENOME {
    tag        "${meta.asmid}"
    label      'busco_genome'
    storeDir   { "${params.genome_stats_outdir}/BUSCO_genome/${hashBucketForType('BUSCO_genome', meta.asmid)}" }

    input:
        tuple val(meta), path(genome)

    output:
        path("${meta.asmid}.BUSCO_summary.${meta.lineage}.txt"), emit: summary

    script:
    """
    module load busco
    export BUSCO_LINEAGES=/srv/projects/db/BUSCO/v12/
    RUNDIR="${meta.asmid}_busco_genome"
    busco -i ${genome} \\
          -l ${meta.lineage} \\
          -m genome \\
          -o \$RUNDIR \\
          -c ${task.cpus} \\
          --offline --download_path \$BUSCO_LINEAGES
    SUMMARY=\$(ls \${RUNDIR}/short_summary.specific.*.txt 2>/dev/null | head -1)
    [ -z "\$SUMMARY" ] && SUMMARY=\$(ls \${RUNDIR}/short_summary*.txt 2>/dev/null | head -1)
    if [ -z "\$SUMMARY" ]; then
        echo "[ERROR] BUSCO genome: no short_summary file found for ${meta.asmid}" >&2
        exit 1
    fi
    cp "\$SUMMARY" "${meta.asmid}.BUSCO_summary.${meta.lineage}.txt"
    rm -rf "\$RUNDIR"
    """

    stub:
    """
    printf '# BUSCO version: 5.x\\n# The lineage dataset is: ${meta.lineage}\\n\\tC:99.0%%[S:98.0%%,D:1.0%%],F:0.5%%,M:0.5%%,n:758\\n' \\
        > "${meta.asmid}.BUSCO_summary.${meta.lineage}.txt"
    """
}
