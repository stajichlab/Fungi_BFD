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
    // BBTools 39.x (bundled in busco/6.1.0) treats any stats.sh argument that
    // contains "xms"/"xmx" (case-insensitive, anywhere in the string) as a Java
    // heap flag. BUSCO passes the absolute input path to stats.sh, so a strain ID
    // such as LIQ201XMSP breaks the JVM launch. Stage the genome under a fixed
    // name that keeps the compression suffix.
    def safe_genome = genome.name.endsWith('.gz') ? 'genome.fa.gz' : 'genome.fa'
    """
    module load busco
    ln -sf ${genome} ${safe_genome}
    export BUSCO_LINEAGES=${params.busco_lineages_dir}
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    : \${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    RUNDIR="${meta.asmid}_busco_genome"
    busco -i ${safe_genome} \\
          -l ${meta.lineage} \\
          -m genome \\
          -o \$RUNDIR \\
          --out_path \$SCRATCH \\
          -c ${task.cpus} \\
          --offline --download_path \$BUSCO_LINEAGES
    SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary.specific.*.txt 2>/dev/null | head -1)
    [ -z "\$SUMMARY" ] && SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary*.txt 2>/dev/null | head -1)
    if [ -z "\$SUMMARY" ]; then
        echo "[ERROR] BUSCO genome: no short_summary file found for ${meta.asmid}" >&2
        exit 1
    fi
    cp "\$SUMMARY" "${meta.asmid}.BUSCO_summary.${meta.lineage}.txt"
    rm -rf "\$SCRATCH/\$RUNDIR"
    """

    stub:
    """
    printf '# BUSCO version: 5.x\\n# The lineage dataset is: ${meta.lineage}\\n\\tC:99.0%%[S:98.0%%,D:1.0%%],F:0.5%%,M:0.5%%,n:758\\n' \\
        > "${meta.asmid}.BUSCO_summary.${meta.lineage}.txt"
    """
}
