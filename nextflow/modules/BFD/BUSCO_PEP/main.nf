include { hashBucketForType } from '../../common/utils.nf'

process BUSCO_PEP {
    tag        "${meta.locustag}"
    label      'busco_pep'
    storeDir   { "${params.genome_stats_outdir}/BUSCO_protein/${hashBucketForType('BUSCO_protein', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.BUSCO_summary.${meta.lineage}.txt"), emit: summary

    script:
    """
    module load busco
    export BUSCO_LINEAGES=/srv/projects/db/BUSCO/v12/
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    SCRATCH=\${SCRATCH:-/tmp}
    RUNDIR="${meta.locustag}_busco_pep"
    busco -i ${proteins} \\
          -l ${meta.lineage} \\
          -m proteins \\
          -o \$RUNDIR \\
          --out_path \$SCRATCH \\
          -c ${task.cpus} \\
          --offline --download_path \$BUSCO_LINEAGES
    SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary.specific.*.txt 2>/dev/null | head -1)
    [ -z "\$SUMMARY" ] && SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary*.txt 2>/dev/null | head -1)
    if [ -z "\$SUMMARY" ]; then
        echo "[ERROR] BUSCO proteins: no short_summary file found for ${meta.locustag}" >&2
        exit 1
    fi
    cp "\$SUMMARY" "${meta.locustag}.BUSCO_summary.${meta.lineage}.txt"
    rm -rf "\$SCRATCH/\$RUNDIR"
    """

    stub:
    """
    printf '# BUSCO version: 5.x\\n# The lineage dataset is: ${meta.lineage}\\n\\tC:98.0%%[S:97.0%%,D:1.0%%],F:1.0%%,M:1.0%%,n:758\\n' \\
        > "${meta.locustag}.BUSCO_summary.${meta.lineage}.txt"
    """
}
