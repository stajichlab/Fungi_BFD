process BUSCO_PEP {
    tag        "${locustag}"
    label      'busco_pep'
    storeDir   "${params.genome_stats_outdir}/BUSCO_protein"

    input:
        tuple val(locustag), val(basename), val(lineage), path(proteins)

    output:
        path("${basename}.BUSCO_summary.${lineage}.txt"), emit: summary

    script:
    """
    module load busco
    export BUSCO_LINEAGES=/srv/projects/db/BUSCO/v12/
    RUNDIR="${basename}_busco_pep"
    busco -i ${proteins} \\
          -l ${lineage} \\
          -m proteins \\
          -o \$RUNDIR \\
          -c ${task.cpus} \\
          --offline --download_path \$BUSCO_LINEAGES
    SUMMARY=\$(ls \${RUNDIR}/short_summary.specific.*.txt 2>/dev/null | head -1)
    [ -z "\$SUMMARY" ] && SUMMARY=\$(ls \${RUNDIR}/short_summary*.txt 2>/dev/null | head -1)
    if [ -z "\$SUMMARY" ]; then
        echo "[ERROR] BUSCO proteins: no short_summary file found for ${basename}" >&2
        exit 1
    fi
    cp "\$SUMMARY" "${basename}.BUSCO_summary.${lineage}.txt"
    """

    stub:
    """
    printf '# BUSCO version: 5.x\\n# The lineage dataset is: ${lineage}\\n\\tC:98.0%%[S:97.0%%,D:1.0%%],F:1.0%%,M:1.0%%,n:758\\n' \\
        > "${basename}.BUSCO_summary.${lineage}.txt"
    """
}
