process ORTHOFINDER_PARSE {
    label      'comparative_orthofinder'
    tag        "parse"
    publishDir "${params.outdir}/${params.project}/orthofinder", mode: 'copy'

    input:
    path orthofinder_out

    output:
    path "Orthogroups.tsv",                    emit: orthogroups
    path "Orthogroups_UnassignedGenes.tsv",     emit: unassigned, optional: true
    path "Statistics_Overall.tsv",              emit: stats,      optional: true

    script:
    """
    OF_DIR="${orthofinder_out}"
    if ls "\${OF_DIR}"/Results_* >/dev/null 2>&1; then
        OF_DIR=\$(ls -d "\${OF_DIR}"/Results_* | sort | tail -1)
    fi
    cp "\${OF_DIR}/Orthogroups/Orthogroups.tsv" .
    [ -f "\${OF_DIR}/Orthogroups/Orthogroups_UnassignedGenes.tsv" ] && \\
        cp "\${OF_DIR}/Orthogroups/Orthogroups_UnassignedGenes.tsv" . || true
    [ -f "\${OF_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" ] && \\
        cp "\${OF_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" . || true
    """

    stub:
    """
    touch Orthogroups.tsv Statistics_Overall.tsv
    """
}
