process REPORT_ANI {
    tag   "${group_name}"
    label 'report'

    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}" }, mode: 'copy'

    input:
        tuple val(group_name), path(ani_tsv), path(names_tsv)

    output:
        path("${group_name}_ANI_report.txt"), emit: report
        path("${group_name}_genome_names.tsv"), emit: names

    script:
    """
    cp ${names_tsv} ${group_name}_genome_names.tsv
    python3 ${projectDir}/bin/report_ani.py \\
        --input    ${ani_tsv} \\
        --names    ${group_name}_genome_names.tsv \\
        --group    "${group_name}" \\
        --level    "${params.compare}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --outlier-threshold ${params.ani_outlier_threshold} \\
        --output   ${group_name}_ANI_report.txt
    """

    stub:
    """
    printf '=== ANI Report: ${group_name} (stub) ===\\n' > ${group_name}_ANI_report.txt
    cp ${names_tsv} ${group_name}_genome_names.tsv
    """
}
