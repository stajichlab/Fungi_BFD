process REPORT_QUERY_ANI {
    tag   "${group_name}"
    label 'report'

    publishDir { "${params.outdir}/skani_query/${params.compare}/${group_name}" }, mode: 'copy'

    input:
        tuple val(group_name), path(ani_tsv), path(names_tsv)

    output:
        path("${group_name}_query_report.txt")
        path("${group_name}_query_calls.tsv")

    script:
    """
    python3 ${projectDir}/bin/report_query_ani.py \\
        --input  ${ani_tsv} \\
        --names  ${names_tsv} \\
        --group  "${group_name}" \\
        --level  "${params.compare}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --outlier-threshold ${params.ani_outlier_threshold} \\
        --top-n  ${params.query_top_n} \\
        --report ${group_name}_query_report.txt \\
        --calls  ${group_name}_query_calls.tsv
    """

    stub:
    """
    printf '=== ${group_name} (stub) ===\\n' > ${group_name}_query_report.txt
    printf 'group\\tlevel\\tquery_asmid\\tquery_species\\tquery_strain\\tbest_ref_asmid\\tbest_ref_genus\\tbest_ref_species\\tbest_ani\\tn_refs_ge_outlier\\ttier\\n' > ${group_name}_query_calls.tsv
    """
}
