process REPORT_ANI {
    tag   "${group_name}"
    label 'report'

    // storeDir removed (2026-08-26): it sat alongside publishDir pointed at the
    // same path, so its bare-existence check silently won -- a group whose
    // membership changed (taxonomy reclassification, new genome added) never
    // regenerated its report as long as *a* report already existed there,
    // regardless of -resume. publishDir + Nextflow's normal input-hash caching
    // (already keyed on ani_tsv/names_tsv, which do change when membership
    // changes) is the correct behavior -- same fix as MERGE_ANI and the
    // pairwise *_COMPARE modules. See
    // nextflow/docs/DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md, gate inventory
    // item 6 / Option 2 extension (KCTC_13826BP/MRD-KRBAY incident, 2026-08-26).
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
