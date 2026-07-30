// ── connected components from mash dist (python, in the report env) ──────────
process MASH_COMPONENTS {
    tag   "${group_name}"
    label 'report'

    input:
        tuple val(group_name), path(mash_dist)

    output:
        tuple val(group_name), path("${group_name}.components.tsv")

    script:
    """
    python3 ${projectDir}/bin/mash_components.py --input ${mash_dist} \\
        --ani ${params.prefilter_ani} --min-size ${params.min_group_size} \\
        > ${group_name}.components.tsv
    """

    stub:
    """
    touch ${group_name}.components.tsv
    """
}
