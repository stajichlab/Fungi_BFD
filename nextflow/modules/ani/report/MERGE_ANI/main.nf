// Concatenate one-or-many partial TSVs into the group's full ANI table.
process MERGE_ANI {
    tag   "${group_name}"
    label 'report'

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}"

    input:
        tuple val(group_name), path(part_tsvs)

    output:
        tuple val(group_name), path("${group_name}.ani.tsv")

    script:
    """
    cat ${part_tsvs} > ${group_name}.ani.tsv
    """

    stub:
    """
    cat ${part_tsvs} > ${group_name}.ani.tsv
    """
}
