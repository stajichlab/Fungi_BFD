// ── mash prefilter: paste sketches, all-vs-all dist (cheap) ───────────────────
include { aniTimeFor; capCpus; capMemGB } from '../../../common/utils.nf'

process MASH_PREFILTER {
    tag   "${group_name} n=${sketches.size()}"
    label 'mash'

    cpus   { capCpus(sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8) }
    memory { capMemGB(sketches.size() > 500 ? 32 : 16).toString() + ' GB' }
    time   { aniTimeFor(sketches.size(), task.attempt) }

    input:
        tuple val(group_name), path(sketches)

    output:
        tuple val(group_name), path("${group_name}.mash_dist.tsv")

    script:
    """
    ls *.msh > sketch_list.txt
    mash paste combined -l sketch_list.txt
    mash dist -p ${task.cpus} combined.msh combined.msh > ${group_name}.mash_dist.tsv
    """

    stub:
    """
    touch ${group_name}.mash_dist.tsv
    """
}
