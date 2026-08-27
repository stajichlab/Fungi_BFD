// Concatenate one-or-many partial TSVs into the group's full ANI table.
process MERGE_ANI {
    tag   "${group_name}"
    label 'report'

    // publishDir (not storeDir): storeDir only checks output-path existence, so it
    // would never re-run a group's merge after group membership changes (e.g. a
    // genome reclassified to/from this species/genus by a samples.csv taxonomy
    // fix) -- it silently keeps serving the stale pre-change merged table forever,
    // regardless of -resume. publishDir + Nextflow's normal input-hash caching
    // correctly detects a changed part_tsvs set and reruns automatically. Same
    // fix already applied to the pairwise *_COMPARE modules (see their comments);
    // this was the last storeDir left in the ANI report path (see
    // nextflow/docs/DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md, gate inventory item 6
    // / Option 2 extension -- found via the KCTC_13826BP/MRD-KRBAY misidentified-
    // representative incident, 2026-08-26).
    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}" }, mode: 'copy'

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
