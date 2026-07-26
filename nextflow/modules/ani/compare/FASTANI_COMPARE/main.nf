// ── fastANI (query list × ref list) ──────────────────────────────────────────
// Used both for full-group single jobs and per-component prefiltered jobs.
// query/ref staged into separate subdirs so identical files don't collide.
// group_size drives CPU/memory scaling.
process FASTANI_COMPARE {
    tag   "${group_name} [${batch_tag}] n=${group_size}"
    label 'fastani'

    cpus   { group_size > 500 ? 64 : group_size > 200 ? 24 : 8 }
    memory { group_size > 500 ? '128 GB' : group_size > 200 ? '48 GB' : '16 GB' }

    // publishDir (not storeDir): storeDir only checks output-path existence, so it
    // would never re-run a group's comparison after new genomes are added to that
    // group later -- it silently keeps serving the stale pre-addition result forever,
    // regardless of -resume. publishDir + Nextflow's normal input-hash caching
    // correctly detects a changed sketch/genome list and reruns automatically
    // (found via a genome added 2026-07-19 that never appeared in the GENUS-level
    // Aspergillus comparison even after later reruns -- see .living/learnings.md).
    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name),
              path(query_genomes, stageAs: 'query/*'),
              path(ref_genomes,   stageAs: 'ref/*'),
              val(batch_tag), val(group_size)

    output:
        tuple val(group_name), path("${group_name}.${batch_tag}.ani.tsv")

    script:
    """
    ls query/* > query_list.txt
    ls ref/*   > ref_list.txt
    fastANI \\
        --ql query_list.txt \\
        --rl ref_list.txt \\
        -o ${group_name}.${batch_tag}.ani.tsv \\
        --fragLen ${params.fastani_fraglen} \\
        -k ${params.fastani_kmer} \\
        -t ${task.cpus}
    """

    stub:
    """
    printf 'query\\treference\\t99.0\\t100\\t100\\n' > ${group_name}.${batch_tag}.ani.tsv
    """
}
