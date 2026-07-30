// ── mash dist over cached sketches ───────────────────────────────────────────
include { aniTimeFor; capCpus; capMemGB } from '../../../common/utils.nf'

process MASH_COMPARE {
    tag   "${group_name} n=${sketches.size()}"
    label 'mash'

    cpus   { capCpus(sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8) }
    memory { capMemGB(sketches.size() > 500 ? 32 : 16).toString() + ' GB' }
    time   { aniTimeFor(sketches.size(), task.attempt) }

    // publishDir (not storeDir): storeDir only checks output-path existence, so it
    // would never re-run a group's comparison after new genomes are added to that
    // group later -- it silently keeps serving the stale pre-addition result forever,
    // regardless of -resume. publishDir + Nextflow's normal input-hash caching
    // correctly detects a changed sketch/genome list and reruns automatically
    // (found via a genome added 2026-07-19 that never appeared in the GENUS-level
    // Aspergillus comparison even after later reruns -- see .living/learnings.md).
    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name), path(sketches)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    ls *.msh > sketch_list.txt
    mash paste combined -l sketch_list.txt
    mash dist -p ${task.cpus} combined.msh combined.msh > mash_raw.tsv

    # mash cols: ref  query  distance  p-value  shared-hashes  ->  ANI = 100*(1-dist)
    awk -F'\\t' '{
        nr=split(\$1,b,"/"); r=b[nr];
        nq=split(\$2,a,"/"); q=a[nq];
        if (q!=r) printf "%s\\t%s\\t%.4f\\n", q, r, (1-\$3)*100
    }' mash_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}
