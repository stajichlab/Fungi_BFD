// ── skani triangle (sparse) over cached sketches ─────────────────────────────
process SKANI_COMPARE {
    tag   "${group_name} n=${sketches.size()}"
    label 'skani'

    cpus   { sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8 }
    memory { sketches.size() > 500 ? '64 GB' : sketches.size() > 200 ? '32 GB' : '16 GB' }

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
    ls *.sketch > sketch_list.txt
    skani triangle --sparse -l sketch_list.txt \\
        --min-af ${params.skani_min_af} -t ${task.cpus} -o skani_raw.tsv

    # skani cols: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query [names]
    awk -F'\\t' 'NR>1 {
        nr=split(\$1,b,"/"); r=b[nr]; sub(/\\.sketch\$/,"",r);
        nq=split(\$2,a,"/"); q=a[nq]; sub(/\\.sketch\$/,"",q);
        if (q!=r) print q"\\t"r"\\t"\$3
    }' skani_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}
