// ── skani triangle (sparse) over cached sketches ─────────────────────────────
process SKANI_COMPARE {
    tag   "${group_name} n=${sketches.size()}"
    label 'skani'

    cpus   { sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8 }
    memory { sketches.size() > 500 ? '64 GB' : sketches.size() > 200 ? '32 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

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
