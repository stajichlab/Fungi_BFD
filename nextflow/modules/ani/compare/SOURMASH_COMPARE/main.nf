// ── sourmash compare (ANI matrix) over cached signatures ─────────────────────
process SOURMASH_COMPARE {
    tag   "${group_name} n=${sigs.size()}"
    label 'sourmash'

    cpus   { sigs.size() > 500 ? 32 : sigs.size() > 200 ? 16 : 8 }
    memory { sigs.size() > 500 ? '128 GB' : sigs.size() > 200 ? '64 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(sigs)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    ls *.sig.zip > sig_list.txt
    sourmash compare --ani --containment --quiet \\
        --from-file sig_list.txt --csv cmp.csv -k ${params.sourmash_kmer}

    python3 ${projectDir}/bin/sourmash_matrix_to_long.py --input cmp.csv --output ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}
