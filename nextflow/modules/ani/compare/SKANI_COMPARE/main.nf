include { skaniPresetFlag; skaniCpusFor; skaniMemoryFor; aniTimeFor } from '../../../common/utils.nf'

process SKANI_COMPARE {
    tag   "${group_name} [${batch_tag}] n=${n_genomes}"
    label 'skani'

    cpus   { skaniCpusFor(n_genomes) }
    memory { skaniMemoryFor(n_genomes, task.attempt) }
    time   { aniTimeFor(n_genomes, task.attempt) }

    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        // batch_tag distinguishes multiple invocations for the same group_name
        // ("full" = whole group, "cN" = prefilter component, "bI_bJ" = batch pair).
        // mode:
        //   'triangle' — query == ref (whole group / component / diagonal batch);
        //                sketches the list once, then skani triangle (upper
        //                triangle, self-pairs stripped downstream).
        //   'dist'     — query != ref (off-diagonal batch pair); skani dist
        //                ql x rl, which sketches on the fly and keeps per-job
        //                memory proportional to the batch, not the whole group.
        tuple val(group_name), val(n_genomes),
              path(query_genomes, stageAs: 'query/*'),
              path(ref_genomes,   stageAs: 'ref/*'),
              val(batch_tag), val(mode)

    output:
        tuple val(group_name), path("${group_name}.${batch_tag}.ani.tsv")

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    def compare
    if (mode == 'dist') {
        compare = """
        ls query/* > query_list.txt
        ls ref/*   > ref_list.txt
        skani dist ${preset} ${cflag} \\
            --ql query_list.txt --rl ref_list.txt \\
            --min-af ${params.skani_min_af} -t ${task.cpus} \\
            -o skani_raw.tsv
        """
    } else {
        compare = """
        ls query/* > genome_list.txt
        skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sketch_db
        skani triangle -l genome_list.txt -E \\
            --min-af ${params.skani_min_af} -t ${task.cpus} \\
            -o skani_raw.tsv
        """
    }
    """
    ${compare}
    # Normalize to name1<TAB>name2<TAB>ANI; strip header & self-comparisons.
    # (skani dist header is Ref_file<Query_file<ANI<...>; triangle header is
    # ref-genome<query-genome<ANI<...> -- both put the pair in \$1/\$2 and ANI
    # in \$3, so one filter works for either.)
    awk -F'\\t' 'NR>1 && \$1 != \$2 {
        print \$1"\\t"\$2"\\t"\$3
    }' skani_raw.tsv > ${group_name}.${batch_tag}.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.${batch_tag}.ani.tsv
    """
}
