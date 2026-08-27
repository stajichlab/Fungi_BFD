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
        // (one per skani_prefilter component, "cN"); plain whole-group calls
        // pass "full", matching the historical filename so existing outputs/
        // downstream globs are unaffected when skani_prefilter is off.
        tuple val(group_name), val(n_genomes), path(genomes), val(batch_tag)

    output:
        tuple val(group_name), path("${group_name}.${batch_tag}.ani.tsv")

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    def genome_list = genomes.collect { it.toString() }.join('\n')
    """
    # skani 0.3.x: sketch everything into one directory, then skani triangle
    # on that directory. No per-genome .sketch files, no chunk merge needed.
    printf '%s\\n' "${genome_list}" > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sketch_db

    # skani triangle with -E (sparse) outputs the same row-by-row format as
    # skani dist: name1 <TAB> name2 <TAB> ANI <TAB> AF_ref <TAB> AF_query
    skani triangle -l genome_list.txt -E \\
        --min-af ${params.skani_min_af} -t ${task.cpus} \\
        -o skani_raw.tsv

    # Strip self-comparisons (diagonal) and extra columns
    awk -F'\\t' 'NR>1 && \$1 != \$2 {
        print \$1"\\t"\$2"\\t"\$3
    }' skani_raw.tsv > ${group_name}.${batch_tag}.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.${batch_tag}.ani.tsv
    """
}
