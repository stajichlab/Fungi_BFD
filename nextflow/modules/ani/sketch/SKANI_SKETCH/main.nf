include { skaniPresetFlag } from '../../../common/utils.nf'

// skani 0.3.x produces a consolidated sketch database per chunk (sketches.db +
// index.db + markers.bin) instead of individual .sketch files. The chunk-level
// databases are merged with `skani sketch --merge` in SKANI_COMPARE.
process SKANI_SKETCH {
    tag   "${group_name} [${genomes.size()} genomes]"
    label 'skani'
    cpus   4
    memory '8 GB'

    storeDir "${params.sketch_cache}/skani/${params.skani_preset}_af${params.skani_min_af}"

    input:
        tuple val(group_name), path(genomes), val(sketch_db)

    output:
        tuple val(group_name), path(sketch_db)

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    """
    # skani sketch errors out if its -o dir already exists; clear any partial
    # output left behind by an interrupted/retried run before re-sketching.
    rm -rf sk_out
    printf '%s\\n' ${genomes} > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sk_out
    mv sk_out/sketches.db ${sketch_db}
    mv sk_out/index.db    ${sketch_db}.index
    mv sk_out/markers.bin ${sketch_db}.markers
    """

    stub:
    """
    touch ${sketch_db} ${sketch_db}.index ${sketch_db}.markers
    """
}
