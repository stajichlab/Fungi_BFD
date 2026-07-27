include { skaniPresetFlag } from '../../../common/utils.nf'

// skani 0.3.x outputs a sketch directory (sketches.db + index.db + markers.bin)
// per chunk instead of individual .sketch files. SKANI_COMPARE passes the chunk
// directories directly to skani dist via -rl/-ql.
process SKANI_SKETCH {
    tag   "${group_name} [${genomes.size()} genomes]"
    label 'skani'
    cpus   4
    memory '8 GB'

    storeDir "${params.sketch_cache}/skani/${params.skani_preset}_af${params.skani_min_af}"

    input:
        tuple val(group_name), path(genomes), val(sketch_dir)

    output:
        tuple val(group_name), path(sketch_dir)

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    """
    # skani sketch errors out if its -o dir already exists; clear any partial
    # output left behind by an interrupted/retried run before re-sketching.
    rm -rf "${sketch_dir}"
    printf '%s\\n' ${genomes} > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o "${sketch_dir}"
    """

    stub:
    """
    mkdir -p "${sketch_dir}"
    touch "${sketch_dir}/sketches.db" "${sketch_dir}/index.db" "${sketch_dir}/markers.bin"
    """
}
