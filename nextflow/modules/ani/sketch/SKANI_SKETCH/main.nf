include { skaniPresetFlag } from '../../../common/utils.nf'

// skani sketching is fast, so genomes are batched --skani_sketch_chunk per job
// (instead of one SLURM job per genome) to cut scheduler churn. Output filenames
// are declared deterministically (one <genome>.sketch per input) so storeDir
// still caches each chunk and reuses it across runs / --compare levels.
process SKANI_SKETCH {
    tag   "${group_name} [${genomes.size()} genomes]"
    label 'skani'
    cpus   4
    memory '8 GB'

    storeDir "${params.sketch_cache}/skani/${params.skani_preset}_af${params.skani_min_af}"

    input:
        tuple val(group_name), path(genomes), val(sketch_names)

    output:
        tuple val(group_name), path(sketch_names)

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    """
    # skani sketch errors out if its -o dir already exists; clear any partial
    # output left behind by an interrupted/retried run before re-sketching.
    rm -rf sk_out
    printf '%s\\n' ${genomes} > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sk_out
    mv sk_out/*.sketch .
    """

    stub:
    """
    for s in ${sketch_names.join(' ')}; do touch \$s; done
    """
}
