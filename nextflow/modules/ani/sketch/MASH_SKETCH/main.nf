process MASH_SKETCH {
    tag   "${genome.name}"
    label 'mash'
    cpus   2
    memory '4 GB'

    storeDir "${params.sketch_cache}/mash/k${params.mash_kmer}_s${params.mash_sketch_size}"

    input:
        tuple val(group_name), path(genome)

    output:
        tuple val(group_name), path("${genome.baseName}.msh")

    script:
    """
    mash sketch -k ${params.mash_kmer} -s ${params.mash_sketch_size} \\
        -o ${genome.baseName} ${genome}
    """

    stub:
    """
    touch ${genome.baseName}.msh
    """
}
