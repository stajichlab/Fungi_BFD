process SOURMASH_SKETCH {
    tag   "${genome.name}"
    label 'sourmash'
    cpus   2
    memory '4 GB'

    storeDir "${params.sketch_cache}/sourmash/k${params.sourmash_kmer}_scaled${params.sourmash_scaled}"

    input:
        tuple val(group_name), path(genome)

    output:
        tuple val(group_name), path("${genome.name}.sig.zip")

    script:
    """
    sourmash sketch dna -p k=${params.sourmash_kmer},scaled=${params.sourmash_scaled} \\
        --name "${genome.name}" ${genome} -o ${genome.name}.sig.zip
    """

    stub:
    """
    touch ${genome.name}.sig.zip
    """
}
