include { tablesDir } from '../../common/utils.nf'

process MERGE_SAMPLES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path(samples)
    path(matched)

    output:
    path "samples.csv.gz", emit: samples
    path "species.csv.gz", emit: species

    script:
    """
    python3 ${params.scripts}/subset_samples.py \\
        --samples ${samples} \\
        --matched ${matched} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${params.scripts}/build_species_table.py \\
        --samples ${samples} \\
        --matched ${matched} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    """

    stub:
    """
    python3 ${params.scripts}/subset_samples.py \\
        --samples ${samples} \\
        --matched ${matched} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${params.scripts}/build_species_table.py \\
        --samples ${samples} \\
        --matched ${matched} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    """
}
