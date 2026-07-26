process CALC_INTERGENIC {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/intergenic_stats"

    input:
    tuple val(locustag), val(basename), path(gff_file)

    output:
    path "${basename}.gene_intergenic_distances.csv.gz", emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_intergenic.py \\
        ${gff_file} -o .
    pigz gene_pairwise_distances.csv
    mv gene_pairwise_distances.csv.gz ${basename}.gene_intergenic_distances.csv.gz
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > ${basename}.gene_intergenic_distances.csv.gz
    """
}
