process PHYLING_FILTER {
    tag   "${markerset}"
    label 'phyling_filter'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(align_dir)

    output:
        tuple val(markerset), val(taxon_slug), path("filter"), emit: filtered

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    module load phyling
    phyling filter \\
        -I ${align_dir} \\
        -n ${params.top_n} \\
        -o filter \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p filter
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > filter/BUSCOmarker001.fa
    """
}
