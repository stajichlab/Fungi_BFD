process PHYLING_TREE {
    tag   "${markerset}"
    label 'phyling_tree'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(filter_dir)

    output:
        path("tree"), emit: tree

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    module load phyling
    phyling tree \\
        -I ${filter_dir} \\
        -M ${params.tree_method} \\
        --concat \\
        --partition \\
        -o tree \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p tree
    echo '(Taxon_A:0.1,Taxon_B:0.1,(Taxon_C:0.05,Taxon_D:0.05):0.1);' > tree/concat.treefile
    touch tree/concat.partition
    """
}
