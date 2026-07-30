process PHYLING_ALIGN {
    tag   "${markerset}"
    label 'phyling_align'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(markerset_dir), path(fastas)

    output:
        tuple val(markerset), val(taxon_slug), path("align"), emit: align

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    mkdir -p staged
    for f in *.fa *.fa.gz; do
        [ -e "\$f" ] || continue
        base="\$(basename "\$f" .gz)"
        base="\$(basename "\$base" .fa)"
        base="\${base%.proteins}"
        base="\${base%.cds-transcripts}"
        ext="\$(echo "\$f" | grep -oE '\\.fa(\\.gz)?\$')"
        ln -sfn "\$(readlink -f "\$f")" "staged/\${base}\${ext}"
    done

    module load phyling
    phyling align \\
        -I staged \\
        -m ${markerset_dir} \\
        -o align \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p align
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > align/BUSCOmarker001.fa
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > align/BUSCOmarker002.fa
    """
}
