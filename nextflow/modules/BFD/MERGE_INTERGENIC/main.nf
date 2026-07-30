include { tablesDir } from '../../common/utils.nf'

process MERGE_INTERGENIC {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_intergenic_distances.csv.gz", emit: csv

    script:
    """
    first=1
    while IFS=\$'\\t' read -r f _mtime _size; do
        [ -n "\$f" ] || continue
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done < ${manifest} | gzip > gene_intergenic_distances.csv.gz
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > gene_intergenic_distances.csv.gz
    """
}
