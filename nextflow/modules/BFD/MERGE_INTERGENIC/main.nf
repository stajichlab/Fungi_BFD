include { tablesDir } from '../../common/utils.nf'

process MERGE_INTERGENIC {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_intergenic_distances.parquet", emit: parquet

    script:
    """
    first=1
    while IFS=\$'\\t' read -r f _mtime _size; do
        [ -n "\$f" ] || continue
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done < ${manifest} | gzip > gene_intergenic_distances.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('gene_intergenic_distances.csv.gz', sample_size=-1)) TO 'gene_intergenic_distances.parquet' (FORMAT PARQUET);"
    rm -f gene_intergenic_distances.csv.gz
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > gene_intergenic_distances.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('gene_intergenic_distances.csv.gz', sample_size=-1)) TO 'gene_intergenic_distances.parquet' (FORMAT PARQUET);"
    rm -f gene_intergenic_distances.csv.gz
    """
}
