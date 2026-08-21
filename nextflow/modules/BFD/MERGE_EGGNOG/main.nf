include { tablesDir } from '../../common/utils.nf'

process MERGE_EGGNOG {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(annotations)

    output:
        path("eggnog.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_eggnog.py -o eggnog.csv ${annotations}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('eggnog.csv', sample_size=-1)) TO 'eggnog.parquet' (FORMAT PARQUET);"
    rm -f eggnog.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,seed_ortholog,evalue,score,eggnog_ogs,cog_category,description,preferred_name,go_terms,ec,kegg_ko,kegg_pathway,cazy,pfams\\n' > eggnog.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('eggnog.csv', sample_size=-1)) TO 'eggnog.parquet' (FORMAT PARQUET);"
    rm -f eggnog.csv
    """
}
