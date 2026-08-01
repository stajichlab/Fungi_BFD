include { tablesDir } from '../../common/utils.nf'

process MERGE_CAZY {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(overviews)
        path(cazymes)

    output:
        path("cazy.overview.parquet"),    emit: overview_parquet
        path("cazy.cazymes_hmm.parquet"), emit: hmm_parquet

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_cazy.py \\
        --overviews ${overviews} \\
        --cazymes   ${cazymes} \\
        --out-overview cazy.overview.csv \\
        --out-hmm      cazy.cazymes_hmm.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.overview.csv', sample_size=-1)) TO 'cazy.overview.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.cazymes_hmm.csv', sample_size=-1)) TO 'cazy.cazymes_hmm.parquet' (FORMAT PARQUET);"
    rm -f cazy.overview.csv cazy.cazymes_hmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,EC,cazyme_fam,sub_fam,diamond_fam,substrate,toolcount\\n' > cazy.overview.csv
    printf 'species_prefix,HMM_id,profile_length,protein_id,protein_length,evalue,q_start,q_end,s_start,s_end,coverage\\n' > cazy.cazymes_hmm.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.overview.csv', sample_size=-1)) TO 'cazy.overview.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.cazymes_hmm.csv', sample_size=-1)) TO 'cazy.cazymes_hmm.parquet' (FORMAT PARQUET);"
    rm -f cazy.overview.csv cazy.cazymes_hmm.csv
    """
}
