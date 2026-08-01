include { tablesDir } from '../../common/utils.nf'

process MERGE_AA_FREQ {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "aa_freq.parquet", emit: parquet

    script:
    """
    first=1
    while IFS=\$'\\t' read -r f _mtime _size; do
        [ -n "\$f" ] || continue
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done < ${manifest} | gzip > aa_freq.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('aa_freq.csv.gz', sample_size=-1)) TO 'aa_freq.parquet' (FORMAT PARQUET);"
    rm -f aa_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > aa_freq.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('aa_freq.csv.gz', sample_size=-1)) TO 'aa_freq.parquet' (FORMAT PARQUET);"
    rm -f aa_freq.csv.gz
    """
}
