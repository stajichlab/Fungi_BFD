include { tablesDir } from '../../common/utils.nf'

process MERGE_CODON_FREQ {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "codon_freq.parquet", emit: parquet

    script:
    """
    first=1
    while IFS=\$'\\t' read -r f _mtime _size; do
        [ -n "\$f" ] || continue
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done < ${manifest} | gzip > codon_freq.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('codon_freq.csv.gz', sample_size=-1)) TO 'codon_freq.parquet' (FORMAT PARQUET);"
    rm -f codon_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,codon,frequency\\n' | gzip > codon_freq.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('codon_freq.csv.gz', sample_size=-1)) TO 'codon_freq.parquet' (FORMAT PARQUET);"
    rm -f codon_freq.csv.gz
    """
}
