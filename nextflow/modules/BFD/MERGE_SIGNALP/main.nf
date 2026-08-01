include { tablesDir } from '../../common/utils.nf'

process MERGE_SIGNALP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("signalp.signal_peptide.parquet"), emit: parquet

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_signalp.py -o signalp.signal_peptide.csv ${gff3s}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('signalp.signal_peptide.csv', sample_size=-1)) TO 'signalp.signal_peptide.parquet' (FORMAT PARQUET);"
    rm -f signalp.signal_peptide.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,peptide_start,peptide_end,probability\\n' > signalp.signal_peptide.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('signalp.signal_peptide.csv', sample_size=-1)) TO 'signalp.signal_peptide.parquet' (FORMAT PARQUET);"
    rm -f signalp.signal_peptide.csv
    """
}
