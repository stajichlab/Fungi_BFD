include { tablesDir } from '../../common/utils.nf'

process MERGE_TELOMERES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest
    path samples

    output:
    path "telomere_summary.parquet", emit: summary_parquet
    path "telomere_tracts.parquet",  emit: tracts_parquet

    script:
    """
    python3 ${params.scripts}/summarize_telomeres.py \
        --manifest    ${manifest} \
        --samples     ${samples} \
        --summary-out telomere_summary.tsv.gz \
        --tracts-out  telomere_tracts.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomere_summary.tsv.gz', delim='\\t', sample_size=-1)) TO 'telomere_summary.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomere_tracts.tsv.gz',  delim='\\t', sample_size=-1)) TO 'telomere_tracts.parquet'  (FORMAT PARQUET);"
    rm -f telomere_summary.tsv.gz telomere_tracts.tsv.gz
    """

    stub:
    """
    printf 'ASMID\tLOCUSTAG\tSPECIES\tSTRAIN\ttelomere_scaffolds\ttelomere_scaffolds_both_ends\ttelomere_tracts\ttelomere_total_length_bp\ttelomere_total_repeats\ttelomere_5prime_count\ttelomere_3prime_count\ttelomere_plus_count\ttelomere_minus_count\ttelomere_terminal_count\ttelomere_internal_count\ttelomere_top_monomer\ttelomere_monomers\n' | gzip > telomere_summary.tsv.gz
    printf 'tract_id\tASMID\tLOCUSTAG\tscaffold\tend_type\tstrand\tmonomer\trepeat_count\ttract_length\tstart\tend_coord\tterminal\tdistance_to_end\ttract_seq\tflank_seq\n' | gzip > telomere_tracts.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomere_summary.tsv.gz', delim='\\t', sample_size=-1)) TO 'telomere_summary.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomere_tracts.tsv.gz',  delim='\\t', sample_size=-1)) TO 'telomere_tracts.parquet'  (FORMAT PARQUET);"
    rm -f telomere_summary.tsv.gz telomere_tracts.tsv.gz
    """
}
