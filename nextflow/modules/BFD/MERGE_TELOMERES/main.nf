include { tablesDir } from '../../common/utils.nf'

process MERGE_TELOMERES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest
    path samples

    output:
    path "telomeres.parquet", emit: parquet

    script:
    """
    python3 ${params.scripts}/summarize_telomeres.py \
        --manifest ${manifest} \
        --samples  ${samples} \
        -o         telomeres.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomeres.tsv.gz', delim='\\t', sample_size=-1)) TO 'telomeres.parquet' (FORMAT PARQUET);"
    rm -f telomeres.tsv.gz
    """

    stub:
    """
    printf 'ASMID\tSPECIES\tSTRAIN\ttelomere_scaffolds\ttelomere_tracts\ttelomere_total_length_bp\ttelomere_total_repeats\ttelomere_5prime_count\ttelomere_3prime_count\ttelomere_plus_count\ttelomere_minus_count\ttelomere_terminal_count\ttelomere_internal_count\ttelomere_top_monomer\ttelomere_monomers\n' | gzip > telomeres.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('telomeres.tsv.gz', delim='\\t', sample_size=-1)) TO 'telomeres.parquet' (FORMAT PARQUET);"
    rm -f telomeres.tsv.gz
    """
}
