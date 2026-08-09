include { tablesDir } from '../../common/utils.nf'

process MERGE_ASM_STATS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest
    path samples

    output:
    path "asm_stats.parquet", emit: parquet

    script:
    """
    python3 ${projectDir}/bin/summarize_asm_stats.py \\
        --manifest ${manifest} \\
        --samples  ${samples} \\
        -o         asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    """

    stub:
    """
    printf 'ASMID\\tSPECIES\\tSTRAIN\\tcontig_count\\ttotal_length_bp\\tmin_contig_bp\\tmax_contig_bp\\tmedian_contig_bp\\tmean_contig_bp\\tL50\\tN50_bp\\tL90\\tN90_bp\\tgc_pct\\tn_gap_count\\ttotal_n_bases\\tmasked_bases\\tmasked_pct\\tt2t_scaffolds\\ttelomere_fwd\\ttelomere_rev\\n' | gzip > asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    """
}
