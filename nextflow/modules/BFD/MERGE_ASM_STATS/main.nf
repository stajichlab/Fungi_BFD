include { tablesDir } from '../../common/utils.nf'

process MERGE_ASM_STATS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "asm_stats.tsv.gz", emit: tsv

    script:
    """
    python3 ${params.scripts}/summarize_asm_stats.py \\
        --manifest ${manifest} \\
        --samples  ${params.samples} \\
        -o         asm_stats.tsv.gz
    """

    stub:
    """
    printf 'ASMID\\tSPECIES\\tSTRAIN\\tcontig_count\\ttotal_length_bp\\tmin_contig_bp\\tmax_contig_bp\\tmedian_contig_bp\\tmean_contig_bp\\tL50\\tN50_bp\\tL90\\tN90_bp\\tgc_pct\\tn_gap_count\\ttotal_n_bases\\tmasked_bases\\tmasked_pct\\tt2t_scaffolds\\ttelomere_fwd\\ttelomere_rev\\n' | gzip > asm_stats.tsv.gz
    """
}
