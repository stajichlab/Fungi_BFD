// MERGE_WGD_KSD — final cross-species merge of all wgd ksd ks.tsv tables into
// tables/wgd.ks.parquet (+ tables/wgd.ks.summary.parquet), following the BFD
// MERGE_* convention (publishDir tablesDir(), label 'merge').
//
// Like the gatedGlobIn idiom, the ks.tsv files are globbed off disk from
// params.outdir (never a live Nextflow channel), so a -resume run that only
// recomputes part of the dataset still merges EVERYTHING published so far.
// The input channel is therefore just a completion gate: it fires once after
// the last WGD_KSD task of this invocation (or immediately, when all ksds are
// cached from earlier waves), and the python merge does the rest.
//
// Emits:
//   wgd.ks.parquet        — one row per (genome, wgd gene pair): pruned to the
//                            analysis columns, keyed by species_prefix (LOCUSTAG)
//                            + genome (sampletag). NULLs for pairs with no valid
//                            Ka/Ks estimate (~43% of rows).
//                            Pairs with dS > params.wgd_max_ds (saturated) and, with
//                            params.wgd_drop_null_ds, pairs with no dS are dropped here;
//                            Pairs with fewer than params.wgd_min_codons codons left after
//                            wgd's family-wide gap stripping are dropped too (dS from 1-2
//                            codons is unreliable). alignmentcoverage, alignmentidentity and
//                            strippedalignmentlength are kept. The raw ks.tsv is not changed.
//   wgd.ks.summary.parquet— per-genome counts: n_pairs, n_pairs_with_ds, n_families
include { tablesDir } from '../../../common/utils.nf'

process MERGE_WGD_KSD {
    tag        'wgd_ksd_merge'
    label      'merge'

    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    val(outdir)  // absolute path of params.outdir; completion gate only
    path(sync)   // completion gate only; data is globbed off disk

    output:
    path("wgd.ks.parquet"),         emit: ks_table
    path("wgd.ks.summary.parquet"), emit: ks_summary

    script:
    def max_ds   = params.wgd_max_ds   ? "--max-ds ${params.wgd_max_ds}" : ''
    def min_cod  = params.wgd_min_codons ? "--min-codons ${params.wgd_min_codons}" : ''
    def drop_nul = params.wgd_drop_null_ds.toString().toBoolean() ? '--drop-null-ds' : ''
    """
    python3 ${projectDir}/bin/merge_wgd_ks.py \
        --glob '${outdir}/wgd_ksd/*/*.ks.tsv' \
        -o wgd.ks.parquet \
        -s wgd.ks.summary.parquet \
        ${max_ds} ${min_cod} ${drop_nul}
    """

    stub:
    """
    touch wgd.ks.parquet wgd.ks.summary.parquet
    """
}
