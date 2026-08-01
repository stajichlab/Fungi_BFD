//
// BFD_MERGE — consolidate per-genome statistics into the DuckDB-loadable tables,
// plus the run's samples/species manifest.
//
// Merge scope mirrors BFD_FUNCTIONAL:
//   glob mode — merge every file in the stats dir, including genomes from earlier
//               runs, gated on this run's outputs so new files are on disk first.
//               Only valid without --taxon: globbing would pull in genomes outside
//               the requested clade.
//   run mode  — merge only what this run covered.
//
// Merge inputs are passed as a manifest file (path, mtime, size per line) rather
// than staged with `path 'inputs/*'`: at ~8k genomes the staging command grew
// enormous (risking "Argument list too long") and created a symlink per file.
// Baking mtime+size into the manifest *content* also gives the staleness check
// for free — if any input is regenerated the manifest content changes, so the
// merge task's input hash changes and -resume re-runs it; if nothing changed the
// manifest is byte-identical and the merge is a cache hit.
//

include { MERGE_AA_FREQ    } from '../../modules/BFD/MERGE_AA_FREQ/main.nf'
include { MERGE_CODON_FREQ } from '../../modules/BFD/MERGE_CODON_FREQ/main.nf'
include { MERGE_INTERGENIC } from '../../modules/BFD/MERGE_INTERGENIC/main.nf'
include { MERGE_GENE_STATS } from '../../modules/BFD/MERGE_GENE_STATS/main.nf'
include { MERGE_ASM_STATS  } from '../../modules/BFD/MERGE_ASM_STATS/main.nf'
include { MERGE_SAMPLES    } from '../../modules/BFD/MERGE_SAMPLES/main.nf'

include { gatedGlobIn; toManifest } from '../../modules/common/utils.nf'

// Root a gated glob in params.genome_stats_outdir.
def gatedGlobStats(sync_ch, String glob) {
    gatedGlobIn(sync_ch, params.genome_stats_outdir, glob)
}

workflow BFD_MERGE {
    take:
    // Taken as individual channels, not as BFD_GENOME_STATS.out: a workflow's
    // .out spreads into one argument per emit, so passing it whole would arrive
    // here as 7 separate arguments.
    aa_cached       // genomes whose aa_freq output predates this run
    aa_csv          // aa_freq outputs produced by this run
    codon_cached
    codon_csv
    intergenic_csv
    gene_stats
    asm_stats
    use_glob

    main:
    // Always the full/unscoped master samples+species table now (T-014 §D.2) --
    // no matched-keys input, no per-taxon subset.
    MERGE_SAMPLES(file(params.samples))

    if (use_glob) {
        // In glob mode a table is rebuilt even when its producer did not run this
        // session, so genomes computed earlier still reach the table. Globs are
        // "**/*.ext", not "*.ext" -- these directories are hash-bucketed one level
        // deep (T-014); a flat glob would silently match nothing for any bucketed
        // genome, the same silent-loss failure class already documented for
        // gated-glob-vs-live-channel elsewhere in this repo (.living/learnings.md).
        MERGE_AA_FREQ(toManifest(gatedGlobStats(
            params.run_aa_freq.toBoolean() ? aa_csv.flatten().collect().ifEmpty([]) : Channel.of(true),
            "aa_freq/**/*.aa_freq.csv.gz"), 'aa_freq.manifest.txt'))

        MERGE_CODON_FREQ(toManifest(gatedGlobStats(
            params.run_codon_freq.toBoolean() ? codon_csv.flatten().collect().ifEmpty([]) : Channel.of(true),
            "codon_freq/**/*.codon_freq.csv.gz"), 'codon_freq.manifest.txt'))

        MERGE_INTERGENIC(toManifest(gatedGlobStats(
            params.run_intergenic.toBoolean() ? intergenic_csv.collect() : Channel.of(true),
            "intergenic_stats/**/*.gene_intergenic_distances.csv.gz"), 'intergenic.manifest.txt'))

        MERGE_GENE_STATS(toManifest(gatedGlobStats(
            params.run_gene_stats.toBoolean() ? gene_stats.collect() : Channel.of(true),
            "gene_stats/**/*.csv.gz"), 'gene_stats.manifest.txt'))

        MERGE_ASM_STATS(toManifest(gatedGlobStats(
            params.run_asm_stats.toBoolean() ? asm_stats.collect() : Channel.of(true),
            "asm_stats/**/*.stats.txt"), 'asm_stats.manifest.txt'))

    } else {
        // Run mode: the frequency merges take the union of genomes already cached
        // before the run and genomes this run produced — both resolved without any
        // existence test racing an in-flight write (plan item 13).
        if (params.run_aa_freq.toBoolean())
            MERGE_AA_FREQ(toManifest(aa_cached.mix(aa_csv.flatten()), 'aa_freq.manifest.txt'))

        if (params.run_codon_freq.toBoolean())
            MERGE_CODON_FREQ(toManifest(codon_cached.mix(codon_csv.flatten()), 'codon_freq.manifest.txt'))

        if (params.run_intergenic.toBoolean())
            MERGE_INTERGENIC(toManifest(intergenic_csv, 'intergenic.manifest.txt'))

        if (params.run_gene_stats.toBoolean())
            MERGE_GENE_STATS(toManifest(gene_stats, 'gene_stats.manifest.txt'))

        if (params.run_asm_stats.toBoolean())
            MERGE_ASM_STATS(toManifest(asm_stats, 'asm_stats.manifest.txt'))
    }
}
