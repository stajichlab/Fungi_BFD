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
include { MERGE_BUSCO_GENOME } from '../../modules/BFD/MERGE_BUSCO_GENOME/main.nf'
include { MERGE_TELOMERES  } from '../../modules/BFD/MERGE_TELOMERES/main.nf'
include { MERGE_SAMPLES    } from '../../modules/BFD/MERGE_SAMPLES/main.nf'
include { BUILD_DUCKDB     } from '../../modules/BFD/BUILD_DUCKDB/main.nf'

include { gatedGlobIn; toManifest } from '../../modules/common/utils.nf'

// Merge input for one stats table, as a manifest (see toManifest()).
//   glob mode: every file under params.genome_stats_outdir matching `glob`,
//              listed once `live` (this run's outputs) has finished, or at once
//              when the producer did not run this session. Globs are "**/*.ext"
//              because the directories are hash-bucketed one level deep (T-014);
//              a flat glob would match nothing.
//   run mode:  the files in `runFiles`; null when the producer did not run, so
//              the caller skips that merge.
// `live` emits nothing when every genome was already cached and the producer
// had no tasks (BATCH_*_FREQ); the glob then never fires and that merge is
// skipped, the same as before this helper.
def statsManifest(boolean useGlob, boolean ran, live, runFiles, String glob, String name) {
    if (useGlob) {
        def sync = ran ? live.collect() : Channel.of(true)
        return toManifest(gatedGlobIn(sync, params.genome_stats_outdir, glob), name)
    }
    return ran ? toManifest(runFiles, name) : null
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
    telomeres
    busco_genome
    functional_sync // BFD_FUNCTIONAL.out.sync -- gates BUILD_DUCKDB on the
                     // pfam/cazy/merops/... merges too, not just genome_stats
                     // (see comment above merge_outs in BFD_FUNCTIONAL.nf)
    use_glob

    main:
    // Always the full/unscoped master samples+species table now (T-014 §D.2) --
    // no matched-keys input, no per-taxon subset.
    MERGE_SAMPLES(file(params.samples))

    // Gating list for BUILD_DUCKDB: collects one output channel from every
    // MERGE_* process actually invoked below (any single emit suffices -- all
    // of a process's emits complete together), so the duckdb rebuild (if
    // enabled) waits for all of this run's merges to publish before it globs
    // tables/ -- same purpose as gatedGlobIn's sync_ch, just fanned in from
    // process outputs instead of an upstream file-existence channel.
    // functional_sync carries the same signal for BFD_FUNCTIONAL's merges
    // (pfam/cazy/merops/swissprot/signalp/tmhmm/targetp/idp/wolfpsort/predgpi)
    // -- without it BUILD_DUCKDB's cache key never depends on those tables.
    def merge_outs = [MERGE_SAMPLES.out.samples, MERGE_SAMPLES.out.species, functional_sync]

    // In glob mode a table is rebuilt even when its producer did not run this
    // session, so genomes computed earlier still reach the table. The frequency
    // tables in run mode take the union of genomes cached before the run and
    // genomes this run computed (BFD_GENOME_STATS.planFreq), so neither list
    // races an in-flight write.
    def run_aa    = params.run_aa_freq.toBoolean()
    def run_codon = params.run_codon_freq.toBoolean()

    def aa_in = statsManifest(use_glob, run_aa, aa_csv, aa_cached.mix(aa_csv.flatten()),
        "aa_freq/**/*.aa_freq.csv.gz", 'aa_freq.manifest.txt')
    if (aa_in) {
        MERGE_AA_FREQ(aa_in)
        merge_outs << MERGE_AA_FREQ.out.parquet
    }

    def codon_in = statsManifest(use_glob, run_codon, codon_csv, codon_cached.mix(codon_csv.flatten()),
        "codon_freq/**/*.codon_freq.csv.gz", 'codon_freq.manifest.txt')
    if (codon_in) {
        MERGE_CODON_FREQ(codon_in)
        merge_outs << MERGE_CODON_FREQ.out.parquet
    }

    def intergenic_in = statsManifest(use_glob, params.run_intergenic.toBoolean(), intergenic_csv, intergenic_csv,
        "intergenic_stats/**/*.gene_intergenic_distances.csv.gz", 'intergenic.manifest.txt')
    if (intergenic_in) {
        MERGE_INTERGENIC(intergenic_in)
        merge_outs << MERGE_INTERGENIC.out.parquet
    }

    def gene_stats_in = statsManifest(use_glob, params.run_gene_stats.toBoolean(), gene_stats, gene_stats,
        "gene_stats/**/*.csv.gz", 'gene_stats.manifest.txt')
    if (gene_stats_in) {
        MERGE_GENE_STATS(gene_stats_in)
        merge_outs << MERGE_GENE_STATS.out.gene_info
    }

    def asm_in = statsManifest(use_glob, params.run_asm_stats.toBoolean(), asm_stats, asm_stats,
        "asm_stats/**/*.stats.txt", 'asm_stats.manifest.txt')
    if (asm_in) {
        MERGE_ASM_STATS(asm_in, file(params.samples))
        merge_outs << MERGE_ASM_STATS.out.parquet
    }

    def telomeres_in = statsManifest(use_glob, params.run_telomeres.toBoolean(), telomeres, telomeres,
        "telomeres/**/*.telomeres.tsv.gz", 'telomeres.manifest.txt')
    if (telomeres_in) {
        MERGE_TELOMERES(telomeres_in, file(params.samples))
        merge_outs << MERGE_TELOMERES.out.summary_parquet
    }

    def busco_in = statsManifest(use_glob, params.run_busco_genome.toBoolean(), busco_genome, busco_genome,
        "BUSCO_genome/**/*.BUSCO_summary.*.txt", 'busco_genome.manifest.txt')
    if (busco_in) {
        MERGE_BUSCO_GENOME(busco_in)
        merge_outs << MERGE_BUSCO_GENOME.out.parquet
    }

    if (params.run_build_duckdb.toBoolean()) {
        // Nextflow's strict parser rejects the Groovy spread operator
        // (`mix(*merge_outs)`), hence the explicit fold.
        def sync_ch = merge_outs.inject(Channel.empty()) { acc, ch -> acc.mix(ch) }
        BUILD_DUCKDB(sync_ch.collect())
    }
}
