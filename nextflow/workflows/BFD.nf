//
// BFD Functional Annotation Pipeline — main workflow definition.
//
// Runs the functional-annotation tools (Pfam, CAZy, MEROPS, SwissProt, SignalP,
// TMHMM, TargetP, IDP, WoLF PSORT, predGPI, eggNOG, Phobius, DeepTMHMM, TCDB,
// antiSMASH, tRNAscan-SE, Infernal, dbCAN CGC; each gated by params.run_<tool>)
// over every species in samples.csv, computes per-genome sequence and
// annotation statistics, and consolidates each into a DuckDB-loadable table.
//
// This file builds the input channels and delegates; the work lives in
// subworkflows/local/BFD_*.nf and modules/BFD/ (see REFACTOR_NEXTFLOW_PLAN.md §2.5).
//
// Key params (defined in nextflow.config / conf/profile_BFD.config):
//   --merge_all   (default true)  Merge ALL result files in the output subdirs, not
//                                 just those produced in this run. Ignored under
//                                 --taxon, which would otherwise pull in genomes
//                                 outside the requested clade.
//   --skip_merge  (default false) Skip all MERGE_* steps entirely.
//

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { INPUT_SETUP       } from '../subworkflows/local/INPUT_SETUP.nf'
include { BFD_FUNCTIONAL    } from '../subworkflows/local/BFD_FUNCTIONAL.nf'
include { BFD_GENOME_STATS  } from '../subworkflows/local/BFD_GENOME_STATS.nf'
include { BFD_MERGE         } from '../subworkflows/local/BFD_MERGE.nf'
include { taxonRowFilter; asmidRowFilter; loadSuppressSet; suppressRowFilter } from '../modules/common/utils.nf'
include { cleanStrain; makeSampleTag; resolveGenomeFile; dirIndex; lookupFile } from '../modules/common/utils.nf'

workflow BFD {

    // ── --help: print parameter documentation and exit ───────────────────────
    if (params.help) {
        log.info paramsHelp(command: "nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD")
        return
    }

    // ── Container cache must come from the environment (env-only, no hardcoded path) ─
    // params.singularity_cache is resolved from the env chain in nextflow.config;
    // this pipeline needs it for every *_sif param + the apptainer cacheDir.
    if (params.singularity_cache == null) {
        error "params.singularity_cache is not set. Export one of NXF_APPTAINER_CACHEDIR / NXF_SINGULARITY_CACHEDIR / APPTAINER_CACHEDIR / SINGULARITY_CACHEDIR pointing at the shared container cache (launch via run_unified.sh / run_funannotate.sh to get it set for you). See HOWTO_singularity.md."
    }

    // ── Validate params and samples.csv structure (fail fast) ────────────────
    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    def taxonFilter    = taxonRowFilter()
    def asmidFilter    = asmidRowFilter()
    def suppressFilter = suppressRowFilter(loadSuppressSet())

    // ── Base sample channel ──────────────────────────────────────────────────
    // Emits one meta map per row.
    // STRAIN: first ';'-delimited token; single quotes stripped.
    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .filter(suppressFilter)
        .map { row ->
            // meta.id is the filesystem-safe SPECIES_STRAIN tag and the primary key
            // every module names its outputs after (plan section 2.4).
            [
                id      : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: ''),
                locustag: row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim(),
                asmid   : row.ASMID?.trim() ?: '',
                species : row.SPECIES?.trim() ?: '',
                strain  : cleanStrain(row.STRAIN?.trim() ?: ''),
            ]
        }
        .take((params.n_test as int) > 0 ? (params.n_test as int) : -1)

    // ── Input setup: symlink predict_results files into input/ subdirs ───────
    def ready_ch
    if (params.run_setup.toBoolean()) {
        INPUT_SETUP(rows_ch)
        ready_ch = INPUT_SETUP.out.done
    } else {
        ready_ch = rows_ch
    }

    // ── Per-tool input channels ──────────────────────────────────────────────
    // Each genome's input files are looked up once, in inputs_ch. Each tool
    // channel then keeps only the genomes that have the files it needs, so a
    // missing file skips that genome for that analysis instead of failing the run.
    //
    // Lookups go through one directory listing per input dir (dirIndex() in
    // modules/common/utils.nf), not a per-genome .exists() call. The listings are
    // built lazily, on first use inside a channel closure: with run_setup=true,
    // SETUP_SYMLINKS creates the input/ links during this run, and a listing taken
    // at workflow definition would miss them (2026-10-21: 99 of ~23,500 genomes
    // ran). The first element reaches a closure only after SETUP_SYMLINKS ends.
    def indexCache = new java.util.concurrent.ConcurrentHashMap()
    def lazyIndex  = { String d -> indexCache.computeIfAbsent(d, { k -> dirIndex(k) }) }
    // Second source for the assembly-level steps: genomes that predict skipped have
    // no predict_results scaffolds, so input/dna holds no usable file for them.
    def cleanIndex = { params.clean_genome_dir ? lazyIndex.call(params.clean_genome_dir as String) : [:] }

    def bfd_debug = params.bfd_debug.toBoolean()
    def inputs_ch = ready_ch.map { meta ->
        def f = [
            pep: lookupFile(lazyIndex.call(params.pep_dir),    "${meta.id}.proteins.fa"),
            cds: lookupFile(lazyIndex.call(params.cds_dir),    "${meta.id}.cds-transcripts.fa"),
            gff: lookupFile(lazyIndex.call(params.gff_dir),    "${meta.id}.gff3"),
            dna: lookupFile(lazyIndex.call(params.genome_dir), "${meta.id}.scaffolds.fa"),
        ]
        if (!f.pep && bfd_debug) log.warn "Skipping ${meta.id} (${meta.locustag}): protein file not found"
        tuple(meta, f)
    }

    def proteins_ch   = inputs_ch.filter { _m, f -> f.pep != null }.map { m, f -> tuple(m, f.pep) }
    def codon_freq_ch = inputs_ch.filter { _m, f -> f.cds != null }.map { m, f -> tuple(m, f.cds) }
    def intergenic_ch = inputs_ch.filter { _m, f -> f.gff != null }.map { m, f -> tuple(m, f.gff) }
    // (meta, gff3, genome): gene stats, and antiSMASH / tRNAscan-SE / Infernal.
    def gene_stats_ch = inputs_ch.filter { _m, f -> f.gff != null && f.dna != null }.map { m, f -> tuple(m, f.gff, f.dna) }
    // (meta, gff3, proteins): dbCAN CGC/PUL.
    def cgc_ch        = inputs_ch.filter { _m, f -> f.gff != null && f.pep != null }.map { m, f -> tuple(m, f.gff, f.pep) }
    def busco_pep_ch  = proteins_ch.map { m, p -> tuple(m + [lineage: params.busco_lineage], p) }

    // Assembly-level input (asm_stats, telomeres, BUSCO genome). These steps name
    // their outputs by ASMID, so rows without one are skipped. resolveGenomeFile()
    // tries <id>.scaffolds.fa, then <asmid>.fa.gz / <asmid>.masked.fasta.gz, in
    // genome_dir; failing that, the same names in clean_genome_dir, and the genome
    // is marked annotated=false.
    def assembly_ch = inputs_ch
        .filter { m, _f -> m.asmid }
        .map { m, _f ->
            def g = resolveGenomeFile(lazyIndex.call(params.genome_dir), m.id, m.asmid)
            def annotated = g != null
            g = g ?: resolveGenomeFile(cleanIndex.call(), m.id, m.asmid)
            if (!g) { log.warn "Skipping ${m.id} (asm_stats/telomeres/BUSCO genome): no ${m.id}.scaffolds.fa or ${m.asmid}.fa.gz/.masked.fasta.gz in ${params.genome_dir} or ${params.clean_genome_dir}" }
            g ? tuple(m + [annotated: annotated], g) : null
        }
        .filter { it != null }
    // BUSCO lineage is params.busco_lineage (default fungi_odb12) for every genome.
    def busco_genome_ch = assembly_ch.map { m, g -> tuple(m + [lineage: params.busco_lineage], g) }

    // ── Run ──────────────────────────────────────────────────────────────────
    // merge_all globs the whole results dir (every genome ever processed). Under
    // --taxon that would pull in genomes outside the requested clade, so fall
    // through to the current-run (already taxon-filtered) path instead.
    def use_glob   = params.merge_all.toBoolean() && !params.taxon
    def skip_merge = params.skip_merge.toBoolean()

    BFD_FUNCTIONAL(proteins_ch, gene_stats_ch, cgc_ch, use_glob, skip_merge)

    BFD_GENOME_STATS(
        proteins_ch, codon_freq_ch, intergenic_ch, gene_stats_ch,
        assembly_ch, assembly_ch, busco_genome_ch, busco_pep_ch
    )

    if (!skip_merge) {
        // MERGE_SAMPLES always builds the full/unscoped master table now (T-014
        // §D.2), so no matched-keys manifest is needed here anymore.
        BFD_MERGE(
            BFD_GENOME_STATS.out.aa_cached,
            BFD_GENOME_STATS.out.aa_csv,
            BFD_GENOME_STATS.out.codon_cached,
            BFD_GENOME_STATS.out.codon_csv,
            BFD_GENOME_STATS.out.intergenic_csv,
            BFD_GENOME_STATS.out.gene_stats,
            BFD_GENOME_STATS.out.asm_stats,
            BFD_GENOME_STATS.out.telomeres,
            BFD_GENOME_STATS.out.busco_genome,
            BFD_FUNCTIONAL.out.sync,
            use_glob
        )
    }
}
