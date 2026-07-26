//
// BFD Functional Annotation Pipeline — main workflow definition.
//
// Runs Pfam, CAZy, MEROPS, SignalP, TMHMM, TargetP, IDP (AIUPred), WoLF PSORT and
// predGPI over every species in samples.csv, computes per-genome sequence and
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
include { taxonRowFilter    } from '../modules/common/utils.nf'

workflow BFD {

    // ── --help: print parameter documentation and exit ───────────────────────
    if (params.help) {
        log.info paramsHelp(command: "nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD")
        return
    }

    // ── Validate params and samples.csv structure (fail fast) ────────────────
    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    def taxonFilter = taxonRowFilter()

    // ── Base sample channel ──────────────────────────────────────────────────
    // Emits tuple(locustag, basename, species, strain) per row.
    // STRAIN: first ';'-delimited token; single quotes stripped.
    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = SampleUtils.cleanStrain(row.STRAIN?.trim() ?: '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            tuple(locustag, basename, species, strain)
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
    // Each drops genomes whose input file is absent, so a missing file skips that
    // genome for that analysis instead of failing the run.
    def proteins_ch = ready_ch
        .map { locustag, basename, species, strain ->
            def prot = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
            if (!prot.exists()) {
                log.warn "Skipping ${basename} (${locustag}): protein file not found"
                return null
            }
            return tuple(locustag, basename, species, strain, prot)
        }
        .filter { it != null }

    def aa_freq_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    def codon_freq_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.cds_dir}/${basename}.cds-transcripts.fa", glob: false)
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    def intergenic_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.gff_dir}/${basename}.gff3", glob: false)
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    def gene_stats_ch = ready_ch.map { locustag, basename, species, strain ->
        def gff = file("${params.gff_dir}/${basename}.gff3", glob: false)
        def dna = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
        (gff.exists() && dna.exists()) ? tuple(locustag, basename, gff, dna) : null
    }.filter { it != null }

    // Assembly-stats input: AAFTF assess reads the genome FASTA and the merged
    // table joins to species on ASMID, so we attach ASMID (absent from ready_ch)
    // by joining a basename→ASMID lookup built from samples.csv.
    def asmid_by_basename_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def basename = SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            tuple(basename, row.ASMID?.trim() ?: '')
        }
        .filter { basename, asmid -> asmid }

    def asm_stats_ch = ready_ch
        .map { locustag, basename, species, strain ->
            def f = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
            f.exists() ? tuple(basename, f) : null
        }
        .filter { it != null }
        .join(asmid_by_basename_ch)            // key = basename → (basename, genome, asmid)
        .map { basename, genome, asmid -> tuple(asmid, basename, genome) }

    // BUSCO lineage: params.busco_lineage (default fungi_odb12) globally. The
    // val(lineage) slot is kept explicit so per-clade overrides can be wired in
    // here later (e.g. from row.BUSCO_LINEAGE).
    def busco_genome_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
        f.exists() ? tuple(locustag, basename, params.busco_lineage, f) : null
    }.filter { it != null }

    def busco_pep_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
        f.exists() ? tuple(locustag, basename, params.busco_lineage, f) : null
    }.filter { it != null }

    // ── Run ──────────────────────────────────────────────────────────────────
    // merge_all globs the whole results dir (every genome ever processed). Under
    // --taxon that would pull in genomes outside the requested clade, so fall
    // through to the current-run (already taxon-filtered) path instead.
    def use_glob   = params.merge_all.toBoolean() && !params.taxon
    def skip_merge = params.skip_merge.toBoolean()

    BFD_FUNCTIONAL(proteins_ch, use_glob, skip_merge)

    BFD_GENOME_STATS(
        aa_freq_ch, codon_freq_ch, intergenic_ch, gene_stats_ch,
        asm_stats_ch, busco_genome_ch, busco_pep_ch
    )

    if (!skip_merge) {
        // Keys come from the post-setup channel, so a --taxon run yields a
        // manifest restricted to that clade.
        def matched_keys_file = ready_ch
            .map { locustag, basename, species, strain -> (locustag ?: '').trim() }
            .filter { it }
            .collectFile(name: 'matched_locustags.txt', newLine: true)

        BFD_MERGE(
            BFD_GENOME_STATS.out.aa_cached,
            BFD_GENOME_STATS.out.aa_csv,
            BFD_GENOME_STATS.out.codon_cached,
            BFD_GENOME_STATS.out.codon_csv,
            BFD_GENOME_STATS.out.intergenic_csv,
            BFD_GENOME_STATS.out.gene_stats,
            BFD_GENOME_STATS.out.asm_stats,
            matched_keys_file,
            use_glob
        )
    }
}
