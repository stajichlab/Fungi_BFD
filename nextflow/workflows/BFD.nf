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
include { taxonRowFilter; asmidRowFilter } from '../modules/common/utils.nf'
include { cleanStrain; makeSampleTag; resolveGenomeFile; dirIndex } from '../modules/common/utils.nf'

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
    def asmidFilter = asmidRowFilter()

    // ── Base sample channel ──────────────────────────────────────────────────
    // Emits one meta map per row.
    // STRAIN: first ';'-delimited token; single quotes stripped.
    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .map { row ->
            // meta.id is the filesystem-safe SPECIES_STRAIN tag and the primary key
            // every module names its outputs after (plan section 2.4).
            [
                id      : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: ''),
                locustag: row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim(),
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
    // Each drops genomes whose input file is absent, so a missing file skips that
    // genome for that analysis instead of failing the run.
    //
    // Every channel below resolves its file(s) against one of these indexes
    // instead of a per-genome .exists() call. Same fix as ANI_SAMPLES.nf's
    // genomeIndex (see dirIndex() in modules/common/utils.nf): pep_dir alone
    // was being checked per-genome by three separate channels here
    // (proteins_ch, aa_freq_ch, busco_pep_ch) for the literal same file, and
    // genome_dir was resolved twice identically (asm_stats_ch,
    // busco_genome_ch) via resolveGenomeFile(). Cheap on today's local
    // pep/cds/gff/genome dirs, but genome_dir already shares input_clean_genomes'
    // S3 convention with the ANI k8s profile, so this is latent, not
    // hypothetical -- and the redundant lookups were pure waste either way.
    def pepIndex    = dirIndex(params.pep_dir)
    def cdsIndex    = dirIndex(params.cds_dir)
    def gffIndex    = dirIndex(params.gff_dir)
    def genomeIndex = dirIndex(params.genome_dir)

    def proteins_ch = ready_ch
        .map { meta ->
            def prot = pepIndex["${meta.id}.proteins.fa"]
            if (!prot) {
                log.warn "Skipping ${meta.id} (${meta.locustag}): protein file not found"
                return null
            }
            return tuple(meta, prot)
        }
        .filter { it != null }

    def aa_freq_ch = ready_ch.map { meta ->
        def f = pepIndex["${meta.id}.proteins.fa"]
        f ? tuple(meta, f) : null
    }.filter { it != null }

    def codon_freq_ch = ready_ch.map { meta ->
        def f = cdsIndex["${meta.id}.cds-transcripts.fa"]
        f ? tuple(meta, f) : null
    }.filter { it != null }

    def intergenic_ch = ready_ch.map { meta ->
        def f = gffIndex["${meta.id}.gff3"]
        f ? tuple(meta, f) : null
    }.filter { it != null }

    def gene_stats_ch = ready_ch.map { meta ->
        def gff = gffIndex["${meta.id}.gff3"]
        def dna = genomeIndex["${meta.id}.scaffolds.fa"]
        (gff && dna) ? tuple(meta, gff, dna) : null
    }.filter { it != null }

    // Assembly-stats input: AAFTF assess reads the genome FASTA and the merged
    // table joins to species on ASMID, so we attach ASMID (absent from ready_ch)
    // by joining a basename→ASMID lookup built from samples.csv. Every genome
    // channel below needs ASMID too: resolveGenomeFile() falls back to
    // <asmid>.fa.gz / <asmid>.masked.fasta.gz when genome_dir holds raw/clean
    // assemblies (e.g. input_clean_genomes/) rather than SETUP_SYMLINKS' output.
    def asmid_by_basename_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .map { row ->
            def basename = makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            tuple(basename, row.ASMID?.trim() ?: '')
        }
        .filter { basename, asmid -> asmid }

    def ready_with_asmid_ch = ready_ch
        .map { meta -> tuple(meta.id, meta) }
        .join(asmid_by_basename_ch)            // inner join: drops rows with no ASMID
        .map { _id, meta, asmid -> meta + [asmid: asmid] }

    def asm_stats_ch = ready_with_asmid_ch
        .map { meta ->
            def f = resolveGenomeFile(genomeIndex, meta.id, meta.asmid)
            if (!f) log.warn "Skipping ${meta.id} (asm_stats): no ${meta.id}.scaffolds.fa or ${meta.asmid}.fa.gz/.masked.fasta.gz in ${params.genome_dir}"
            f ? tuple(meta, f) : null
        }
        .filter { it != null }

    // BUSCO lineage: params.busco_lineage (default fungi_odb12) globally. The
    // val(lineage) slot is kept explicit so per-clade overrides can be wired in
    // here later (e.g. from row.BUSCO_LINEAGE).
    def busco_genome_ch = ready_with_asmid_ch.map { meta ->
        def f = resolveGenomeFile(genomeIndex, meta.id, meta.asmid)
        if (!f) log.warn "Skipping ${meta.id} (busco_genome): no ${meta.id}.scaffolds.fa or ${meta.asmid}.fa.gz/.masked.fasta.gz in ${params.genome_dir}"
        f ? tuple(meta + [lineage: params.busco_lineage], f) : null
    }.filter { it != null }

    def busco_pep_ch = ready_ch.map { meta ->
        def f = pepIndex["${meta.id}.proteins.fa"]
        f ? tuple(meta + [lineage: params.busco_lineage], f) : null
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
            use_glob
        )
    }
}
