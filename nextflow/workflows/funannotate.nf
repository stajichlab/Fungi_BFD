//
// funannotate — genome cleaning, RNA-seq acquisition, gene prediction and
// functional annotation.
//
// Process definitions live in modules/funannotate/{setup,genome,rnaseq,predict,function}/.
// Helper functions live in modules/funannotate/utils.nf.
//
// Usage (from project root):
//   nextflow run nextflow/main.nf --pipeline funannotate \
//       -c nextflow/nextflow.config -profile funannotate -resume
//

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { makeSampleTag; loadSuppressSet; suppressRowFilter; loadMicrosporidiaOutSet } from '../modules/common/utils.nf'
include { loadAbinitioReuseMap; loadRnaseqRepresentativeOverride; loadHybridParentage } from '../modules/funannotate/utils.nf'

// ── Subworkflows ────────────────────────────────────────────────────────────
include { FUNANNOTATE_GENOME_PREP }       from '../subworkflows/local/FUNANNOTATE_GENOME_PREP.nf'
include { FUNANNOTATE_RNASEQ }            from '../subworkflows/local/FUNANNOTATE_RNASEQ.nf'
include { FUNANNOTATE_PREDICTION }        from '../subworkflows/local/FUNANNOTATE_PREDICTION.nf'
include { FUNANNOTATE_ANNOTATION }        from '../subworkflows/local/FUNANNOTATE_ANNOTATION.nf'


workflow FUNANNOTATE {
    // ── --help: print parameter documentation and exit ──────────────────────────
    if (params.help) {
        log.info paramsHelp(command: "nextflow run nextflow/main.nf --pipeline funannotate -c nextflow/nextflow.config -profile funannotate")
        return
    }
    // ── Container cache must come from the environment (env-only, no hardcoded path) ─
    // params.singularity_cache is resolved from the env chain in nextflow.config;
    // this pipeline needs it for every *_sif param + the apptainer cacheDir.
    if (params.singularity_cache == null) {
        throw new Exception("params.singularity_cache is not set. Export one of NXF_APPTAINER_CACHEDIR / NXF_SINGULARITY_CACHEDIR / APPTAINER_CACHEDIR / SINGULARITY_CACHEDIR pointing at the shared container cache (launch via run_funannotate.sh to get it set for you). See HOWTO_singularity.md.")
    }
    // ── Validate params and samples.csv structure (fail fast) ───────────────────
    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    // target/training_target get interpolated directly into task bash scripts
    // (FUNANNOTATE_PREDICT, FUNANNOTATE_TRAIN, BACKFILL_ABINITIO_PARAMS, ...),
    // where a relative path resolves against that task's own isolated work
    // directory, not launchDir -- silently writing real prediction output
    // somewhere that `cleanup = true` deletes after the run. A -params-file
    // (unlike this file's own defaults, which interpolate ${launchDir}) can
    // only set a literal string, so normalize to absolute here regardless of
    // what any params file provided.
    params.target          = file(params.target as String).toAbsolutePath().toString()
    params.training_target = file(params.training_target as String).toAbsolutePath().toString()

    def suppressSet    = loadSuppressSet()
    def suppressFilter = suppressRowFilter(suppressSet)

    def forceIndependentSet = (params.force_independent_species as String)
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .toSet()
    if (forceIndependentSet) {
        log.info "force_independent_species: ${forceIndependentSet.join(', ')}"
    }

    // Finer-grained than forceIndependentSet (species-keyed): forces GeneMark
    // specifically to run independently for these OUT values, while leaving
    // AUGUSTUS/SNAP reuse (forceIndependentSet) untouched for the same strain.
    // See nextflow/docs/GENEMARK_RUN_DESIGN.md.
    def forceIndependentGenemarkSet = (params.genemark_force_independent_strains as String)
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .toSet()
    if (forceIndependentGenemarkSet) {
        log.info "genemark_force_independent_strains: ${forceIndependentGenemarkSet.join(', ')}"
    }

    // Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    // out values that get a Prodigal track alongside GeneMark-ES when also
    // below predict_min_asm_bp.
    def microsporidiaSet = loadMicrosporidiaOutSet()

    // ── Taxonomy filter ───────────────────────────────────────────────────────
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { row -> true }
    }

    // ── ASMID filter ──────────────────────────────────────────────────────────
    def asmidFilter = params.asmid
        ? { row -> row.ASMID?.trim() == (params.asmid as String).trim() }
        : { row -> true }
    if (params.asmid) {
        log.info "ASMID filter: processing only '${params.asmid}'"
    }

    // ── Build jobs channel ────────────────────────────────────────────────────
    def jobs = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .filter(suppressFilter)
        .map { row ->
            def species       = (row.SPECIES?.trim() ?: '').replaceAll(/['"]/, '')
            def strain        = (row.STRAIN?.trim() ?: '').replaceAll(/['"]/, '').replaceAll(/;.*$/, '').trim().replace(':', ' ')
            def out           = makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def asmid         = row.ASMID?.trim()
            def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco         = row.BUSCO_LINEAGE?.trim()
            def header_length = 24
            def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
            def taxonid       = row.NCBI_TAXONID?.trim()
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid ->
            def gz = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, gz, _tid ->
            if (!gz.exists()) {
                log.warn "Missing genome for ${out} (asmid=${asmid}): ${gz}"
                return false
            }
            if (params.debug.toBoolean()) {
                log.info "Queuing ${out}: genome=${gz} (${gz.size()} bytes)"
            }
            return true
        }

    if (params.debug.toBoolean()) {
        jobs.view { t -> "[CHANNEL] Submitting: out=${t[0]}, asmid=${t[1]}, transl_table=${t[7]}, gz=${t[8]}" }
    }

    // ── Stages ────────────────────────────────────────────────────────────────
    FUNANNOTATE_GENOME_PREP(jobs)

    if (params.stop_after_genome_prep.toBoolean()) {
        log.info "Stopping after genome prep (stop_after_genome_prep) -- " +
                  "skipping RNA-seq acquisition, prediction, and annotation entirely."
        return
    }

    // Loaded before FUNANNOTATE_RNASEQ so its species-level "which assembly
    // anchors the shared Trinity-GG build" pick can use the same ANI/BUSCO
    // representative as FUNANNOTATE_PREDICTION's ab-initio reuse decision,
    // instead of an unrelated, arbitrary pick (see .github issue #3).
    def abinitioReuseMap = loadAbinitioReuseMap()

    if (params.run_ani_reuse.toBoolean() && abinitioReuseMap.isEmpty()) {
        error "run_ani_reuse=true but no ab-initio reuse assignments found at " +
              "${params.abinitio_reuse_csv}. Run compare_ANI (--pipeline compare_ani, " +
              "--compare SPECIES --run_ani_reuse true) first, or set --run_ani_reuse false " +
              "to train all strains independently."
    }

    // species_tag -> out overrides for which strain RNASEQ_PREPARE builds the shared
    // Trinity-GG assembly against, when abinitioReuseMap's ANI+BUSCO pick assembles
    // fine but its RNA-seq doesn't align (see loadRnaseqRepresentativeOverride()).
    // Populated by scripts/pick_rnaseq_representative_override.py, not by this pipeline.
    def rnaseqRepOverride = loadRnaseqRepresentativeOverride()

    // hybrid_species_tag -> [parent_species_tag, ...] for interspecific hybrid-cross
    // species -- routes them to composite parent-transcript RNA-seq evidence instead
    // of the normal per-species representative pick (see
    // nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md). Empty map (no hybrid_parentage.csv
    // yet) means every species behaves exactly as before this feature existed.
    def hybridParentage = loadHybridParentage()

    FUNANNOTATE_RNASEQ(FUNANNOTATE_GENOME_PREP.out.predict_genome, abinitioReuseMap, rnaseqRepOverride, hybridParentage)

    // FUNANNOTATE_RNASEQ only assigns its predict_input emission past the
    // stop_after_sra_fetch/stop_after_sra_query gates (see its "end if" markers);
    // with run_sra_fetch on and either gate set, that channel is never populated,
    // so calling PREDICTION/ANNOTATION here would hand them an undefined channel
    // instead of actually stopping the pipeline early as those params intend.
    if (params.run_sra_fetch.toBoolean() &&
        (params.stop_after_sra_fetch.toBoolean() || params.stop_after_sra_query.toBoolean())) {
        log.info "Stopping after RNA-seq acquisition (stop_after_sra_fetch/stop_after_sra_query) " +
                  "-- skipping prediction and annotation."
        return
    }

    FUNANNOTATE_PREDICTION(FUNANNOTATE_RNASEQ.out.predict_input, abinitioReuseMap, forceIndependentSet, forceIndependentGenemarkSet, microsporidiaSet)

    FUNANNOTATE_ANNOTATION(FUNANNOTATE_PREDICTION.out.metadata, FUNANNOTATE_RNASEQ.out.reads, taxonFilter, asmidFilter, suppressFilter)
}
