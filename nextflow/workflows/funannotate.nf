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

include { makeSampleTag }       from '../modules/common/utils.nf'
include { loadAbinitioReuseMap } from '../modules/funannotate/utils.nf'

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
    // ── Validate params and samples.csv structure (fail fast) ───────────────────
    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    def suppressSet = (params.suppress && file(params.suppress).exists())
        ? file(params.suppress).readLines()
              .collect { it.trim().split(',')[0].trim() }
              .findAll { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }

    def forceIndependentSet = (params.force_independent_species as String)
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .toSet()
    if (forceIndependentSet) {
        log.info "force_independent_species: ${forceIndependentSet.join(', ')}"
    }

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
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid ->
            if (suppressSet.contains(asmid)) {
                log.info "Suppressing ${out} (asmid=${asmid})"
                return false
            }
            return true
        }
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

    FUNANNOTATE_RNASEQ(FUNANNOTATE_GENOME_PREP.out.predict_genome)

    def abinitioReuseMap = loadAbinitioReuseMap()

    if (params.run_ani_reuse.toBoolean() && abinitioReuseMap.isEmpty()) {
        error "run_ani_reuse=true but no ab-initio reuse assignments found at " +
              "${params.abinitio_reuse_csv}. Run compare_ANI (--pipeline compare_ani, " +
              "--compare SPECIES --run_ani_reuse true) first, or set --run_ani_reuse false " +
              "to train all strains independently."
    }

    FUNANNOTATE_PREDICTION(FUNANNOTATE_RNASEQ.out.predict_input, abinitioReuseMap, forceIndependentSet)

    FUNANNOTATE_ANNOTATION(FUNANNOTATE_PREDICTION.out.metadata, FUNANNOTATE_RNASEQ.out.reads, taxonFilter, asmidFilter)
}
