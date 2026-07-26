//
// Shared helper functions for BFD pipelines.
//
// These are plain Nextflow functions, deliberately NOT a class under lib/.
// Groovy classes on the lib/ classpath are invisible to Nextflow's strict
// language spec — `nextflow lint` reports "`Utils` is not defined" for every
// reference — and they cannot see the implicit `params`, `log` or `workflow`
// variables, which forces those to be threaded through as extra arguments at
// every call site. Functions in an included .nf have neither problem.
//
// Import what you need:
//   include { tablesDir; clearIfStale } from '../common/utils.nf'
//

// Resolve the output subdirectory under params.tables for the MERGE steps.
//   no --taxon         → ${params.tables}/All_Taxa   (the full dataset)
//   --taxon RANK:VALUE → ${params.tables}/<sanitised VALUE>
// Used by every MERGE_* process so all merged tables land in a subfolder and
// nothing is written loose at the top level of tables/.
def tablesDir() {
    params.taxon
        ? "${params.tables}/${(params.taxon as String).split(':',2)[1].replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
        : "${params.tables}/All_Taxa"
}

// ── Taxonomy ────────────────────────────────────────────────────────────────

// Ranks recognised in samples.csv, ordered broadest → narrowest. Order is
// meaningful: query_ANI requires --query_rank to be strictly narrower than
// --compare, which is an index comparison in this list.
def taxonomicRanks() {
    ['PHYLUM','SUBPHYLUM','CLASS','SUBCLASS','ORDER','FAMILY','GENUS']
}

// Validate a rank parameter and return it upper-cased.
def assertRank(String value, String paramName) {
    def rank = (value as String).toUpperCase()
    if (!(rank in taxonomicRanks())) {
        error "--${paramName} must be one of: ${taxonomicRanks().join(', ')}"
    }
    rank
}

// Build the samples.csv row filter for --taxon RANK:VALUE. Returns an
// always-true closure when --taxon is unset, so callers can filter
// unconditionally instead of branching.
def taxonRowFilter() {
    if (!params.taxon) {
        return { _row -> true }
    }
    def parts = (params.taxon as String).split(':', 2)
    if (parts.size() != 2 || !parts[0] || !parts[1]) {
        error "--taxon must be RANK:VALUE, e.g. --taxon PHYLUM:Ascomycota"
    }
    def taxRank  = parts[0].toUpperCase()
    def taxValue = parts[1]
    log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
    return { row -> row[taxRank]?.trim() == taxValue }
}

// ── ANI ─────────────────────────────────────────────────────────────────────

// skani preset → CLI flag (medium is the tool default → no flag).
// Shared by every skani process; the sketch and the compare step MUST pass the
// same preset or the cached sketches are unusable.
def skaniPresetFlag(String p) {
    def flags = [fast: '--fast', medium: '', slow: '--slow']
    def key   = (p ?: 'medium').toLowerCase()
    if (!flags.containsKey(key)) {
        error "--skani_preset must be one of: fast, medium, slow"
    }
    flags[key]
}

// Resolve the genome filename stem for a samples.csv row under --genome_name_style.
def genomeStem(row, String nameStyle) {
    nameStyle == 'asmid'
        ? row.ASMID?.trim()
        : SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
}

// Write a group's genome-names lookup TSV — the authoritative genome universe
// for the report step, so genomes with no surviving ANI pair still appear (as
// outliers) rather than vanishing. Only genomes actually present on disk are
// listed: neither this file nor the genome list fed to the compare step may
// reference a genome missing from --genome_dir.
//
// `role` adds query_ANI's sixth column; pass null for compare_ANI's five.
def writeNamesTsv(String group_name, List metas, String prefix = 'names') {
    def withRole = metas.any { m -> m.containsKey('role') }
    def header   = "filename\tasmid\tgenus\tspecies\tstrain" + (withRole ? "\trole\n" : "\n")
    def rows     = metas.collect { m ->
        def base = "${m.id}\t${m.asmid}\t${m.genus}\t${m.species}\t${m.strain}"
        withRole ? "${base}\t${m.role}" : base
    }
    def f = file("${workflow.workDir}/${prefix}_${group_name}.tsv")
    f.text = header + rows.join('\n') + "\n"
    f
}

// ── storeDir staleness ──────────────────────────────────────────────────────

// Delete storeDir output files that are older than inputFile so Nextflow re-runs
// the process instead of using the stale cached result. Only deletes files that
// actually exist and are strictly older; missing or equal-timestamp files are left
// alone. Logs every deletion.
def clearIfStale(inputFile, List storedOutputs) {
    // Never touch real storeDir outputs during a stub run — the stub block would
    // otherwise overwrite a deleted real result with a tiny placeholder.
    if (workflow.stubRun) return
    if (!inputFile.exists()) return
    def inputMtime = inputFile.lastModified()   // follows symlinks via Java NIO default
    def stale = storedOutputs.findAll { f -> f.exists() && f.lastModified() < inputMtime }
    if (!stale.isEmpty()) {
        stale.each { f ->
            f.delete()
            log.info "Deleted stale cached output (input newer): ${f}"
        }
    }
}
