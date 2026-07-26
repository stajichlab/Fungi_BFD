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
