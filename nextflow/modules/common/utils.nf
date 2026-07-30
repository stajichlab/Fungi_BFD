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

// ── Sample identity ─────────────────────────────────────────────────────────
//
// Moved here from lib/SampleUtils.groovy: Groovy classes on the lib/ classpath
// are invisible to Nextflow's strict language spec (§2.1). The Python equivalents
// are collect_busco_stats.py::build_basename_map() and
// fix_low_trinity.py::species_key_from_row() — keep them in sync.

// Canonicalise a raw strain value into a human-readable strain name.
// 
// - strips leading/trailing whitespace and quote characters (', ")
// - takes only the first semicolon-delimited token (some rows list
// multiple synonymous strains separated by ';')
// - replaces colons with spaces (colons appear as ' colon ' separators)
// - handles asterisks: a '*' at the start or end of the strain is removed
// entirely; a '*' between two other words is replaced with '-'
// 
// Examples:
// cleanStrain("Af293; CBS 101") → "Af293"
// cleanStrain("T-34*")          → "T-34"
// cleanStrain("ABC*DEF")        → "ABC-DEF"
// cleanStrain("ARSEF * 2860")   → "ARSEF-2860"
def cleanStrain(String rawStrain) {
    return (rawStrain ?: '').trim()
                .replaceAll(/['"]/, '')
                .split(';')[0]
                .trim()
                .replace(':', ' ')
                .replaceAll(/^\s*\*+/, '')      // '*' at the start of the strain -> removed
                .replaceAll(/\*+\s*$/, '')      // '*' at the end of the strain   -> removed
                .replaceAll(/\s*\*+\s*/, '-')   // '*' between two words          -> '-'
                .trim()
}

// Canonicalise a raw strain value into a human-readable strain name.
// 
// - strips leading/trailing whitespace and quote characters (', ")
// - takes only the first semicolon-delimited token (some rows list
// multiple synonymous strains separated by ';')
// - replaces colons with spaces (colons appear as ' colon ' separators)
// - handles asterisks: a '*' at the start or end of the strain is removed
// entirely; a '*' between two other words is replaced with '-'
// 
// Examples:
// cleanStrain("Af293; CBS 101") → "Af293"
// cleanStrain("T-34*")          → "T-34"
// cleanStrain("ABC*DEF")        → "ABC-DEF"
// cleanStrain("ARSEF * 2860")   → "ARSEF-2860"
// static String cleanStrain(String rawStrain) {
// return (rawStrain ?: '').trim()
// .replaceAll(/['"]/, '')
// .split(';')[0]
// .trim()
// .replace(':', ' ')
// .replaceAll(/^\s*\*+/, '')      // '*' at the start of the strain -> removed
// .replaceAll(/\*+\s*$/, '')      // '*' at the end of the strain   -> removed
// .replaceAll(/\s*\*+\s*/, '-')   // '*' between two words          -> '-'
// .trim()
// }
// Build a filesystem-safe "{species}_{strain}" tag from raw samples.csv values.
// 
// Canonicalises:
// - strips leading/trailing whitespace and quote characters (', ") from both fields
// - cleans the strain via {@link #cleanStrain} (first ';' token, colon→space,
// asterisk handling)
// - collapses runs of whitespace, /, #, [, ], ?, {, } into single underscores
// 
// Examples:
// makeSampleTag("Saccharomyces cerevisiae", "CBS 1171")    → "Saccharomyces_cerevisiae_CBS_1171"
// makeSampleTag("Aspergillus fumigatus", "Af293; CBS 101")  → "Aspergillus_fumigatus_Af293"
// makeSampleTag("Fusarium oxysporum", "")                   → "Fusarium_oxysporum"
// makeSampleTag("Beauveria bassiana", "ARSEF 2860*")        → "Beauveria_bassiana_ARSEF_2860"
def makeSampleTag(String rawSpecies, String rawStrain) {
    def sp = (rawSpecies ?: '').trim().replaceAll(/['"]/, '')
    def st = cleanStrain(rawStrain)
    return [sp, st].findAll { it }
                   .join('_')
                   .replaceAll(/[\s\/\#\[\]\?\{\}]+/, '_')
}

// Canonicalise a raw strain value into a human-readable strain name.
// 
// - strips leading/trailing whitespace and quote characters (', ")
// - takes only the first semicolon-delimited token (some rows list
// multiple synonymous strains separated by ';')
// - replaces colons with spaces (colons appear as ' colon ' separators)
// - handles asterisks: a '*' at the start or end of the strain is removed
// entirely; a '*' between two other words is replaced with '-'
// 
// Examples:
// cleanStrain("Af293; CBS 101") → "Af293"
// cleanStrain("T-34*")          → "T-34"
// cleanStrain("ABC*DEF")        → "ABC-DEF"
// cleanStrain("ARSEF * 2860")   → "ARSEF-2860"
// static String cleanStrain(String rawStrain) {
// return (rawStrain ?: '').trim()
// .replaceAll(/['"]/, '')
// .split(';')[0]
// .trim()
// .replace(':', ' ')
// .replaceAll(/^\s*\*+/, '')      // '*' at the start of the strain -> removed
// .replaceAll(/\*+\s*$/, '')      // '*' at the end of the strain   -> removed
// .replaceAll(/\s*\*+\s*/, '-')   // '*' between two words          -> '-'
// .trim()
// }
// Build a filesystem-safe "{species}_{strain}" tag from raw samples.csv values.
// 
// Canonicalises:
// - strips leading/trailing whitespace and quote characters (', ") from both fields
// - cleans the strain via {@link #cleanStrain} (first ';' token, colon→space,
// asterisk handling)
// - collapses runs of whitespace, /, #, [, ], ?, {, } into single underscores
// 
// Examples:
// makeSampleTag("Saccharomyces cerevisiae", "CBS 1171")    → "Saccharomyces_cerevisiae_CBS_1171"
// makeSampleTag("Aspergillus fumigatus", "Af293; CBS 101")  → "Aspergillus_fumigatus_Af293"
// makeSampleTag("Fusarium oxysporum", "")                   → "Fusarium_oxysporum"
// makeSampleTag("Beauveria bassiana", "ARSEF 2860*")        → "Beauveria_bassiana_ARSEF_2860"
// static String makeSampleTag(String rawSpecies, String rawStrain) {
// def sp = (rawSpecies ?: '').trim().replaceAll(/['"]/, '')
// def st = cleanStrain(rawStrain)
// return [sp, st].findAll { it }
// .join('_')
// .replaceAll(/[\s\/\#\[\]\?\{\}]+/, '_')
// }
// Sanitise a single free-text taxonomic value (e.g. a samples.csv SPECIES,
// GENUS, or FAMILY field) into a filesystem-safe tag, for use as a
// `group_name` in file/directory names (e.g. `${group_name}.full.ani.tsv`
// in compare_ANI.nf/query_ANI.nf). GENUS/FAMILY/etc. rarely need this (no
// spaces), but SPECIES values ("Aspergillus fumigatus") do -- a raw,
// unsanitised group_name with an embedded space is a real risk for
// downstream shell globbing/quoting in the same processes that consume it.
// 
// Uses the same whitespace/quote/special-char handling as {@link
// #makeSampleTag}'s final step, so a sanitised SPECIES value here matches
// the species-only prefix of makeSampleTag's own output for the same
// species string.
// 
// Examples:
// sanitizeTag("Aspergillus fumigatus") → "Aspergillus_fumigatus"
// sanitizeTag("Fusarium/oxysporum complex") → "Fusarium_oxysporum_complex"
def sanitizeTag(String raw) {
    return (raw ?: '').trim()
                .replaceAll(/['"]/, '')
                .replaceAll(/[\s\/\#\[\]\?\{\}]+/, '_')
}

// ── Taxonomy ────────────────────────────────────────────────────────────────

// Ranks recognised in samples.csv, ordered broadest → narrowest. Order is
// meaningful: query_ANI requires --query_rank to be strictly narrower than
// --compare, which is an index comparison in this list.
def taxonomicRanks() {
    ['PHYLUM','SUBPHYLUM','CLASS','SUBCLASS','ORDER','FAMILY','GENUS','SPECIES']
}

// Validate a rank parameter and return it upper-cased.
def assertRank(String value, String paramName) {
    def rank = (value as String).toUpperCase()
    if (!(rank in taxonomicRanks())) {
        error "--${paramName} must be one of: ${taxonomicRanks().join(', ')}"
    }
    rank
}

// Given a meta map (may contain empty strings for undefined ranks) and a
// query_rank (e.g. GENUS), walk up the hierarchy from the rank ABOVE
// query_rank to find the first non-empty rank. Returns the sanitised fallback
// rank name, or null if no defined rank is found above query_rank.
//
// Example: query_rank=GENUS, meta.genus='', meta.family='Lasiosphaeriaceae'
//          → returns 'Lasiosphaeriaceae'
// Example: query_rank=GENUS, meta.genus='', meta.family='', meta.order='Sordariales'
//          → returns 'Sordariales'
def fallbackRank(Map meta, String queryRank) {
    def ranks  = taxonomicRanks()
    def qi     = ranks.indexOf(queryRank)
    if (qi < 0) return null
    def above  = ranks[0..<qi]
    def found  = above.findResult { r ->
        def key = r.toLowerCase()
        def val = r == 'CLASS' ? (meta['class'] ?: '') : (meta[key] ?: '')
        val.trim() ? sanitizeTag(val) : null
    }
    found
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

// Build the samples.csv row filter for --asmid. Returns an always-true
// closure when --asmid is unset, so callers can filter unconditionally
// instead of branching. Mirrors taxonRowFilter() above; previously only
// funannotate.nf and earlgrey_mask.nf applied this filter, so an
// --asmid-restricted funannotate run still computed ANI/BUSCO
// representative-picking (compare_ANI.nf, BFD.nf genome_stats) against the
// full unfiltered sample set.
def asmidRowFilter() {
    if (!params.asmid) {
        return { _row -> true }
    }
    def asmidValue = (params.asmid as String).trim()
    log.info "ASMID filter: processing only '${asmidValue}'"
    return { row -> row.ASMID?.trim() == asmidValue }
}

// ── ANI ─────────────────────────────────────────────────────────────────────

// skani preset → CLI flag (medium is the tool default → no flag).
// Shared by every skani process; the sketch and the compare step MUST pass the
// same preset or the cached sketches are unusable.
def skaniCpusFor(int n) {
    n > 500 ? 32 : n > 200 ? 16 : 8
}

def skaniMemoryFor(int n, int attempt = 1) {
    def baseGB  = n > 500 ? 48 : n > 200 ? 32 : 16
    def floorGB = attempt >= 3 ? 100 : attempt == 2 ? 24 : 6
    Math.max(baseGB, floorGB).toString() + ' GB'
}

// Wall-clock ceiling for the ANI compare processes (skani/fastani/mash/sourmash),
// scaled by group size and bumped on retry — same shape as skaniMemoryFor.
// Mainly matters on k8s: `time` maps to the pod's activeDeadlineSeconds there,
// and (unlike the HPC profile, which sets a static per-label `time` via config
// that overrides this) nothing else caps it. Floors at 1h even for the
// smallest group on the first attempt.
//
// Capped at 6h: verified live (2026-07-30) that the ucr-stajichlab namespace's
// opportunistic priority class gets activeDeadlineSeconds=21600 (6h) injected
// on EVERY pod, task pods and the long-lived head pod alike, regardless of
// what's requested — this is a cluster-side policy, not something a namespace
// user can raise. Requesting more than 6h here would just be a no-op (or
// worse, ambiguous whether the cluster clamps or rejects it), so don't.
// The real consequence is on the *head* pod (k8/head-pod.yaml), which has no
// `time`/`activeDeadlineSeconds` of its own and still gets killed at 6h --
// any run (or batch of runs via ani-suite.sh) spanning longer than that loses
// its in-flight nextflow head process(es) and, since the launch dir isn't on
// the PVC (see k8/README_compare_ani.md), its -resume cache too. Recreate the
// head pod (`kubectl apply -f k8/head-pod.yaml`) and re-launch with -resume;
// already-published per-group *.ani.tsv on S3 are safe, but expect already-
// finished groups to be recomputed since the resume cache is gone.
def aniTimeFor(int n, int attempt = 1) {
    def baseH  = n > 500 ? 6 : n > 200 ? 3 : 1
    def floorH = attempt >= 2 ? 6 : 1
    Math.min(6, Math.max(baseH, floorH))
        .toString() + ' h'
}

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
        : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
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
    def header   = "filename\tasmid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain" + (withRole ? "\trole" : "") + "\n"
    def rows     = metas.collect { m ->
        def base = "${m.id}\t${m.asmid}\t${m.phylum}\t${m['class'] ?: ''}\t${m.order}\t${m.family}\t${m.genus}\t${m.species}\t${m.strain}"
        withRole ? "${base}\t${m.role}" : base
    }
    def f = file("${workflow.workDir}/${prefix}_${group_name}.tsv")
    f.text = header + rows.join('\n') + "\n"
    f
}

// Clean genomes are stored gzip-compressed (.fa.gz) to save space (see
// funannotate.nf GENOME_CLEAN and earlgrey_mask.nf EARLGREY_BUILD_LIB).
// Given the uncompressed base path, return the existing file, preferring the
// compressed form, else the plain path (so .exists() reports missing when
// neither is present).
def genomeFile(String base) {
    def gz = file("${base}.gz", glob: false)
    if (gz.exists() && gz.size() > 0) return gz
    return file(base, glob: false)
}

// Build a filename -> Path index for a directory with ONE listing, for
// callers that would otherwise run a .exists() check per item in a
// samples.csv-scale loop (each one a real network round trip on an S3-backed
// dir — see ANI_SAMPLES.nf's genomeIndex, same fix, verified live 2026-07-30
// against a phylum-scale run: 20-30+ min of serial S3 latency collapsed into
// one listing). BFD.nf's genome/pep/cds/gff dirs are local paths today, so
// this is currently just removing redundant stat() calls -- but genome_dir
// already shares input_clean_genomes' S3 convention with the ANI k8s
// profile, so this is latent, not hypothetical.
def dirIndex(String dir) {
    def dirPath = file(dir)
    if (!dirPath.exists()) {
        error "directory does not exist: ${dir}"
    }
    dirPath.listFiles().collectEntries { p -> [(p.getName()): p] }
}

// Resolve a genome FASTA for BFD's genome-stats channels, which need to work
// both post-annotation (genome_dir holds SETUP_SYMLINKS' <meta.id>.scaffolds.fa)
// and pre-annotation (genome_dir = input_clean_genomes/, holding raw/clean
// assemblies named by ASMID — see ANI_SAMPLES.nf, which uses this same
// <asmid>.fa.gz / <asmid>.masked.fasta.gz convention). Tries the annotated-name
// form first so an already-annotated genome keeps using its canonical output;
// falls back to the ASMID form so genome_stats can run directly off assemblies
// with no annotation yet. Returns null (not a non-existent file) when nothing
// on disk matches either convention.
//
// Takes the pre-built dirIndex(genomeDir) rather than genomeDir itself —
// callers that need this per-genome (BFD.nf's asm_stats_ch/busco_genome_ch)
// build the index once and pass it in, instead of paying for a listing (or
// worse, three .exists() calls) on every single genome.
def resolveGenomeFile(Map genomeIndex, String metaId, String asmid) {
    def byId = genomeIndex["${metaId}.scaffolds.fa"]
    if (byId) return byId
    if (!asmid) return null
    return genomeIndex["${asmid}.fa.gz"] ?: genomeIndex["${asmid}.masked.fasta.gz"]
}

// ── Gated globs and manifests ───────────────────────────────────────────────

// Gather files by filesystem glob (not by live channel .collect()), gated behind
// a completion signal.
//
// A live channel only contains items emitted *during this run's execution*. On a
// -resume run where most work is cached, only the actively-recomputed items flow
// through it, so a merge/combine step collecting off that channel silently loses
// everything else. Globbing the published files on disk instead makes the true
// input "everything ever computed", independent of what one invocation saw.
// (Found twice independently: BFD's MERGE_* steps, 2026-06-25, and
// COMBINE_ANI_TABLE's ani.duckdb, 2026-07-23 — see .living/learnings.md.)
def gatedGlobIn(sync_ch, String baseDir, String glob) {
    sync_ch
        .filter { it }
        .flatMap { files("${baseDir}/${glob}") }
}

// Collapse a channel of input files into one sorted manifest, TAB-delimited:
//     <absolute_path>\t<mtime_ms>\t<size_bytes>
//
// Consuming processes read the path from field 1 directly off its stable
// published location, instead of staging thousands of inputs via `path 'in/*'` —
// at ~8k genomes that staging command grew enormous (risking "Argument list too
// long") and created a symlink per file.
//
// Baking mtime+size into the manifest *content* also gives the staleness check
// for free: if any input is regenerated the manifest content changes, so the
// task's input hash changes and -resume re-runs it; if nothing changed the
// manifest is byte-identical and the task is a cache hit.
def toManifest(ch, String name) {
    ch.map { f -> "${f}\t${file(f).lastModified()}\t${file(f).size()}" }
      .collectFile(name: name, newLine: true, sort: true)
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
