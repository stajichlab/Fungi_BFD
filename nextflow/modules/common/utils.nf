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

// Resolve the output directory for the MERGE steps: always params.tables.
// Per T-014 §D.2, MERGE_* always builds one full/unscoped master table set now
// -- a --taxon run still only *processes* that clade (taxonRowFilter() at the
// base genome channel), it no longer also produces a separate
// tables/<Taxon>/ subset. (Previously branched on params.taxon into
// tables/All_Taxa vs tables/<sanitised taxon value>; retired along with
// MERGE_SAMPLES's per-taxon --matched filtering.) Post-hoc taxon-scoped
// extraction from the master DB is EXTRACT_BFD_TAXONOMIC_SUBSET (separate
// tool, not this function).
def tablesDir() {
    "${params.tables}"
}

// Resolve the output directory for BUILD_DUCKDB's BFD.duckdb (params.db,
// default "${launchDir}/db" -- see nextflow.config).
def dbDir() {
    "${params.db}"
}

// ── Sample identity ─────────────────────────────────────────────────────────
//
// Moved here from lib/SampleUtils.groovy: Groovy classes on the lib/ classpath
// are invisible to Nextflow's strict language spec (§2.1). The Python equivalent
// is fix_low_trinity.py::species_key_from_row() — keep them in sync.
// (collect_busco_stats.py::build_basename_map() used to mirror this too, but was
// removed when BUSCO_genome/BUSCO_protein switched to ASMID/LOCUSTAG-keyed
// filenames -- neither needs a Genus_species_strain basename lookup anymore,
// see todo/genome_stats_storage_reorg.md T-014.)

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

// Load params.suppress once into a Set of ASMIDs. One ASMID per line; blank
// lines and '#'-prefixed comments are ignored, and only the first
// comma-separated field is used so annotated lines (e.g. "GCA_000,reason")
// still parse. Returns an empty Set (never null) when params.suppress is
// unset or the file doesn't exist, so callers can use the result unconditionally.
def loadSuppressSet() {
    if (!params.suppress) {
        return ([] as Set)
    }
    def f = file(params.suppress as String)
    if (!f.exists()) {
        return ([] as Set)
    }
    def suppressSet = f.readLines()
        .collect { it.trim().split(',')[0].trim() }
        .findAll { it && !it.startsWith('#') }
        .toSet()
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }
    return suppressSet
}

// Build the set of `out` values (see makeSampleTag()) whose samples.csv
// PHYLUM matches params.microsporidia_phylum. Mirrors loadSuppressSet()'s
// synchronous-read pattern: read once, return a plain Set for O(1)
// .contains() checks in FUNANNOTATE_PREDICTION.nf's per-row closures.
//
// MUST key off SPECIES/STRAIN (via makeSampleTag), matching exactly how
// workflows/funannotate.nf builds `out` for the jobs channel (line 110) --
// NOT SPECIES_IN, a separate raw/unnormalized column that produces a
// different, non-matching tag.
//
// Plain comma-split (not a full CSV parser) matches every other synchronous
// loader in this file (loadSuppressSet, loadAbinitioReuseMap in
// modules/funannotate/utils.nf) -- if embedded commas in quoted fields ever
// bite here, fix it in all of them together, not just this one.
def loadMicrosporidiaOutSet() {
    def phylum = (params.microsporidia_phylum ?: '').trim()
    if (!phylum) {
        return ([] as Set)
    }
    def f = file(params.samples as String)
    if (!f.exists()) {
        return ([] as Set)
    }
    def lines = f.readLines()
    if (lines.size() < 2) {
        return ([] as Set)
    }
    def header  = lines[0].split(',', -1)*.trim()
    def iSpecies = header.indexOf('SPECIES')
    def iStrain  = header.indexOf('STRAIN')
    def iPhylum  = header.indexOf('PHYLUM')
    if (iSpecies < 0 || iStrain < 0 || iPhylum < 0) {
        log.warn "loadMicrosporidiaOutSet: samples.csv missing SPECIES/STRAIN/PHYLUM column; microsporidia_phylum gate is a no-op"
        return ([] as Set)
    }
    def outSet = lines.drop(1)
        .collect { it.split(',', -1) }
        .findAll { f2 -> f2.size() > [iSpecies, iStrain, iPhylum].max() }
        .findAll { f2 -> f2[iPhylum].trim().equalsIgnoreCase(phylum) }
        .collect { f2 -> makeSampleTag(f2[iSpecies].trim(), f2[iStrain].trim()) }
        .toSet()
    if (outSet) {
        log.info "loadMicrosporidiaOutSet: ${outSet.size()} genome(s) tagged PHYLUM=${phylum}"
    }
    return outSet
}

// Build the samples.csv row filter for the suppress list produced by
// loadSuppressSet(). Mirrors taxonRowFilter()/asmidRowFilter() above: pass
// its result straight to .filter() and it's a no-op when the set is empty.
def suppressRowFilter(Set suppressSet) {
    if (!suppressSet) {
        return { _row -> true }
    }
    return { row ->
        def asmid = row.ASMID?.trim()
        if (asmid && suppressSet.contains(asmid)) {
            log.info "Suppressing ASMID ${asmid} (suppress list)"
            return false
        }
        return true
    }
}

// ── ANI ─────────────────────────────────────────────────────────────────────

// Clamp cpus/memory(GB) to the k8s admission-webhook ceiling, when one is
// configured (0 = uncapped, the default everywhere except
// k8/profile_compare_ani_k8s.config). Verified live (2026-07-30): Nautilus's
// ucr-stajichlab namespace rejects any bare (uncontrolled) pod -- what nf-k8s
// submits for every task -- above 16 cores / 32 GB RAM ("admission webhook
// pod.nrp-nautilus.io denied the request"). Hit this for real on
// Saccharomyces (1547 genomes): skaniCpusFor/skaniMemoryFor wanted 32
// cores/100GB, the pod was rejected outright (HTTP 400, no exit code to even
// retry on), and the whole run aborted. HPC has no such restriction -- Slurm
// happily gives skani/fastani/etc. their full 32-64 core / 100-128GB
// requests on the epyc/highmem queues -- so this only clamps when
// k8s_max_cpus/k8s_max_mem_gb are actually set (i.e. under -profile
// compare_ani_k8s), never on Slurm.
def capCpus(int n) {
    def cap = (params.k8s_max_cpus ?: 0) as int
    cap > 0 ? Math.min(n, cap) : n
}

def capMemGB(int gb) {
    def cap = (params.k8s_max_mem_gb ?: 0) as int
    cap > 0 ? Math.min(gb, cap) : gb
}

// skani preset → CLI flag (medium is the tool default → no flag).
// Shared by every skani process; the sketch and the compare step MUST pass the
// same preset or the cached sketches are unusable.
def skaniCpusFor(int n) {
    capCpus(n > 500 ? 32 : n > 200 ? 16 : 8)
}

def skaniMemoryFor(int n, int attempt = 1) {
    def baseGB  = n > 500 ? 48 : n > 200 ? 32 : 16
    def floorGB = attempt >= 3 ? 100 : attempt == 2 ? 24 : 6
    capMemGB(Math.max(baseGB, floorGB)).toString() + ' GB'
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
    // toFile().listFiles(), not Path.listFiles() -- the latter is Nextflow's
    // Path extension and does not follow a symlinked directory (treats it as
    // "not a directory" and returns a 1-element list containing just the
    // symlink). See ANI_SAMPLES.nf's genomeIndex, same fix.
    dirPath.toFile().listFiles().collectEntries { p -> [(p.getName()): file(p.toString(), glob: false)] }
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
        .flatMap { files("${baseDir}/${glob}", type: 'file') }
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
    // glob: false -- f is already a resolved Path, not a pattern. Some asmid/
    // strain names contain glob metacharacters (e.g. a literal '*' from a
    // taxonomy field, seen in F0AEDD87_AH-3-21*), and plain file(f) silently
    // re-globs the string, returning a LinkedList instead of the file itself
    // ("Unknown method invocation `lastModified` on LinkedList type").
    ch.map { f -> "${f}\t${file(f, glob: false).lastModified()}\t${file(f, glob: false).size()}" }
      .toSortedList()
      .flatMap { it }
      .collectFile(name: name, newLine: true)
}

// ── genome_stats/function hash-bucket fan-out ───────────────────────────────
//
// genome_stats/function outputs are stored one file per genome (or one file
// per genome per data column), which reaches tens of thousands of files in a
// single flat directory (see todo/genome_stats_storage_reorg.md, T-014) --
// slow to list/back up on this cluster's NFS-backed storage. hashBucket()
// fans a stable key (ASMID or LOCUSTAG, never the derived Genus_species_strain
// tag -- see the plan for why) out into a hex subdirectory so no single
// directory holds more than a few hundred files.
//
// The Python equivalent is nextflow/bin/genome_stats_paths.py::hash_bucket().
// Both must produce identical output for identical (key, width) -- a mismatch
// here is a silent path miss, not a crash. Verified in sync via
// nextflow/tests/test_hash_bucket_parity.py, which runs this function through
// a throwaway Nextflow script and compares against the Python implementation
// for the same key set.
def hashBucket(String key, int width) {
    def digest = java.security.MessageDigest.getInstance('SHA-1').digest(key.getBytes('UTF-8'))
    def hex = digest.collect { b -> String.format('%02x', (b as int) & 0xff) }.join('')
    hex.substring(0, width)
}

// Per-type bucket width (todo/genome_stats_storage_reorg.md §A): most
// genome_stats types get 256 buckets (2 hex chars); gene_stats and
// intergenic_stats get 4096 buckets (3 hex chars) because they carry ~7x more
// files per genome than the rest. Unknown types (e.g. any function/* subfolder
// not yet migrated) default to width 2. Must stay in sync with
// nextflow/bin/genome_stats_paths.py::BUCKET_WIDTH.
def bucketWidthFor(String type) {
    def widths = [
        asm_stats       : 2,
        BUSCO_genome    : 2,
        BUSCO_protein   : 2,
        aa_freq         : 2,
        codon_freq      : 2,
        gene_stats      : 3,
        intergenic_stats: 3,
    ]
    widths.containsKey(type) ? widths[type] : 2
}

// Convenience wrapper: look up the right width for `type` and hash `key`.
def hashBucketForType(String type, String key) {
    hashBucket(key, bucketWidthFor(type))
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
