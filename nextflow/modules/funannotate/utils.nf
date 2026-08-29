//
// Helper functions shared across the funannotate pipeline.
//
// funannotate-specific (they read params.target, launchDir/rnaseq_reads, the
// shared ab-initio parameter store), so they live here rather than in
// modules/common/utils.nf, which is cross-pipeline.
//
// Plain functions, not a lib/ class -- see REFACTOR_NEXTFLOW_PLAN.md section 2.1.
//

// A funannotate step's GenBank output may be stored uncompressed (.gbk) or
// gzip-compressed (.gbk.gz) so completed folders can be compressed to save space.
// Returns the existing non-empty file (preferring .gbk), or null if neither exists.
// Use this for completion/skip gating so a compressed result still counts as "done".
def gbkResult(String dir, String out) {
    def plain = file("${dir}/${out}.gbk")
    if (plain.exists() && plain.size() > 0) return plain
    def gz = file("${dir}/${out}.gbk.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return null
}

// Clean/masked genomes in input_clean_genomes may be stored gzip-compressed (.gz) to
// save space. Given the uncompressed base path (e.g. .../<asmid>.fa or .../<asmid>.masked.fasta),
// returns the existing non-empty file, preferring the compressed form. Falls back to the
// plain path object when neither exists, so callers' .exists() checks still report missing.
def genomeFile(String base) {
    def gz = file("${base}.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return file(base)
}

def staleRnaseq(String out, String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def gbk = gbkResult("${params.target}/${out}/predict_results", out)
    if (gbk == null) return false  // predict hasn't run yet; normal path handles it
    // Ordinary species: real reads + representative-built Trinity-GG under the paths
    // RNASEQ_PREPARE/SRA_FETCH use. Hybrid-cross species (see
    // nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md) never populate these -- they
    // use rnaseq_reads/hybrid_empty/ (always 0-byte, so *_newer below is always false
    // for them, which is correct: empty placeholders never go stale) and a
    // composite-parents.trinity-GG.fasta instead of a plain one. Checking both path
    // sets unconditionally (rather than threading a hybrid flag through every one of
    // this function's several call sites) is safe and cheap -- for any given species
    // only one set will ever actually exist on disk.
    def r1      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_R1.fastq.gz")
    def r2      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_R2.fastq.gz")
    def se      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_SE.fastq.gz")
    def trinity = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
    def compositeTrinity = file("${launchDir}/rnaseq_data/${species_tag}.composite-parents.trinity-GG.fasta")
    def r1_newer        = r1.exists()      && r1.size() > 0      && r1.lastModified()      > gbk.lastModified()
    def r2_newer        = r2.exists()      && r2.size() > 0      && r2.lastModified()      > gbk.lastModified()
    def se_newer        = se.exists()      && se.size() > 0      && se.lastModified()      > gbk.lastModified()
    def trinity_newer   = trinity.exists() && trinity.size() > 0 && trinity.lastModified() > gbk.lastModified()
    def composite_newer = compositeTrinity.exists() && compositeTrinity.size() > 0 && compositeTrinity.lastModified() > gbk.lastModified()
    if (r1_newer || r2_newer || se_newer || trinity_newer || composite_newer) {
        log.info "stale prediction for ${out}: rnaseq/trinity newer than GBK — scheduling retrain+repredict"
        return true
    }
    return false
}

// Absolute path to a genome's source assembly archive, e.g.
// ${params.source}/GCA_000000000.1/GCA_000000000.1_genomic.fna.gz -- mirrors the
// construction in workflows/funannotate.nf's jobs channel so staleGenome() checks
// the same file that FUNANNOTATE_TRAIN/FUNANNOTATE_PREDICT actually consume.
def genomeSourceFile(String asmid) {
    file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
}

// A genome whose source assembly has been swapped/updated (e.g. a new assembly
// version dropped into the same asmid path) needs both retraining and repredicting --
// otherwise a re-assembled genome silently keeps a GBK annotated against stale
// coordinates. Mirrors staleRnaseq's shape/purpose but keys off the assembly file
// itself rather than the rnaseq inputs.
def staleGenome(String out, String asmid) {
    def gbk = gbkResult("${params.target}/${out}/predict_results", out)
    if (gbk == null) return false  // predict hasn't run yet; normal path handles it
    def gfa = genomeSourceFile(asmid)
    if (gfa.exists() && gfa.size() > 0 && gfa.lastModified() > gbk.lastModified()) {
        log.info "stale prediction for ${out}: genome assembly newer than GBK — scheduling retrain+repredict"
        return true
    }
    return false
}

// ── Species-level ab-initio parameter reuse (todo/species_level_abinitio_reuse.md) ──
// Backfilled/refreshed out-of-band by nextflow/bin/species_reuse_clusters.py, which
// writes both abinitio_reuse_csv and the per-species
// params.gene_prediction_shared_abinitio/<species_tag>/parameters.json stores this
// reads. Not generated by this workflow. Lives as a top-level sibling of
// params.target (not nested inside it) since it's shared across every project tree
// that annotates the same species, not scoped to one genome_annotation/ run.

// Absolute path to a species' shared ab-initio parameters.json, or null if the
// species has no backfilled store yet (species_reuse_clusters.py hasn't run for it,
// or its representative's own PREDICT never completed -- see plan S4.3 "representative
// failure/skip fallback": no store on disk is the correct, safe signal to fall back
// to independent training, not an error).
def sharedParamsJsonFor(String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def json = file("${params.gene_prediction_shared_abinitio}/${species_tag}/parameters.json")
    return (json.exists() && json.size() > 0) ? json : null
}

// Same existence/non-empty check as sharedParamsJsonFor(), pointed at the
// species' shared GeneMark .mod instead of parameters.json -- used by
// GENEMARK_RUN to decide fast-reuse (--predict_with) vs fresh --ES training.
// Backfilled by backfill_abinitio_params.py into the same per-species store,
// same file this pipeline already produced (predict_misc/ab_initio_parameters/
// before GENEMARK_RUN existed; GENEMARK_RUN.out.mod now, see design doc).
def sharedGenemarkModFor(String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def mod = file("${params.gene_prediction_shared_abinitio}/${species_tag}/${species_tag}.genemark.mod")
    return (mod.exists() && mod.size() > 0) ? mod : null
}

// A strain's FUNANNOTATE_TRAIN-produced transcript-to-genome alignment BAM --
// GENEMARK_RUN's ET mode derives its RNA-seq-informed intron hints from this
// (bam2hints -> filterIntronsFindStrand.pl -> join_mult_hints.pl -> gmes_petap.pl
// --ET; see nextflow/docs/GENEMARK_RUN_DESIGN.md's "ET mode" section for why
// this file specifically, not the raw RNA-seq read BAM). Returns '' (not null)
// when absent/empty -- GENEMARK_RUN's own script checks for a non-empty string,
// and an absent file (no RNA-seq for this genome, or FUNANNOTATE_TRAIN hasn't
// run) is the correct, safe signal to fall back to ES, not an error.
def trainingTranscriptBamFor(String out) {
    def bam = file("${params.training_target}/${out}/training/transcript.alignments.bam")
    return (bam.exists() && bam.size() > 0) ? bam.toString() : ''
}

// A reuse_eligible strain's GBK must be considered stale if the shared parameters.json
// it was predicted with has since been refreshed (representative re-annotated, or the
// ANI/reuse assignment changed) -- mirrors staleRnaseq's shape/purpose.
def staleSharedParams(String out, def sharedJson) {
    if (sharedJson == null) return false
    def gbk = gbkResult("${params.target}/${out}/predict_results", out)
    if (gbk == null) return false  // predict hasn't run yet; normal path handles it
    if (sharedJson.lastModified() > gbk.lastModified()) {
        log.info "stale prediction for ${out}: shared ab-initio parameters newer than GBK — scheduling repredict"
        return true
    }
    return false
}

// Eagerly loads abinitio_reuse_csv into out -> [species, reuse_eligible,
// is_representative]. Returns an empty map if the feature is off or the CSV doesn't
// exist yet (e.g. species_reuse_clusters.py hasn't been run for this dataset) --
// every row then falls through to today's independent-training behavior unchanged.
def loadAbinitioReuseMap() {
    def m = [:]
    if (!params.share_abinitio_params.toBoolean()) return m
    def csv = file(params.abinitio_reuse_csv as String)
    if (!csv.exists()) {
        log.warn "share_abinitio_params=true but abinitio_reuse_csv not found: ${csv} — all strains train independently"
        return m
    }
    def lines = csv.readLines()
    if (lines.size() < 2) return m
    def header = lines[0].split(',', -1)*.trim()
    def iSpecies      = header.indexOf('species')
    def iOut          = header.indexOf('out')
    def iEligible     = header.indexOf('reuse_eligible')
    def iIsRep        = header.indexOf('is_representative')
    def iAni          = header.indexOf('ani_to_representative')
    lines.drop(1).each { line ->
        def f = line.split(',', -1)
        if (f.size() <= [iSpecies, iOut, iEligible].max()) return
        // ani_to_representative is blank for the representative's own row (implicitly
        // 100%) and, notably, for any strain skani/mash/sourmash couldn't confidently
        // compare against the representative -- for sketch-based ANI tools that usually
        // means the true identity fell below their reliable detection floor (~80-85%),
        // i.e. a blank value is itself a signal of exceptional divergence, not "not yet
        // computed". Parsed as null here so callers can distinguish "no signal" from 0.0.
        def aniStr = iAni >= 0 && f.size() > iAni ? f[iAni].trim() : ''
        def isRep = iIsRep >= 0 ? f[iIsRep].trim().toLowerCase() == 'true' : false
        m[f[iOut].trim()] = [
            species         : f[iSpecies].trim(),
            reuse_eligible  : f[iEligible].trim().toLowerCase() == 'true',
            is_representative: isRep,
            ani_to_representative: isRep ? 100.0 : (aniStr ? aniStr as Double : null),
        ]
    }
    log.info "Loaded ${m.size()} ab-initio reuse assignments from ${csv} " +
        "(${m.count { it.value.reuse_eligible }} reuse_eligible, " +
        "${m.count { it.value.is_representative }} representative)"
    return m
}

// Eagerly loads rnaseq_representative_override_csv into species_tag -> out. Lets
// scripts/pick_rnaseq_representative_override.py force RNASEQ_PREPARE to build the
// shared Trinity-GG assembly against a specific strain instead of whatever
// abinitioReuseMap's ANI+BUSCO pick happens to be. Exists because that ANI+BUSCO pick
// optimizes for genome assembly quality, not for which genome the RNA-seq reads
// actually align to -- e.g. Ascochyta_rabiei's ANI-picked representative (a
// pks1-deletion construct genome, BUSCO 98.1%) is a perfectly good assembly that the
// species' real RNA-seq (SRR330019xx) barely aligns to (HISAT2 -> 9 Trinity-GG
// clusters -> 7 transcripts total), while GCF_004011695.2 (Me14, BUSCO 98.4%, the
// RefSeq reference for this species) is not the ANI pick but works fine. Returns an
// empty map when no override file exists -- every species then keeps using
// abinitioReuseMap's pick unchanged.
def loadRnaseqRepresentativeOverride() {
    def m = [:]
    def csv = file(params.rnaseq_representative_override_csv as String)
    if (!csv.exists()) return m
    def lines = csv.readLines()
    if (lines.size() < 2) return m
    def header = lines[0].split(',', -1)*.trim()
    def iTag = header.indexOf('species_tag')
    def iOut = header.indexOf('out')
    if (iTag < 0 || iOut < 0) {
        log.warn "rnaseq_representative_override_csv ${csv} missing species_tag/out columns -- ignoring"
        return m
    }
    lines.drop(1).each { line ->
        if (!line.trim()) return
        def f = line.split(',', -1)
        if (f.size() <= [iTag, iOut].max()) return
        m[f[iTag].trim()] = f[iOut].trim()
    }
    log.info "Loaded ${m.size()} RNA-seq representative overrides from ${csv}"
    return m
}

// Eagerly loads hybrid_parentage_csv into hybrid_species_tag -> [parent_species, ...].
// See nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md. Narrow/long format -- one row
// per (hybrid, parent) pair, so 2-way through N-way crosses need no schema change
// (parent_species is a plain "Genus species" string, matched against other species'
// species_tag via the same spaces->underscore rule used everywhere else in this
// codebase). Presence of a species_tag as a key IS the "this species is a hybrid"
// signal FUNANNOTATE_RNASEQ.nf acts on -- samples.csv carries no separate hybrid
// marker column (hybrids are <0.5% of rows; scripts/bootstrap_hybrid_metadata.py
// detects them by scanning SPECIES directly instead). Returns an empty map if the
// CSV doesn't exist yet (e.g. that bootstrap script hasn't been run) --
// every species then falls through to today's non-hybrid behavior unchanged.
def loadHybridParentage() {
    // Plain map, NOT [:].withDefault{[]} -- withDefault auto-vivifies on ANY read,
    // not just the append below, so a future bare `hybridParentage[someOtherTag]`
    // read anywhere in the codebase (as opposed to the `.containsKey()` check
    // FUNANNOTATE_RNASEQ.nf actually uses to test hybrid-ness) would silently
    // insert that tag as an empty-parent-list entry, reclassifying an ordinary
    // species as hybrid. `?: []` at each read call site already covers the
    // "not present" case without needing the map itself to auto-vivify.
    def m = [:]
    def csv = file(params.hybrid_parentage_csv as String)
    if (!csv.exists()) return m
    def lines = csv.readLines()
    if (lines.size() < 2) return m
    def header = lines[0].split(',', -1)*.trim()
    def iTag    = header.indexOf('hybrid_species_tag')
    def iParent = header.indexOf('parent_species')
    if (iTag < 0 || iParent < 0) {
        log.warn "hybrid_parentage_csv ${csv} missing hybrid_species_tag/parent_species columns -- ignoring"
        return m
    }
    lines.drop(1).each { line ->
        if (!line.trim()) return
        def f = line.split(',', -1)
        if (f.size() <= [iTag, iParent].max()) return
        def parentSpecies = f[iParent].trim()
        def parentTag = parentSpecies.replaceAll(/\s+/, '_')
        def hybridTag = f[iTag].trim()
        m.computeIfAbsent(hybridTag) { [] }
        m[hybridTag] << parentTag
    }
    log.info "Loaded hybrid parentage for ${m.size()} hybrid species from ${csv}"
    return m
}
