#!/usr/bin/env nextflow

include { makeSampleTag; sanitizeTag } from './modules/common/utils.nf'

/*
 * query_ANI — asymmetric orphan-vs-reference ANI search
 *
 * Complements compare_ANI.nf's symmetric all-vs-all triangle. Where a full
 * clade (e.g. all 5827 Sordariomycetes) is too large to triangulate just to
 * place a handful of taxonomically-unassigned genomes, this workflow computes
 * ANI ONLY between "query" genomes (samples.csv rows missing --query_rank,
 * default GENUS) and "reference" genomes (everyone else in the same
 * --compare group) via `skani dist`. Cost is O(queries x references) instead
 * of O(references^2) — the reason skani triangle over a whole CLASS/PHYLUM is
 * disproportionate just to place a few unassigned genomes.
 *
 * A --compare group is processed only if it contains >=1 query genome AND
 * >=min_group_size reference genomes; groups with no orphans are skipped
 * automatically, so no --taxon restriction is required to scope a run.
 *
 * Sketches are cached in the same --sketch_cache directory (and same
 * skani-parameter subpath) as compare_ANI.nf, so a genome sketched by either
 * workflow is reused by the other.
 *
 * Usage (from project root):
 *   nextflow run nextflow/query_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --compare FAMILY --samples samples.csv -resume
 *
 *   # scope to one clade while testing
 *   nextflow run nextflow/query_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --compare CLASS --taxon CLASS:Sordariomycetes -resume
 *
 *   # samples missing FAMILY (but with CLASS) against the rest of their class
 *   nextflow run nextflow/query_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --compare CLASS --query_rank FAMILY -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/query_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       -stub-run --n_test 3
 */

// ════════════════════════════════════════════════════════════════════════════
// PARAMETER DEFAULTS
// ════════════════════════════════════════════════════════════════════════════

params.samples                = "samples.csv"
params.genome_dir             = "${launchDir}/input/dna"
params.genome_suffix          = '.scaffolds.fa'
params.genome_name_style      = 'species'

// Group genomes for the search by this taxonomic rank (must be one step
// broader than --query_rank, e.g. compare CLASS while query_rank is FAMILY).
params.compare                = 'CLASS'
// A sample is a "query" (orphan) genome if this rank is empty in samples.csv.
// Everything else in the --compare group with this rank filled is a candidate
// reference.
params.query_rank             = 'GENUS'
params.taxon                  = ''

// Minimum number of REFERENCE genomes required in a group for the search to
// be worth running (query count only needs to be >=1).
params.min_group_size         = 2
// Keep only the top N reference hits per query genome (skani dist -n).
params.query_top_n            = 10

params.ani_cluster_threshold  = 95.0
params.ani_outlier_threshold  = 90.0

params.outdir                 = "${launchDir}/results/ANI"
// Shared with compare_ANI.nf — a genome sketched by either workflow is reused.
params.sketch_cache            = "${launchDir}/work/ANI/sketch_cache"
params.skani_preset            = 'medium'
params.skani_min_af            = 15
params.skani_compression       = 0
params.skani_sketch_chunk      = 50

params.n_test                  = 0

// ════════════════════════════════════════════════════════════════════════════
// HELPERS
// ════════════════════════════════════════════════════════════════════════════

def skaniPresetFlag(String p) {
    def flags = [fast: '--fast', medium: '', slow: '--slow']
    def key   = (p ?: 'medium').toLowerCase()
    if (!flags.containsKey(key)) {
        error "--skani_preset must be one of: fast, medium, slow"
    }
    flags[key]
}

// ════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ════════════════════════════════════════════════════════════════════════════

// Identical sketching contract to compare_ANI.nf's SKANI_SKETCH, including the
// storeDir path, so sketches are shared between the two workflows.
process SKANI_SKETCH {
    tag   "${group_key} [${genomes.size()} genomes]"
    label 'skani'
    cpus   4
    memory '8 GB'

    storeDir "${params.sketch_cache}/skani/${params.skani_preset}_af${params.skani_min_af}"

    input:
        tuple val(group_key), path(genomes), val(sketch_names)

    output:
        tuple val(group_key), path(sketch_names)

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    """
    rm -rf sk_out
    printf '%s\\n' ${genomes} > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sk_out
    mv sk_out/*.sketch .
    """

    stub:
    """
    for s in ${sketch_names.join(' ')}; do touch \$s; done
    """
}

// Asymmetric query-vs-reference ANI: O(queries x references), not O(refs^2).
process SKANI_DIST_QUERY {
    tag   "${group_name} q=${query_sketches.size()} r=${ref_sketches.size()}"
    label 'skani'

    cpus   { ref_sketches.size() > 2000 ? 32 : ref_sketches.size() > 500 ? 16 : 8 }
    memory { ref_sketches.size() > 2000 ? '64 GB' : ref_sketches.size() > 500 ? '32 GB' : '16 GB' }

    // publishDir (not storeDir): storeDir only checks output-path existence, so it
    // would never re-run a group's comparison after new query/reference genomes are
    // added later -- see compare_ANI.nf's SKANI_COMPARE for the same fix and the
    // .living/learnings.md entry that found this.
    publishDir { "${params.outdir}/skani_query/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name), path(query_sketches, stageAs: 'query/*'), path(ref_sketches, stageAs: 'ref/*')

    output:
        tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    """
    ls query/*.sketch > query_list.txt
    ls ref/*.sketch   > ref_list.txt
    skani dist --ql query_list.txt --rl ref_list.txt \\
        --min-af ${params.skani_min_af} -n ${params.query_top_n} -t ${task.cpus} \\
        -o skani_raw.tsv

    # skani cols: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query [names]
    awk -F'\\t' 'NR>1 {
        nr=split(\$1,b,"/"); r=b[nr]; sub(/\\.sketch\$/,"",r);
        nq=split(\$2,a,"/"); q=a[nq]; sub(/\\.sketch\$/,"",q);
        print q"\\t"r"\\t"\$3
    }' skani_raw.tsv > ${group_name}.query.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.query.ani.tsv
    """
}

process REPORT_QUERY_ANI {
    tag   "${group_name}"
    label 'report'

    publishDir { "${params.outdir}/skani_query/${params.compare}/${group_name}" }, mode: 'copy'

    input:
        tuple val(group_name), path(ani_tsv), path(names_tsv)

    output:
        path("${group_name}_query_report.txt")
        path("${group_name}_query_calls.tsv")

    script:
    """
    python3 ${projectDir}/bin/report_query_ani.py \\
        --input  ${ani_tsv} \\
        --names  ${names_tsv} \\
        --group  "${group_name}" \\
        --level  "${params.compare}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --outlier-threshold ${params.ani_outlier_threshold} \\
        --top-n  ${params.query_top_n} \\
        --report ${group_name}_query_report.txt \\
        --calls  ${group_name}_query_calls.tsv
    """

    stub:
    """
    printf '=== ${group_name} (stub) ===\\n' > ${group_name}_query_report.txt
    printf 'group\\tlevel\\tquery_asmid\\tquery_species\\tquery_strain\\tbest_ref_asmid\\tbest_ref_genus\\tbest_ref_species\\tbest_ani\\tn_refs_ge_outlier\\ttier\\n' > ${group_name}_query_calls.tsv
    """
}

// Merge every group's calls into one master orphan-classification table.
process COMBINE_QUERY_CALLS {
    label 'report'

    publishDir "${params.outdir}/skani_query/${params.compare}", mode: 'copy'

    input:
        path calls_tsvs, stageAs: 'calls_in/*'

    output:
        path("orphan_calls.csv")

    script:
    """
    python3 ${projectDir}/bin/combine_query_calls.py --calls-dir calls_in --output orphan_calls.csv
    """

    stub:
    """
    touch orphan_calls.csv
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {

    def validRanks = ['PHYLUM','SUBPHYLUM','CLASS','SUBCLASS','ORDER','FAMILY','GENUS','SPECIES']
    def compareRank = (params.compare as String).toUpperCase()
    if (!(compareRank in validRanks)) {
        error "--compare must be one of: ${validRanks.join(', ')}"
    }
    def queryRank = (params.query_rank as String).toUpperCase()
    if (!(queryRank in validRanks)) {
        error "--query_rank must be one of: ${validRanks.join(', ')}"
    }
    if (validRanks.indexOf(queryRank) <= validRanks.indexOf(compareRank)) {
        error "--query_rank (${queryRank}) must be a narrower rank than --compare (${compareRank})"
    }

    def nameStyle = (params.genome_name_style as String).toLowerCase()
    if (!(nameStyle in ['species', 'asmid'])) {
        error "--genome_name_style must be 'species' or 'asmid'"
    }

    log.info "query_ANI: grouping by ${compareRank}, querying genomes missing ${queryRank}"

    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be RANK:VALUE, e.g. --taxon CLASS:Sordariomycetes"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { _row -> true }
    }

    // ── Build sample channel ──────────────────────────────────────────────────
    // Emits: tuple(group_key, locustag, species, genome_path, genus, strain, asmid, is_query)
    def sample_ch = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def species  = row.SPECIES?.trim() ?: ''
            def genus    = row.GENUS?.trim() ?: ''
            def strain   = row.STRAIN?.trim() ?: ''
            def asmid    = row.ASMID?.trim() ?: ''
            // Sanitised for use as group_name in file/dir names (e.g. ${group_name}.full.ani.tsv)
            // -- GENUS/FAMILY/etc. rarely need this, but SPECIES values ("Aspergillus
            // fumigatus") contain spaces that must not land raw in a shell-globbed filename.
            def groupKey = sanitizeTag(row[compareRank])
            def isQuery  = !(row[queryRank]?.trim())
            def stem     = nameStyle == 'asmid'
                              ? row.ASMID?.trim()
                              : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def genome   = file("${params.genome_dir}/${stem}${params.genome_suffix}", glob: false)

            if (!groupKey) {
                log.warn "Skipping ${locustag}: empty ${compareRank} field"
                return null
            }
            if (!genome.exists()) {
                log.warn "Skipping ${locustag}: genome not found at ${genome}"
                return null
            }
            tuple(groupKey, locustag, species, genome, genus, strain, asmid, isQuery)
        }
        .filter { item -> item != null }

    // ── Group by taxonomic rank, split into query/reference, drop groups with
    //    no orphans or too few references to search against ────────────────────
    def grouped_ch = sample_ch
        .map { groupKey, locustag, species, genome, genus, strain, asmid, isQuery ->
            tuple(groupKey, tuple(locustag, species, genome, genus, strain, asmid, isQuery))
        }
        .groupTuple()
        .map { group_name, members ->
            def queries = members.findAll { m -> m[6] }
            def refs    = members.findAll { m -> !m[6] }
            tuple(group_name, queries, refs)
        }
        .filter { _gn, queries, refs -> queries.size() >= 1 && refs.size() >= (params.min_group_size as int) }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    grouped_ch.subscribe { gn, q, r ->
        log.info "${gn}: ${q.size()} query genome(s) vs ${r.size()} reference genome(s)"
    }

    // ── Per-group genome lists + names TSV (adds a role column) ───────────────
    def prepared_ch = grouped_ch
        .map { group_name, queries, refs ->
            def qGenomes = queries.collect { m -> m[2] }
            def rGenomes = refs.collect    { m -> m[2] }
            // names TSV columns: filename, asmid, genus, species, strain, role
            def rows = queries.collect { m -> "${m[2].getName()}\t${m[5]}\t${m[3]}\t${m[1]}\t${m[4]}\tquery" } +
                       refs.collect    { m -> "${m[2].getName()}\t${m[5]}\t${m[3]}\t${m[1]}\t${m[4]}\treference" }
            def nameFile = file("${workflow.workDir}/query_names_${group_name}.tsv")
            nameFile.text = "filename\tasmid\tgenus\tspecies\tstrain\trole\n" + rows.join('\n') + "\n"
            tuple(group_name, qGenomes, rGenomes, nameFile)
        }

    def prepared_split = prepared_ch.multiMap { group_name, qGenomes, rGenomes, nameFile ->
        query: tuple(group_name, qGenomes)
        ref:   tuple(group_name, rGenomes)
        names: tuple(group_name, nameFile)
    }

    // ── Sketch queries and references separately (grouped by [group,role] so
    //    the two roles can be recombined per-group after sketching) ────────────
    def chunk = Math.max(1, params.skani_sketch_chunk as int)

    def query_sketch_in = prepared_split.query.flatMap { gn, genomes ->
        genomes.collate(chunk).collect { sub -> tuple(tuple(gn, 'query'), sub, sub.collect { g -> "${g.name}.sketch" }) }
    }
    def ref_sketch_in = prepared_split.ref.flatMap { gn, genomes ->
        genomes.collate(chunk).collect { sub -> tuple(tuple(gn, 'ref'), sub, sub.collect { g -> "${g.name}.sketch" }) }
    }

    def sketched_ch = SKANI_SKETCH(query_sketch_in.mix(ref_sketch_in))
        .groupTuple()
        .map { key, sketch_lists -> tuple(key[0], key[1], sketch_lists.flatten()) }

    def sketched_split = sketched_ch.branch {
        query: it[1] == 'query'
        ref:   it[1] == 'ref'
    }
    def query_sketches_ch = sketched_split.query.map { gn, _role, sk -> tuple(gn, sk) }
    def ref_sketches_ch   = sketched_split.ref.map   { gn, _role, sk -> tuple(gn, sk) }

    def dist_in = query_sketches_ch.join(ref_sketches_ch, by: 0)

    def ani_tsv_ch = SKANI_DIST_QUERY(dist_in)

    // ── Report + combine ──────────────────────────────────────────────────────
    def report_out = REPORT_QUERY_ANI(ani_tsv_ch.join(prepared_split.names, by: 0))
    COMBINE_QUERY_CALLS(report_out[1].collect())
}
