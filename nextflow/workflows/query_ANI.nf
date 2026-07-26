//
// query_ANI — asymmetric orphan-vs-reference ANI search.
//
// Complements compare_ANI's symmetric all-vs-all triangle. Where a full clade
// (e.g. all 5827 Sordariomycetes) is too large to triangulate just to place a
// handful of taxonomically-unassigned genomes, this workflow computes ANI ONLY
// between "query" genomes (samples.csv rows missing --query_rank, default GENUS)
// and "reference" genomes (everyone else in the same --compare group) via
// `skani dist`. Cost is O(queries x references) instead of O(references^2).
//
// A --compare group is processed only if it contains >=1 query genome AND
// >=min_group_size reference genomes; groups with no orphans are skipped
// automatically, so no --taxon restriction is required to scope a run.
//
// Sketches are cached in the same --sketch_cache directory (and same
// skani-parameter subpath) as compare_ANI, so a genome sketched by either
// workflow is reused by the other — SKANI_SKETCH is literally the same module.
//
// Usage (from project root):
//   nextflow run nextflow/main.nf -entry QUERY_ANI \
//       -c nextflow/nextflow.config -profile ani_query --compare FAMILY -resume
//
// Parameter defaults live in conf/profile_ANI.config; the ani_query profile
// layers on the query-specific defaults (--compare CLASS, --query_rank GENUS).
//

include { assertRank; writeNamesTsv } from '../modules/common/utils.nf'
include { taxonomicRanks }            from '../modules/common/utils.nf'
include { ANI_SAMPLES }               from '../subworkflows/local/ANI_SAMPLES.nf'
include { SKANI_SKETCH }              from '../modules/ani/sketch/SKANI_SKETCH/main.nf'
include { SKANI_DIST_QUERY }          from '../modules/ani/compare/SKANI_DIST_QUERY/main.nf'
include { REPORT_QUERY_ANI }          from '../modules/ani/report/REPORT_QUERY_ANI/main.nf'
include { COMBINE_QUERY_CALLS }       from '../modules/ani/report/COMBINE_QUERY_CALLS/main.nf'

// Fan a per-group genome list out into SKANI_SKETCH chunk jobs, keyed by
// (group, role) so query and reference sketches can be recombined per group
// after sketching. Declared at file scope: a closure assigned to a local
// variable inside a workflow body is rejected by the strict language spec.
def sketchInputs(ch, String role) {
    def chunk = Math.max(1, params.skani_sketch_chunk as int)
    ch.flatMap { gn, genomes ->
        genomes.collate(chunk).collect { sub ->
            tuple(tuple(gn, role), sub, sub.collect { g -> "${g.name}.sketch" })
        }
    }
}

workflow QUERY_ANI {

    // ── Validate params ───────────────────────────────────────────────────────
    def compareRank = assertRank(params.compare as String,    'compare')
    def queryRank   = assertRank(params.query_rank as String, 'query_rank')
    def ranks       = taxonomicRanks()
    if (ranks.indexOf(queryRank) <= ranks.indexOf(compareRank)) {
        error "--query_rank (${queryRank}) must be a narrower rank than --compare (${compareRank})"
    }

    log.info "query_ANI: grouping by ${compareRank}, querying genomes missing ${queryRank}"

    // ── Samples → groups, split into query / reference ────────────────────────
    ANI_SAMPLES(params.samples, compareRank, queryRank)

    def grouped_ch = ANI_SAMPLES.out.samples
        .groupTuple()
        .map { group_name, metas ->
            tuple(group_name, metas.findAll { m -> m.is_query }, metas.findAll { m -> !m.is_query })
        }
        .filter { _gn, queries, refs ->
            queries.size() >= 1 && refs.size() >= (params.min_group_size as int)
        }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    grouped_ch.subscribe { gn, q, r ->
        log.info "${gn}: ${q.size()} query genome(s) vs ${r.size()} reference genome(s)"
    }

    // ── Per-group genome lists + names TSV (adds a role column) ───────────────
    def prepared_ch = grouped_ch.map { group_name, queries, refs ->
        def labelled = queries.collect { m -> m + [role: 'query'] } +
                       refs.collect    { m -> m + [role: 'reference'] }
        tuple(group_name,
              queries.collect { m -> m.genome },
              refs.collect    { m -> m.genome },
              writeNamesTsv(group_name, labelled, 'query_names'))
    }

    def prepared_split = prepared_ch.multiMap { group_name, qGenomes, rGenomes, nameFile ->
        query: tuple(group_name, qGenomes)
        ref:   tuple(group_name, rGenomes)
        names: tuple(group_name, nameFile)
    }

    // ── Sketch queries and references separately (see sketchInputs above) ────
    def sketched_ch = SKANI_SKETCH(
            sketchInputs(prepared_split.query, 'query')
                .mix(sketchInputs(prepared_split.ref, 'ref'))
        )
        .groupTuple()
        .map { key, sketch_lists -> tuple(key[0], key[1], sketch_lists.flatten()) }

    def sketched_split = sketched_ch.branch {
        query: it[1] == 'query'
        ref:   it[1] == 'ref'
    }

    def dist_in = sketched_split.query.map { gn, _role, sk -> tuple(gn, sk) }
        .join(sketched_split.ref.map { gn, _role, sk -> tuple(gn, sk) }, by: 0)

    def ani_tsv_ch = SKANI_DIST_QUERY(dist_in)

    // ── Report + combine ──────────────────────────────────────────────────────
    def report_out = REPORT_QUERY_ANI(ani_tsv_ch.join(prepared_split.names, by: 0))
    COMBINE_QUERY_CALLS(report_out[1].collect())
}
