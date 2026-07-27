//
// query_ANI — asymmetric orphan-vs-reference ANI search.
//
// Complements compare_ANI's symmetric all-vs-all triangle. Where a full clade
// is too large to triangulate just to place a handful of taxonomically-unassigned
// genomes, this workflow computes ANI ONLY between "query" genomes (samples.csv
// rows missing --query_rank, default GENUS) and "reference" genomes (everyone
// else in the same --compare group) via `skani dist`. Cost is O(queries x
// references) instead of O(references^2).
//
// A --compare group is processed only if it contains >=1 query genome AND
// >=min_group_size reference genomes; groups with no orphans are skipped
// automatically, so no --taxon restriction is required.
//
// --fallback mode:
// For genomes with no defined --query_rank (e.g., blank GENUS), walks up the
// taxonomy to the next defined rank (FAMILY → ORDER → CLASS …). Groups all
// orphans sharing the same fallback rank, sketches the combined set, runs
// skani triangle, then filters output to keep only orphan-vs-reference pairs.
// This places unassigned genomes by searching against all reference genomes
// within their nearest defined taxonomic group.
//
// Usage (from project root):
//   nextflow run nextflow/main.nf --pipeline query_ani \
//       -c nextflow/nextflow.config -profile ani --compare FAMILY -resume
//
//   # Fallback mode — genomes missing GENUS are placed by FAMILY/ORDER/etc.:
//   nextflow run nextflow/main.nf --pipeline query_ani \
//       -c nextflow/nextflow.config -profile ani --compare FAMILY --fallback true -resume
//
// Parameter defaults live in conf/profile_ANI.config; the ani_query profile
// layers on the query-specific defaults (--compare CLASS, --query_rank GENUS).
//

include { assertRank; writeNamesTsv; fallbackRank } from '../modules/common/utils.nf'
include { taxonomicRanks }                            from '../modules/common/utils.nf'
include { ANI_SAMPLES }               from '../subworkflows/local/ANI_SAMPLES.nf'
include { SKANI_DIST_QUERY }          from '../modules/ani/compare/SKANI_DIST_QUERY/main.nf'
include { SKANI_COMPARE }             from '../modules/ani/compare/SKANI_COMPARE/main.nf'
include { FILTER_FALLBACK_PAIRS }     from '../modules/ani/compare/FILTER_FALLBACK_PAIRS/main.nf'
include { REPORT_QUERY_ANI }          from '../modules/ani/report/REPORT_QUERY_ANI/main.nf'
include { COMBINE_QUERY_CALLS }       from '../modules/ani/report/COMBINE_QUERY_CALLS/main.nf'

workflow QUERY_ANI {

    def compareRank = assertRank(params.compare as String,    'compare')
    def queryRank   = assertRank(params.query_rank as String, 'query_rank')
    def useFallback = params.fallback as Boolean
    def ranks       = taxonomicRanks()
    if (ranks.indexOf(queryRank) <= ranks.indexOf(compareRank)) {
        error "--query_rank (${queryRank}) must be a narrower rank than --compare (${compareRank})"
    }

    log.info "query_ANI: grouping by ${compareRank}, querying genomes missing ${queryRank}" +
        (useFallback ? " [fallback mode: grouping by next defined rank]" : "")

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

    if (useFallback) {
        def fallback_ch = grouped_ch
            .flatMap { group_name, orphans, refs ->
                orphans.collect { m ->
                    def fbRank = fallbackRank(m, queryRank)
                    fbRank ? tuple(fbRank, m, group_name, refs) : null
                }
            }
            .filter { item -> item != null }
            .groupTuple()
            .map { fb_name, orphan_list, _gn_list, ref_lists ->
                def group_refs = ref_lists.flatten().unique { mm -> mm.id }
                tuple(fb_name, orphan_list, group_refs)
            }
            .filter { _fb_name, orphans, refs ->
                orphans.size() >= 1 && refs.size() >= (params.min_group_size as int)
            }

        fallback_ch.subscribe { gn, q, r ->
            log.info "  fallback ${gn}: ${q.size()} orphan(s) × ${r.size()} reference(s)"
        }

        def compare_in = fallback_ch.map { fb_name, orphans, refs ->
            def all      = orphans + refs
            def labelled = orphans.collect { mm -> mm + [role: 'query'] } +
                           refs.collect    { mm -> mm + [role: 'reference'] }
            tuple(
                fb_name,
                all.collect { mm -> mm.genome },
                writeNamesTsv(fb_name, labelled, 'query_names'),
                orphans
            )
        }

        def ani_raw_ch = SKANI_COMPARE(compare_in.map { rec -> tuple(rec[0], rec[1]) })
        def filtered_ch = FILTER_FALLBACK_PAIRS(
            ani_raw_ch.join(compare_in.map { rec -> tuple(rec[0], rec[3]) }, by: 0)
        )

        def report_in = filtered_ch
            .map { fb_name, pairs ->
                def names_file = compare_in
                    .filter { rec -> rec[0] == fb_name }
                    .map { rec -> rec[2] }
                    .first()
                tuple(fb_name, pairs, names_file)
            }

        def report_out = REPORT_QUERY_ANI(report_in)
        COMBINE_QUERY_CALLS(report_out[1].collect())

    } else {
        def prepared_ch = grouped_ch.map { group_name, queries, refs ->
            def qGenomes = queries.collect { m -> m.genome }
            def rGenomes = refs.collect    { m -> m.genome }
            def labelled = queries.collect { m -> m + [role: 'query'] } +
                           refs.collect    { m -> m + [role: 'reference'] }
            tuple(group_name, qGenomes, rGenomes, writeNamesTsv(group_name, labelled, 'query_names'))
        }

        def dist_in = prepared_ch.map { group_name, qGenomes, rGenomes, _nameFile ->
            tuple(group_name, qGenomes, rGenomes)
        }

        def ani_tsv_ch = SKANI_DIST_QUERY(dist_in)
        def report_out = REPORT_QUERY_ANI(ani_tsv_ch.join(prepared_ch.map { rec -> tuple(rec[0], rec[3]) }, by: 0))
        COMBINE_QUERY_CALLS(report_out[1].collect())
    }
}
