//
// FILTER_FALLBACK_PAIRS — from a fallback group's full all-vs-all ANI matrix,
// keep only orphan-vs-reference pairs (skip orphan-orphan and reference-reference).
//
// skani triangle (via SKANI_COMPARE) produces a group-wide matrix. In fallback mode
// the group contains both orphan genomes (no defined query_rank) and reference genomes
// (have defined query_rank). This process filters to the cross-type pairs only,
// then emits the same 3-column format (query, reference, ANI) as the standard
// query_ANI skani dist output.
//
process FILTER_FALLBACK_PAIRS {
    tag   "${group_name} [${orphan_metas.size()} orphans]"
    label 'report'

    input:
        tuple val(group_name), path(ani_tsv), val(orphan_metas)

    output:
        tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    def orphan_ids = orphan_metas.collect { m -> m.id }
    """
    awk -F'\\t' -v orphans="${orphan_ids.join(',')}" '
    BEGIN { n = split(orphans, arr, ","); for (i=1;i<=n;i++) set[arr[i]] = 1 }
    NR > 1 {
        q = \$1; r = \$2
        is_orphan_q = (q in set)
        is_orphan_r = (r in set)
        # Keep pairs where exactly one side is an orphan (xor)
        if (is_orphan_q != is_orphan_r) print
    }' "${ani_tsv}" > "${group_name}.query.ani.tsv"
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > "${group_name}.query.ani.tsv"
    """
}
