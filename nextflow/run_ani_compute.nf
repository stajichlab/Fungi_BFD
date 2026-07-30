// run_ani_compute.nf — the part meant to be distributed (e.g. across k8s).
//
// Lives at the repo's nextflow/ root, sibling to main.nf — not in k8/ — so
// that ${projectDir} inside REPORT_ANI/COMBINE_ANI_TABLE's bin/*.py calls
// resolves the same way it does when main.nf launches them (verified live:
// placing this one directory deeper broke it).
//
// Mirrors workflows/compare_ANI.nf's ANI_SAMPLES -> ANI_COMPARE_METHOD chain
// exactly (same subworkflows, so results are identical), but stops there: no
// REPORT_ANI, no COMBINE_ANI_TABLE. Each compare-method module already
// publishDirs its per-group *.ani.tsv, so that's the only side effect. Run
// run_ani_gather.nf afterward (from anywhere with S3 access — no k8s
// required) once all groups you care about have finished here.
//
// Usage (from a repo checkout on the PVC, inside the head pod):
//   nextflow run run_ani_compute.nf \
//       -c nextflow.config -profile compare_ani_k8s \
//       -params-file /path/to/params_ani.yaml -resume

include { assertRank; asmidRowFilter } from './modules/common/utils.nf'
include { ANI_SAMPLES }        from './subworkflows/local/ANI_SAMPLES.nf'
include { ANI_COMPARE_METHOD } from './subworkflows/local/ANI_COMPARE_METHOD.nf'

workflow {
    def compareRank = assertRank(params.compare as String, 'compare')
    def method = (params.ani_method as String).toLowerCase()
    if (!(method in ['skani', 'mash', 'sourmash', 'fastani'])) {
        error "--ani_method must be one of: skani, mash, sourmash, fastani"
    }

    log.info "ANI method: ${method}"
    log.info "Genome dir: ${params.genome_dir}"

    ANI_SAMPLES(params.samples, compareRank, '', asmidRowFilter())

    def grouped_ch = ANI_SAMPLES.out.samples
        .groupTuple()
        .filter { _gname, metas -> metas.size() >= params.min_group_size as int }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    def genome_ch = grouped_ch.map { group_name, metas -> tuple(group_name, metas.collect { m -> m.genome }) }

    ANI_COMPARE_METHOD(genome_ch, method)
}
