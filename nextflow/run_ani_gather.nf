// run_ani_gather.nf — run after run_ani_compute.nf, from anywhere with S3
// access. Deliberately NOT k8s-specific: this is cheap (CSV parsing + small
// text files), so there's no reason to run it on the cluster at all — a
// laptop or any server with `nextflow`, the same samples.csv, and the S3
// credentials/endpoint config works fine (plain 'local' executor).
//
// Lives at the repo's nextflow/ root, sibling to main.nf — not in k8/ — for
// the same ${projectDir}/bin/*.py reason as run_ani_compute.nf (verified
// live: REPORT_ANI/COMBINE_ANI_TABLE's bin/ calls resolve to the wrong
// directory one level down).
//
// Re-derives each group's genome list fresh from samples.csv (same
// ANI_SAMPLES subworkflow the compute phase used — cheap, no genome file
// dependency) rather than depending on anything the compute phase published
// beyond the *.ani.tsv files themselves, then runs REPORT_ANI +
// COMBINE_ANI_TABLE against whatever's actually landed in
// ${outdir}/${ani_method}/${compare}/ so far.
//
// Usage:
//   nextflow run run_ani_gather.nf -c nextflow.config \
//       -params-file /path/to/params_ani.yaml
// (add `-profile compare_ani_k8s` only if you want this step's own S3
// endpoint/credentials config too — it doesn't need the k8s executor or PVC
// parts of that profile, just the aws{} block.)

include { assertRank; writeNamesTsv; gatedGlobIn; toManifest; asmidRowFilter } from './modules/common/utils.nf'
include { ANI_SAMPLES }       from './subworkflows/local/ANI_SAMPLES.nf'
include { REPORT_ANI }        from './modules/ani/report/REPORT_ANI/main.nf'
include { COMBINE_ANI_TABLE } from './modules/ani/report/COMBINE_ANI_TABLE/main.nf'

workflow {
    def compareRank = assertRank(params.compare as String, 'compare')

    ANI_SAMPLES(params.samples, compareRank, '', asmidRowFilter())

    def grouped_ch = ANI_SAMPLES.out.samples
        .groupTuple()
        .filter { _gname, metas -> metas.size() >= params.min_group_size as int }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    def names_map = grouped_ch.map { group_name, metas -> tuple(group_name, writeNamesTsv(group_name, metas)) }

    // Pick up whatever run_ani_compute.nf already published for this group —
    // glob rather than a fixed filename, since the exact name varies by
    // --ani_method (skani: <group>.full.ani.tsv under batches/; others differ).
    // Groups run_ani_compute.nf hasn't finished yet are silently skipped here
    // (that's the point — you can gather partial results and rerun later).
    def ani_tsv_ch = names_map
        .map { group_name, _names -> tuple(group_name, files("${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/**/*.ani.tsv")) }
        .filter { _gn, matches -> matches }
        .map { gn, matches -> tuple(gn, matches[0]) }

    REPORT_ANI(ani_tsv_ch.join(names_map, by: 0))

    def ani_sync   = ani_tsv_ch.collect().map { true }.ifEmpty(true)
    def names_sync = REPORT_ANI.out.names.collect().map { true }.ifEmpty(true)

    def ani_all   = gatedGlobIn(ani_sync, params.outdir, "${params.ani_method}/${params.compare}/**/*.ani.tsv")
    def names_all = gatedGlobIn(names_sync, params.outdir, "${params.ani_method}/${params.compare}/**/*_genome_names.tsv")

    COMBINE_ANI_TABLE(
        toManifest(ani_all,   'ani_tsv.manifest.txt'),
        toManifest(names_all, 'names_tsv.manifest.txt')
    )
}
