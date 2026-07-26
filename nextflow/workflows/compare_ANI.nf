//
// compare_ANI — symmetric all-vs-all ANI within taxonomic groups.
//
// Groups samples.csv by --compare rank, then sketches and compares every genome
// in each group with --ani_method (skani | mash | sourmash | fastani).
// See workflows/query_ANI.nf for the asymmetric orphan-placement counterpart.
//
// Usage (from project root):
//   nextflow run nextflow/main.nf -entry COMPARE_ANI \
//       -c nextflow/nextflow.config -profile ani --compare GENUS -resume
//
// Parameter defaults live in conf/profile_ANI.config.
//

include { assertRank; writeNamesTsv } from '../modules/common/utils.nf'
include { ANI_SAMPLES }               from '../subworkflows/local/ANI_SAMPLES.nf'
include { ANI_COMPARE_METHOD }        from '../subworkflows/local/ANI_COMPARE_METHOD.nf'
include { REPORT_ANI }                from '../modules/ani/report/REPORT_ANI/main.nf'
include { COMBINE_ANI_TABLE }         from '../modules/ani/report/COMBINE_ANI_TABLE/main.nf'

workflow COMPARE_ANI {

    // ── Validate params ───────────────────────────────────────────────────────
    def compareRank = assertRank(params.compare as String, 'compare')

    def method = (params.ani_method as String).toLowerCase()
    if (!(method in ['skani','mash','sourmash','fastani'])) {
        error "--ani_method must be one of: skani, mash, sourmash, fastani"
    }

    log.info "ANI method: ${method}"
    log.info "Genome filename style: ${params.genome_name_style} (suffix: ${params.genome_suffix})"
    if (method == 'fastani' && (params.fastani_prefilter as boolean)) {
        log.info "fastANI mash-prefilter cascade ENABLED (prefilter_ani=${params.prefilter_ani}%)"
    }

    // ── Samples → groups ──────────────────────────────────────────────────────
    // query_rank is '' here: compare_ANI treats every genome as a reference.
    ANI_SAMPLES(params.samples, compareRank, '')

    // n_test limits *groups* (applied after groupTuple so --n_test 3 = 3 groups).
    def grouped_ch = ANI_SAMPLES.out.samples
        .groupTuple()
        .filter { _gname, metas -> metas.size() >= params.min_group_size as int }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // ── Per-group genome list + names TSV ─────────────────────────────────────
    // ANI_SAMPLES already dropped genomes absent from disk, so every meta here
    // has a real file; the names TSV and the genome list stay in agreement.
    def prepared_ch = grouped_ch.map { group_name, metas ->
        tuple(group_name, metas.collect { m -> m.genome }, writeNamesTsv(group_name, metas))
    }

    def prepared_split = prepared_ch.multiMap { group_name, genomes, nameFile ->
        ani:   tuple(group_name, genomes)
        names: tuple(group_name, nameFile)
    }
    def genome_ch = prepared_split.ani     // tuple(group_name, [genome files])
    def names_map = prepared_split.names   // tuple(group_name, names_file)

    // ── Sketch + compare ──────────────────────────────────────────────────────
    ANI_COMPARE_METHOD(genome_ch, method)
    def ani_tsv_ch = ANI_COMPARE_METHOD.out.ani_tsv

    // ── Reports ───────────────────────────────────────────────────────────────
    REPORT_ANI(ani_tsv_ch.join(names_map, by: 0))

    // ── Combine every group's pairs + labels into one queryable table ─────────
    COMBINE_ANI_TABLE(
        ani_tsv_ch.map { _gn, tsv -> tsv }.collect(),
        names_map.map  { _gn, tsv -> tsv }.collect()
    )
}
