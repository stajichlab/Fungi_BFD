//
// ANI_SAMPLES — canonical samples.csv → per-genome channel for the ANI workflows.
//
// compare_ANI and query_ANI previously carried near-identical copies of this
// parsing block, differing only in whether an `is_query` flag was computed.
// The flag is always computed here; compare_ANI simply ignores it.
//
// Emits one entry per genome that exists on disk:
//     tuple( group_key, meta )
//
// where meta is a map (REFACTOR_NEXTFLOW_PLAN.md §2.4):
//     [ id, locustag, species, genus, strain, asmid, genome, is_query ]
//
// meta.genome is the resolved genome Path. It travels inside the map rather
// than as a separate `path` element because these are plain channel values —
// nothing is staged until a process input claims it — and keeping one payload
// per genome means group/ungroup operations don't have to re-zip parallel lists.
//

include { assertRank; taxonRowFilter; genomeStem } from '../../modules/common/utils.nf'

workflow ANI_SAMPLES {
    take:
    samples_csv     // path to samples.csv
    compare_rank    // validated, upper-case rank to group by
    query_rank      // upper-case rank whose absence marks a query genome, or '' to disable

    main:
    def nameStyle = (params.genome_name_style as String).toLowerCase()
    if (!(nameStyle in ['species', 'asmid'])) {
        error "--genome_name_style must be 'species' or 'asmid'"
    }

    def rowFilter = taxonRowFilter()

    samples_ch = channel
        .fromPath(samples_csv)
        .splitCsv(header: true)
        .filter(rowFilter)
        .map { row ->
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def groupKey = row[compare_rank]?.trim() ?: ''
            def stem     = genomeStem(row, nameStyle)
            def genome   = file("${params.genome_dir}/${stem}${params.genome_suffix}", glob: false)

            if (!groupKey) {
                log.warn "Skipping ${locustag}: empty ${compare_rank} field"
                return null
            }
            if (!genome.exists()) {
                log.warn "Skipping ${locustag}: genome not found at ${genome}"
                return null
            }

            def meta = [
                id      : genome.getName(),
                locustag: locustag,
                species : row.SPECIES?.trim() ?: '',
                genus   : row.GENUS?.trim()   ?: '',
                strain  : row.STRAIN?.trim()  ?: '',
                asmid   : row.ASMID?.trim()   ?: '',
                genome  : genome,
                // With no query_rank configured every genome is a reference.
                is_query: query_rank ? !(row[query_rank]?.trim()) : false,
            ]
            tuple(groupKey, meta)
        }
        .filter { item -> item != null }

    emit:
    samples = samples_ch
}
