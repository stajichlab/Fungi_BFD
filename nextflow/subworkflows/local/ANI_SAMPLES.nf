//
// ANI_SAMPLES — canonical samples.csv → per-genome channel for the ANI workflows.
//
// compare_ANI and query_ANI previously carried near-identical copies of this
// parsing block, differing only in whether an `is_query` flag was computed.
// The flag is always computed here; compare_ANI simply ignores it.
//
// Genome lookup: reads ASMID from samples.csv and finds the matching .fa.gz or
// .masked.fasta.gz in params.genome_dir (defaults to input_clean_genomes/).
// skani reads gzipped FASTA directly — no decompression needed.
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

include { assertRank; taxonRowFilter } from '../../modules/common/utils.nf'
include { sanitizeTag } from '../../modules/common/utils.nf'

workflow ANI_SAMPLES {
    take:
    samples_csv     // path to samples.csv
    compare_rank    // validated, upper-case rank to group by
    query_rank      // upper-case rank whose absence marks a query genome, or '' to disable

    main:
    def rowFilter = taxonRowFilter()

    samples_ch = channel
        .fromPath(samples_csv)
        .splitCsv(header: true)
        .filter(rowFilter)
        .map { row ->
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def asmid    = row.ASMID?.trim() ?: ''
            // Sanitised for use as group_name in file/dir names (e.g.
            // ${group_name}.full.ani.tsv) — GENUS/FAMILY rarely need it, but SPECIES
            // values ("Aspergillus fumigatus") contain spaces that must not land raw
            // in a shell-globbed filename.
            def groupKey = sanitizeTag(row[compare_rank])

            if (!groupKey) {
                log.warn "Skipping ${locustag}: empty ${compare_rank} field"
                return null
            }
            if (!asmid) {
                log.warn "Skipping ${locustag}: empty ASMID field"
                return null
            }

            // Look for <asmid>.fa.gz or <asmid>.masked.fasta.gz in genome_dir.
            // skani reads .gz directly — no decompression step needed.
            def gz      = file("${params.genome_dir}/${asmid}.fa.gz",           glob: false)
            def masked  = file("${params.genome_dir}/${asmid}.masked.fasta.gz", glob: false)
            def genome  = gz.exists() ? gz : (masked.exists() ? masked : null)

            if (!genome) {
                log.warn "Skipping ${locustag} (${asmid}): no .fa.gz or .masked.fasta.gz in ${params.genome_dir}"
                return null
            }

            def meta = [
                id      : genome.getName(),
                locustag: locustag,
                phylum  : row.PHYLUM?.trim()    ?: '',
                'class' : row.CLASS?.trim()     ?: '',
                order   : row.ORDER?.trim()     ?: '',
                family  : row.FAMILY?.trim()    ?: '',
                genus   : row.GENUS?.trim()     ?: '',
                species : row.SPECIES?.trim()   ?: '',
                strain  : row.STRAIN?.trim()    ?: '',
                asmid   : asmid,
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
