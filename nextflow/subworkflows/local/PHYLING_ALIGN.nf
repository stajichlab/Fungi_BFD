include { MARKERSET_PREPARE } from '../../modules/phyling/markerset/MARKERSET_PREPARE/main.nf'
include { PHYLING_ALIGN }     from '../../modules/phyling/align/PHYLING_ALIGN/main.nf'
include { PHYLING_FILTER }    from '../../modules/phyling/filter/PHYLING_FILTER/main.nf'
include { PHYLING_TREE }      from '../../modules/phyling/tree/PHYLING_TREE/main.nf'
include { cleanStrain; makeSampleTag } from '../../modules/common/utils.nf'

workflow PHYling {
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
    }

    def taxon_slug = params.taxon
        ? (params.taxon as String).split(':', 2)[1].replaceAll(/\s+/, '_')
        : 'all'

    log.info "Output directory: ${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/"

    def input_dir = params.seq_type == 'cds' ? params.cds_dir : params.pep_dir
    def suffix    = params.seq_type == 'cds' ? '.cds-transcripts.fa' : '.proteins.fa'

    def taxonFilter
    if (params.taxon) {
        def parts2 = (params.taxon as String).split(':', 2)
        def taxRank  = parts2[0].toUpperCase()
        def taxValue = parts2[1]
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { row -> true }
    }

    def fasta_ch = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = cleanStrain(row.STRAIN?.trim() ?: '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def fasta    = file("${input_dir}/${basename}${suffix}", glob: false)
            if (!fasta.exists()) {
                log.warn "PHYling: skipping ${basename} (${locustag}): ${basename}${suffix} not found"
                return null
            }
            return fasta
        }
        .filter { it != null }
        .take(params.n_test > 0 ? params.n_test as int : -1)
        .collect()

    def markerset_ch = channel.fromList(
        params.markerset.tokenize(',').collect { it.trim() }
    )

    MARKERSET_PREPARE(markerset_ch)

    def align_input = MARKERSET_PREPARE.out.ready
        .map { markerset, markerset_dir -> tuple(markerset, taxon_slug, markerset_dir) }
        .combine(fasta_ch)

    PHYLING_ALIGN(align_input)
    PHYLING_FILTER(PHYLING_ALIGN.out.align)
    PHYLING_TREE(PHYLING_FILTER.out.filtered)
}
