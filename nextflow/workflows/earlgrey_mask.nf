include { SELECT_REPS }         from '../modules/earlgrey/SELECT_REPS/main.nf'
include { EARLGREY_BUILD_LIB }  from '../modules/earlgrey/EARLGREY_BUILD_LIB/main.nf'
include { REPEATMASK_STRAIN }   from '../modules/earlgrey/REPEATMASK_STRAIN/main.nf'
include { DELIVER_MASK }        from '../modules/earlgrey/DELIVER_MASK/main.nf'
include { genomeFile }          from '../modules/common/utils.nf'

workflow EARLGREY_MASK {
    def reps = SELECT_REPS(
        file(params.samples,             glob: false),
        file(params.asm_stats_parquet,   glob: false),
        file(params.busco_genome_parquet, glob: false),
    )

    def asmidFilter = params.asmid
        ? { row ->
            def target  = (params.asmid as String).trim()
            def rep     = row.REP_ASMID?.trim()
            def members = (row.MEMBER_ASMIDS ?: '').split(';').collect { it?.trim() }
            rep == target || members.contains(target)
          }
        : { row -> true }
    if (params.asmid) {
        log.info "ASMID filter: processing only the species containing '${params.asmid}'"
    }

    def records = reps
        .splitCsv(header: true)
        .filter(asmidFilter)
        .take(params.n_test > 0 ? params.n_test as int : -1)
        .map { row ->
            tuple(row.SPECIES?.trim(), row.REP_ASMID?.trim(), (row.MEMBER_ASMIDS ?: '').trim())
        }

    def build_in = records.map { species, rep_asmid, _members ->
        def g = genomeFile("${params.genome_dir}/${rep_asmid}${params.genome_suffix}")
        if (!g.exists() && !workflow.stubRun) {
            log.warn "Skipping ${species}: representative genome not found at ${g}"
            return null
        }
        tuple(species, rep_asmid, g)
    }.filter { it != null }

    EARLGREY_BUILD_LIB(build_in)

    def members_ch = records.flatMap { species, _rep, members ->
        (members ? members.split(';') : [])
            .findAll { it?.trim() }
            .collect { asm -> tuple(species, asm.trim()) }
    }

    def lib_by_species = EARLGREY_BUILD_LIB.out.lib
        .map { species, _rep, lib -> tuple(species, lib) }

    def mask_in = members_ch
        .combine(lib_by_species, by: 0)
        .map { species, asmid, library ->
            def g = genomeFile("${params.genome_dir}/${asmid}${params.genome_suffix}")
            if (!g.exists() && !workflow.stubRun) {
                log.warn "Skipping member ${asmid} (${species}): genome not found at ${g}"
                return null
            }
            tuple(asmid, species, g, library)
        }
        .filter { it != null }

    REPEATMASK_STRAIN(mask_in)

    def masked_ch = EARLGREY_BUILD_LIB.out.masked
        .mix(REPEATMASK_STRAIN.out)

    DELIVER_MASK(masked_ch)
}
