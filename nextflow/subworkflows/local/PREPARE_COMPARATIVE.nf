include { STAGE_FILES }   from '../../modules/comparative/stage_files/STAGE_FILES/main.nf'
include { makeSampleTag } from '../../modules/common/utils.nf'

workflow PREPARE_COMPARATIVE {
    main:
    def group_map = [:]
    if (params.group) {
        def gf = file(params.group)
        def first = true
        gf.eachLine { line ->
            if (first) { first = false; return }
            def cols = line.trim().split(',', 2)
            if (cols.size() == 2 && cols[0]) group_map[cols[0].trim()] = cols[1].trim()
        }
    }

    def ignored = [] as Set
    if (params.ignore) {
        file(params.ignore).eachLine { line ->
            def lt = line.trim()
            if (lt) ignored << lt
        }
    }

    def taxon_filters = []
    if (params.taxon) {
        params.taxon.split(',').each { spec ->
            def kv = spec.trim().split(':', 2)
            if (kv.size() == 2) taxon_filters << [rank: kv[0].trim().toUpperCase(), value: kv[1].trim()]
        }
    }

    def species_ch = channel.fromPath(params.samples)
        .splitCsv(header: true, strip: true)
        .filter { row ->
            def lt = row['LOCUSTAG']
            if (!lt || ignored.contains(lt)) return false
            def in_group = group_map.containsKey(lt)
            def in_taxon = taxon_filters.any { f -> row[f.rank]?.equalsIgnoreCase(f.value) }
            return in_group || in_taxon
        }
        .map { row ->
            def lt  = row['LOCUSTAG']
            def grp = group_map.containsKey(lt)
                        ? group_map[lt]
                        : (taxon_filters.find { f -> row[f.rank]?.equalsIgnoreCase(f.value) }?.value ?: 'default')
            def base = makeSampleTag(row['SPECIES'] ?: '', row['STRAIN'] ?: '')
            [lt, grp, base]
        }

    def manifest_ch = species_ch
        .collectFile(
            name:    "${params.project}.manifest.csv",
            seed:    "LOCUSTAG,GROUP,BASENAME\n",
            newLine: true
        ) { lt, grp, base -> "${lt},${grp},${base}" }

    STAGE_FILES(manifest_ch, params.pep_dir, params.cds_dir)

    emit:
    proteins_dir = STAGE_FILES.out.proteins_dir
    cds_dir      = STAGE_FILES.out.cds_dir
    manifest     = STAGE_FILES.out.manifest
}
