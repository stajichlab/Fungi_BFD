include { PARALOGOSCOPE } from '../subworkflows/local/PARALOGOSCOPE.nf'
include { cleanStrain; makeSampleTag; dirIndex } from '../modules/common/utils.nf'

//
// paralogoscope — per-species whole-genome duplication dating via wgd
// (dmd -> ksd [+ syn]). Inputs follow the BFD structure:
//   input/pep/{meta.id}.proteins.fa      (not used directly by wgd dmd)
//   input/cds/{meta.id}.cds-transcripts.fa   (wgd dmd source + ksd source)
//   input/gff3/{meta.id}.gff3                   (wgd syn, when enabled)
// meta.id = makeSampleTag(SPECIES, STRAIN), the same primary key BFD.nf uses.
//
// Selection mirrors PREPARE_COMPARATIVE: --group / --taxon / --ignore,
// then --n_test limits the first N rows.
//
workflow PARALOGOSCOPE_RUN {
    main:
    // ── Selection: group / taxon / ignore (OR + exclusion, as comparative) ──
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
        if (!taxon_filters) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon CLASS:Dothideomycetes"
        }
    }

    // ── Resolve per-species input files (indexed dirs, single listing) ─────
    def cdsIndex = dirIndex(params.cds_dir)
    def gffIndex = (params.run_wgd_syn.toBoolean() && params.gff_dir && file(params.gff_dir).exists())
                        ? dirIndex(params.gff_dir) : [:]

    def species_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, strip: true)
        .filter { row ->
            def lt = row['LOCUSTAG']?.replaceAll(/[\r\n]/, '')?.trim()
            if (!lt || ignored.contains(lt)) return false
            if (!taxon_filters) return true
            return taxon_filters.any { f -> row[f.rank]?.equalsIgnoreCase(f.value) }
        }
        .map { row ->
            def lt  = row['LOCUSTAG']?.replaceAll(/[\r\n]/, '')?.trim()
            def meta = [
                id      : makeSampleTag(row['SPECIES']?.trim() ?: '', row['STRAIN']?.trim() ?: ''),
                locustag: lt,
                species : row['SPECIES']?.trim() ?: '',
                strain  : cleanStrain(row['STRAIN']?.trim() ?: ''),
            ]
            def cds = cdsIndex["${meta.id}.cds-transcripts.fa"]
            if (!cds) {
                log.warn "paralogoscope: skipping ${meta.id} (${lt}): no ${meta.id}.cds-transcripts.fa in ${params.cds_dir}"
                return null
            }
            return tuple(meta, cds)
        }
        .filter { it != null }
        .take((params.n_test as int) > 0 ? (params.n_test as int) : -1)

    def gff_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, strip: true)
        .filter { row ->
            if (!params.run_wgd_syn.toBoolean()) return false
            def lt = row['LOCUSTAG']?.replaceAll(/[\r\n]/, '')?.trim()
            return lt && !ignored.contains(lt)
        }
        .map { row ->
            def lt  = row['LOCUSTAG']?.replaceAll(/[\r\n]/, '')?.trim()
            // identical shape to species_ch's meta so .join() keys match
            def meta = [
                id      : makeSampleTag(row['SPECIES']?.trim() ?: '', row['STRAIN']?.trim() ?: ''),
                locustag: lt,
                species : row['SPECIES']?.trim() ?: '',
                strain  : cleanStrain(row['STRAIN']?.trim() ?: ''),
            ]
            def gff = gffIndex["${meta.id}.gff3"]
            if (!gff) {
                log.warn "paralogoscope: no GFF3 for ${meta.id}, wgd syn will skip it"
                return null
            }
            return tuple(meta, gff)
        }
        .filter { it != null }

    PARALOGOSCOPE(species_ch, gff_ch)

    emit:
    families = PARALOGOSCOPE.out.families
    ks       = PARALOGOSCOPE.out.ks
}
