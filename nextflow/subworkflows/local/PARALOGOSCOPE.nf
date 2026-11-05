include { WGD_DMD } from '../../modules/comparative/wgd/WGD_DMD/main.nf'
include { WGD_KSD } from '../../modules/comparative/wgd/WGD_KSD/main.nf'
include { WGD_PAIR_KS } from '../../modules/comparative/wgd/WGD_PAIR_KS/main.nf'
include { WGD_SYN } from '../../modules/comparative/wgd/WGD_SYN/main.nf'
include { MERGE_WGD_KSD } from '../../modules/comparative/wgd/MERGE_WGD_KSD/main.nf'

workflow PARALOGOSCOPE {
    take:
    species_ch   // channel: tuple(meta, cds) per species
    gff_ch       // channel: tuple(meta, gff) per species (may be empty)

    main:
    WGD_DMD(species_ch)
    def dmd_fams = WGD_DMD.out.families

    // WGD_DMD echoes (meta, cds, families); strip cds for the next step.
    def fam_only = dmd_fams
        .map { meta, cds, families -> tuple(meta, families) }

    def fam_cds = fam_only.join(species_ch)
        .map { meta, families, cds -> tuple(meta, families, cds) }

    def ks_ch
    // params.wgd_ks_method: 'pair' (WGD_PAIR_KS, pair-level codeml; default), 'ksd_pairwise'
    // (wgd ksd --pairwise) or 'ksd' (wgd ksd default, family-level gap stripping).
    if (!(params.wgd_ks_method in ['pair', 'ksd_pairwise', 'ksd'])) {
        error "--wgd_ks_method must be pair, ksd_pairwise or ksd (got '${params.wgd_ks_method}')"
    }
    if (params.wgd_ks_method == 'pair') {
        // bundle params.wgd_pair_bundle genomes per task (about 1-1.5 h per task)
        def bundled = fam_cds
            .buffer(size: (params.wgd_pair_bundle as int), remainder: true)
            .map { lst -> tuple(lst.collect { it[0] }, lst.collect { it[1] }, lst.collect { it[2] }) }
        WGD_PAIR_KS(bundled)
        ks_ch = WGD_PAIR_KS.out.ks.flatten()
    } else {
        WGD_KSD(fam_cds)
        ks_ch = WGD_KSD.out.ks
    }

    if (params.run_wgd_syn.toBoolean()) {
        WGD_SYN(
            fam_only.join(gff_ch)
                .map { meta, families, gff -> tuple(meta, families, gff) }
        )
    }

    // Final merge: after all of this invocation's ksds finish (or immediately
    // when all are cached), glob every ks.tsv published so far off disk into
    // tables/wgd.ks.parquet. outdir must be absolute (task workdirs don't share
    // the launch dir's relative path).
    def ks_outdir = params.outdir.startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"
    MERGE_WGD_KSD(ks_outdir, ks_ch.collect().map { it }, species_ch.count())

    MERGE_WGD_KSD.out.check.subscribe { f ->
        def t = f.text.trim()
        if (t != 'OK') log.error "PARALOGOSCOPE merge check found problems:\n${t}"
    }

    

    emit:
    families = WGD_DMD.out.families
    ks       = ks_ch
    anchors  = (params.run_wgd_syn.toBoolean()) ? WGD_SYN.out.anchors : Channel.empty()
    ks_table = MERGE_WGD_KSD.out.ks_table
}
