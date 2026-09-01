include { WGD_DMD } from '../../modules/comparative/wgd/WGD_DMD/main.nf'
include { WGD_KSD } from '../../modules/comparative/wgd/WGD_KSD/main.nf'
include { WGD_SYN } from '../../modules/comparative/wgd/WGD_SYN/main.nf'
include { MERGE_WGD_KSD } from '../../modules/comparative/wgd/MERGE_WGD_KSD/main.nf'

workflow PARALOGOSCOPE {
    take:
    species_ch   // channel: tuple(meta, cds) per species
    gff_ch       // channel: tuple(meta, gff) per species (may be empty)

    main:
    WGD_DMD(species_ch)

    // WGD_DMD echoes (meta, cds, families); strip cds for the next step.
    def fam_only = WGD_DMD.out.families
        .map { meta, cds, families -> tuple(meta, families) }

    WGD_KSD(
        fam_only.join(species_ch)
            .map { meta, families, cds -> tuple(meta, families, cds) }
    )

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
    MERGE_WGD_KSD(ks_outdir, WGD_KSD.out.ks.collect().map { it })

    emit:
    families = WGD_DMD.out.families
    ks       = WGD_KSD.out.ks
    anchors  = (params.run_wgd_syn.toBoolean()) ? WGD_SYN.out.anchors : Channel.empty()
    ks_table = MERGE_WGD_KSD.out.ks_table
}
