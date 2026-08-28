include { WGD_DMD } from '../../modules/comparative/wgd/WGD_DMD/main.nf'
include { WGD_KSD } from '../../modules/comparative/wgd/WGD_KSD/main.nf'
include { WGD_SYN } from '../../modules/comparative/wgd/WGD_SYN/main.nf'

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

    emit:
    families = WGD_DMD.out.families
    ks       = WGD_KSD.out.ks
    anchors  = (params.run_wgd_syn.toBoolean()) ? WGD_SYN.out.anchors : Channel.empty()
}
