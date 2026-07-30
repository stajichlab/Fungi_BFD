include { PREPARE_COMPARATIVE } from '../subworkflows/local/PREPARE_COMPARATIVE.nf'
include { CLUSTER_MMSEQS2 }     from '../subworkflows/local/CLUSTER_MMSEQS2.nf'
include { CLUSTER_MCL }         from '../subworkflows/local/CLUSTER_MCL.nf'
include { CLUSTER_ORTHOFINDER } from '../subworkflows/local/CLUSTER_ORTHOFINDER.nf'

workflow COMPARATIVE {
    if (!params.project) {
        error "Required parameter --project not specified.\n" +
              "  Example: --project MyComparison"
    }
    if (!params.group && !params.taxon) {
        error "At least one of --group or --taxon must be specified to select organisms.\n" +
              "  Examples: --group group.csv\n" +
              "            --taxon CLASS:Dothideomycetes\n" +
              "            --taxon PHYLUM:Ascomycota,ORDER:Hypocreales"
    }

    PREPARE_COMPARATIVE()

    if (params.run_mmseqs2.toBoolean()) {
        CLUSTER_MMSEQS2(PREPARE_COMPARATIVE.out.proteins_dir)
    }

    if (params.run_mcl.toBoolean()) {
        CLUSTER_MCL(PREPARE_COMPARATIVE.out.proteins_dir)
    }

    if (params.run_orthofinder.toBoolean()) {
        CLUSTER_ORTHOFINDER(PREPARE_COMPARATIVE.out.proteins_dir)
    }
}
