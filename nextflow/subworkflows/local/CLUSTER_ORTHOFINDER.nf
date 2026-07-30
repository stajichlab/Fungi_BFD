include { ORTHOFINDER_RUN }   from '../../modules/comparative/orthofinder/ORTHOFINDER_RUN/main.nf'
include { ORTHOFINDER_PARSE } from '../../modules/comparative/orthofinder/ORTHOFINDER_PARSE/main.nf'

workflow CLUSTER_ORTHOFINDER {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    ORTHOFINDER_RUN(proteins_dir)
    ORTHOFINDER_PARSE(ORTHOFINDER_RUN.out.out_dir)

    emit:
    orthogroups = ORTHOFINDER_PARSE.out.orthogroups
}
