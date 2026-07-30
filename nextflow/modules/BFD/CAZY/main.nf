process RUN_CAZY {
    tag        "${meta.locustag}"
    label      'cazy'
    storeDir   { "${params.outdir}/cazy/${meta.id}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.overview.tsv.gz"),   emit: overview
        path("${meta.id}.cazymes.tsv.gz"),    emit: cazymes
        path("${meta.id}.substrates.tsv.gz"), emit: substrates

    script:
    """
    module load dbcanlight
    mkdir -p ${meta.id}
    dbcanlight search -i ${proteins} -m cazyme -o ${meta.id} -t ${task.cpus}
    dbcanlight search -i ${proteins} -m sub    -o ${meta.id} -t ${task.cpus}
    dbcanlight conclude ${meta.id}
    pigz -f ${meta.id}/cazymes.tsv ${meta.id}/substrates.tsv ${meta.id}/overview.tsv
    mv ${meta.id}/overview.tsv.gz   ${meta.id}.overview.tsv.gz
    mv ${meta.id}/cazymes.tsv.gz    ${meta.id}.cazymes.tsv.gz
    mv ${meta.id}/substrates.tsv.gz ${meta.id}.substrates.tsv.gz
    """

    stub:
    """
    printf 'Gene_ID\\tEC\\tcazyme_fam\\tsub_fam\\tdiamond_fam\\tSubstrate\\t#ofTools\\n' | gzip > ${meta.id}.overview.tsv.gz
    printf 'HMM_Profile\\tProfile_Length\\tGene_ID\\tGene_Length\\tEvalue\\tProfile_Start\\tProfile_End\\tGene_Start\\tGene_End\\tCoverage\\n' | gzip > ${meta.id}.cazymes.tsv.gz
    printf 'Gene_ID\\tSubstrate\\n' | gzip > ${meta.id}.substrates.tsv.gz
    """
}
