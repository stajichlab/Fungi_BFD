process RUN_CAZY {
    tag        "${locustag}"
    label      'cazy'
    storeDir   { "${params.outdir}/cazy/${basename}" }

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.overview.tsv.gz"),   emit: overview
        path("${basename}.cazymes.tsv.gz"),    emit: cazymes
        path("${basename}.substrates.tsv.gz"), emit: substrates

    script:
    """
    module load dbcanlight
    mkdir -p ${basename}
    dbcanlight search -i ${proteins} -m cazyme -o ${basename} -t ${task.cpus}
    dbcanlight search -i ${proteins} -m sub    -o ${basename} -t ${task.cpus}
    dbcanlight conclude ${basename}
    pigz -f ${basename}/cazymes.tsv ${basename}/substrates.tsv ${basename}/overview.tsv
    mv ${basename}/overview.tsv.gz   ${basename}.overview.tsv.gz
    mv ${basename}/cazymes.tsv.gz    ${basename}.cazymes.tsv.gz
    mv ${basename}/substrates.tsv.gz ${basename}.substrates.tsv.gz
    """

    stub:
    """
    printf 'Gene_ID\\tEC\\tcazyme_fam\\tsub_fam\\tdiamond_fam\\tSubstrate\\t#ofTools\\n' | gzip > ${basename}.overview.tsv.gz
    printf 'HMM_Profile\\tProfile_Length\\tGene_ID\\tGene_Length\\tEvalue\\tProfile_Start\\tProfile_End\\tGene_Start\\tGene_End\\tCoverage\\n' | gzip > ${basename}.cazymes.tsv.gz
    printf 'Gene_ID\\tSubstrate\\n' | gzip > ${basename}.substrates.tsv.gz
    """
}
