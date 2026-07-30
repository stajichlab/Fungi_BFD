process DELIVER_MASK {
    tag    "${asmid}"
    label  'earlgrey_deliver'

    publishDir path: { workflow.stubRun ? "${launchDir}/work/stub_masked" : params.masked_dir },
               mode: 'copy', overwrite: true

    input:
        tuple val(asmid), path(masked)

    output:
        path("${asmid}.masked.fasta.gz")

    script:
    'true'

    stub:
    'true'
}
