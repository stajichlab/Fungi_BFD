process BATCH_CODON_FREQ {
    label      'genestats'
    publishDir "${params.genome_stats_outdir}/codon_freq", mode: 'copy'

    input:
    tuple val(locustags), val(basenames), path(prots)

    output:
    path "*.codon_freq.csv.gz", emit: csv

    script:
    def baseList = basenames instanceof List ? basenames : [basenames]
    def protList = prots     instanceof List ? prots     : [prots]
    def cmds = [baseList, protList].transpose().collect { basename, prot ->
        "python3 ${params.scripts}/calculate_codon_freq.py ${prot.name} -o ${basename}.codon_freq.csv.gz"
    }.join('\n')
    """
    module load biopython
    ${cmds}
    """

    stub:
    def baseList = basenames instanceof List ? basenames : [basenames]
    """
    ${ baseList.collect { "printf 'species_prefix,codon,frequency\\n' | gzip > ${it}.codon_freq.csv.gz" }.join('\n') }
    """
}
