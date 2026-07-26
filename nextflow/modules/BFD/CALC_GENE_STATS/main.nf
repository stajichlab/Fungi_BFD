process CALC_GENE_STATS {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/gene_stats"

    input:
    tuple val(locustag), val(basename), path(gff_file), path(dna_file)

    output:
    path "${basename}.gene_info.csv.gz",        emit: gene_info
    path "${basename}.gene_exons.csv.gz",       emit: gene_exons
    path "${basename}.gene_CDS.csv.gz",         emit: gene_CDS
    path "${basename}.gene_introns.csv.gz",     emit: gene_introns
    path "${basename}.gene_transcripts.csv.gz", emit: gene_transcripts
    path "${basename}.gene_trnas.csv.gz",       emit: gene_trnas
    path "${basename}.gene_proteins.csv.gz",    emit: gene_proteins

    script:
    """
    set -euo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load biopython
    module load bedtools/2.30.0
    python3 ${params.scripts}/build_genestats_table.py \\
        ${gff_file} \\
        -d . \\
        -o .
    pigz gene_info.csv gene_exons.csv gene_CDS.csv gene_introns.csv \\
         gene_transcripts.csv gene_trnas.csv gene_proteins.csv
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        mv \${f}.csv.gz ${basename}.\${f}.csv.gz
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > ${basename}.\${f}.csv.gz
    done
    """
}
