process CALC_GENE_STATS {
    label    'genestats'
    tag      "${meta.locustag}"
    storeDir "${params.genome_stats_outdir}/gene_stats"

    input:
    tuple val(meta), path(gff_file), path(dna_file)

    output:
    path "${meta.id}.gene_info.csv.gz",        emit: gene_info
    path "${meta.id}.gene_exons.csv.gz",       emit: gene_exons
    path "${meta.id}.gene_CDS.csv.gz",         emit: gene_CDS
    path "${meta.id}.gene_introns.csv.gz",     emit: gene_introns
    path "${meta.id}.gene_transcripts.csv.gz", emit: gene_transcripts
    path "${meta.id}.gene_trnas.csv.gz",       emit: gene_trnas
    path "${meta.id}.gene_proteins.csv.gz",    emit: gene_proteins

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
        mv \${f}.csv.gz ${meta.id}.\${f}.csv.gz
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > ${meta.id}.\${f}.csv.gz
    done
    """
}
