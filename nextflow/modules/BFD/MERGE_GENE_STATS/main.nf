include { tablesDir } from '../../common/utils.nf'

process MERGE_GENE_STATS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_info.csv.gz",        emit: gene_info
    path "gene_exons.csv.gz",       emit: gene_exons
    path "gene_CDS.csv.gz",         emit: gene_CDS
    path "gene_introns.csv.gz",     emit: gene_introns
    path "gene_transcripts.csv.gz", emit: gene_transcripts
    path "gene_trnas.csv.gz",       emit: gene_trnas
    path "gene_proteins.csv.gz",    emit: gene_proteins

    script:
    """
    set -euo pipefail
    for type in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        first=1
        while IFS=\$'\\t' read -r f _mtime _size; do
            case "\$f" in *.\${type}.csv.gz) ;; *) continue ;; esac
            if [ "\$first" = "1" ]; then zcat "\$f"; first=0
            else zcat "\$f" | tail -n +2; fi
        done < ${manifest} | pigz -p ${task.cpus} > \${type}.csv.gz
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > \${f}.csv.gz
    done
    """
}
