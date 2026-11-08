include { tablesDir } from '../../common/utils.nf'

process MERGE_GENE_STATS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_info.parquet",        emit: gene_info
    path "gene_exons.parquet",       emit: gene_exons
    path "gene_CDS.parquet",         emit: gene_CDS
    path "gene_introns.parquet",     emit: gene_introns
    path "gene_transcripts.parquet", emit: gene_transcripts
    path "gene_trnas.parquet",       emit: gene_trnas
    path "gene_proteins.parquet",    emit: gene_proteins

    script:
    // DuckDB reads the per-genome .csv.gz files listed in the manifest directly
    // (bin/manifest_to_parquet.py). The former zcat-concatenate loop +
    // read_csv_auto(sample_size=-1) was 4.5x slower on 1,020 genomes and failed
    // at ~23k genomes with "Trying to seek a compressed CSV File" (2026-11-05).
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    # Sequence text is not stored in the tables (2026-11-06: CDS_sequence 169 GB,
    # peptide 115 GB). It is served by bgzip + samtools faidx copies keyed by
    # LOCUSTAG in Fungi_BFD_runs/seqdb/{proteins,cds}/ (bin/seqdb_build.sh).
    # Intron seq (23.7 GB, dropped 2026-11-07) comes from the bgzip + faidx genome:
    # samtools faidx predict_results/<out>.scaffolds.fa.gz chrom:start-end (-i for strand -1).
    for type in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        case \$type in
            gene_transcripts) excl=CDS_sequence ;;
            gene_proteins)    excl=peptide ;;
            gene_introns)     excl=seq ;;
            *)                excl= ;;
        esac
        python3 ${projectDir}/bin/manifest_to_parquet.py --manifest ${manifest} --exclude "\$excl" \\
            --suffix .\${type}.csv.gz \\
            --out \${type}.parquet --threads ${task.cpus} \\
            --memory-limit ${(task.memory.toGiga() * 0.75) as int}GB --temp-dir \${SCRATCH:-\$PWD}
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > \${f}.csv.gz
        module load duckdb 2>/dev/null || true
        duckdb -c "COPY (SELECT * FROM read_csv_auto('\${f}.csv.gz', sample_size=-1)) TO '\${f}.parquet' (FORMAT PARQUET);"
        rm -f \${f}.csv.gz
    done
    """
}
