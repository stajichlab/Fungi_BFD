include { tablesDir; codeMd5 } from '../../common/utils.nf'

process BUILD_SWISSPROT_ANNOT {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(dat)

    output:
        path("swissprot_annot.parquet"),   emit: annot_parquet
        path("swissprot_subject.parquet"), emit: subject

    // One-shot table build from the UniProtKB/Swiss-Prot flatfile:
    //   swissprot_annot.parquet   per-accession annotation (protein name, gene,
    //                             organism, GO, EC, InterPro, Pfam, length, taxid)
    //   swissprot_subject.parquet the same plus an integer sp_id and org_code;
    //                             swissprot_hits / swissprot_best reference sp_id
    //                             (sql/merge/swissprot_subject.sql).
    // sp_id follows accession order, so it is stable only within one release:
    // a new flatfile makes MERGE_SWISSPROT stale (swissprot_dat is one of its
    // staleness inputs in BFD_FUNCTIONAL.nf).
    script:
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    python3 ${projectDir}/bin/parse_uniprot_sprot.py ${dat} -o swissprot_annot.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot_annot.csv', sample_size=-1)) TO 'swissprot_annot.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);"
    python3 ${projectDir}/bin/duckdb_merge.py --sql ${projectDir}/sql/merge/swissprot_subject.sql \\
        --files /dev/null --out swissprot_subject.parquet --threads ${task.cpus} \\
        --param annot="'swissprot_annot.csv'"
    rm -f swissprot_annot.csv
    """

    stub:
    """
    printf 'accession,entry_name,protein_name,gene_name,gene_synonyms,organism,go_terms,ec_numbers,interpro,pfam,length,ncbi_taxid\\n' > swissprot_annot.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot_annot.csv', sample_size=-1)) TO 'swissprot_annot.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT 1::INTEGER AS sp_id, 'P00000' AS accession) TO 'swissprot_subject.parquet' (FORMAT PARQUET);"
    rm -f swissprot_annot.csv
    """
}
