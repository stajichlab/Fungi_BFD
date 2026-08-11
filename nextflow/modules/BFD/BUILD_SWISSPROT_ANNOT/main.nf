include { tablesDir } from '../../common/utils.nf'

process BUILD_SWISSPROT_ANNOT {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(dat)

    output:
        path("swissprot_annot.parquet"), emit: annot_parquet

    // One-shot table build: parse the UniProtKB/Swiss-Prot flatfile into a
    // per-accession annotation table (protein name, gene, organism, GO terms,
    // EC numbers, InterPro and Pfam cross-refs). Keyed by primary accession;
    // join to swissprot.parquet on swissprot_acc. Rebuilt whenever this
    // process runs (no storeDir) so an updated flatfile always propagates.
    script:
    """
    python3 ${projectDir}/bin/parse_uniprot_sprot.py ${dat} -o swissprot_annot.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot_annot.csv', sample_size=-1)) TO 'swissprot_annot.parquet' (FORMAT PARQUET);"
    rm -f swissprot_annot.csv
    """

    stub:
    """
    printf 'accession,entry_name,protein_name,gene_name,gene_synonyms,organism,go_terms,ec_numbers,interpro,pfam\\n' > swissprot_annot.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot_annot.csv', sample_size=-1)) TO 'swissprot_annot.parquet' (FORMAT PARQUET);"
    rm -f swissprot_annot.csv
    """
}
