include { tablesDir } from '../../common/utils.nf'

// Per T-014 §D.2: always builds the single, full/unscoped master samples+species
// table, regardless of --taxon. taxonRowFilter() already restricts which genomes
// are *processed* under --taxon (at the base genome channel); MERGE_SAMPLES no
// longer also restricts the table it writes, so a --taxon run no longer produces
// a separate tables/<Taxon>/ subset -- see also the retired `matched` input this
// process used to take.
process MERGE_SAMPLES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path(samples)

    output:
    path "samples.parquet", emit: samples
    path "species.parquet", emit: species

    script:
    """
    python3 ${params.scripts}/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${params.scripts}/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    """

    stub:
    """
    python3 ${params.scripts}/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${params.scripts}/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    """
}
