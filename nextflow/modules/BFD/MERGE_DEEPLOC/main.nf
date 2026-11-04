include { tablesDir } from '../../common/utils.nf'

// One row per protein (DeepLoc predicts every protein), so this table is large
// at full BFD scale. bin/merge_deeploc.py reads the gz CSVs with DuckDB and
// writes Parquet directly, with no CSV intermediate.
process MERGE_DEEPLOC {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(results)

    output:
        path("deeploc.parquet"), emit: parquet

    script:
    // DuckDB sizes itself from node RAM, not the SLURM limit. Give it ~75% of
    // the task memory; it spills the ORDER BY to duckdb_tmp/ in the work dir.
    def mem_gb = Math.max(1, (task.memory.toGiga() * 0.75) as int)
    """
    module load duckdb 2>/dev/null || true
    python3 ${projectDir}/bin/merge_deeploc.py -o deeploc.parquet \\
        --memory-limit ${mem_gb}GB --threads ${task.cpus} ${results}
    rm -rf duckdb_tmp
    """

    stub:
    """
    module load duckdb 2>/dev/null || true
    printf 'Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome,Peripheral,Transmembrane,Lipid anchor,Soluble\\n' | gzip > STUB.deeploc_fast.csv.gz
    python3 ${projectDir}/bin/merge_deeploc.py -o deeploc.parquet STUB.deeploc_fast.csv.gz
    """
}
