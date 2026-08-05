//
// PICK_REPRESENTATIVE_STRAIN — select one representative strain per species from
// ANI data, write abinitio_reuse_assignments.csv, and backfill the shared
// ab-initio parameter store from already-predicted representatives.
//
// Input channels:
//   ani_tsv         — merged all-pairs TSV (from CONCAT_ANI_TSVS or /dev/null)
//   predict_input   — predict_input written to TSV for strain metadata
//   samples_csv     — samples.csv for phylum/order/family column lookup
//   tables_dir      — tables/ (busco_genome.parquet + asm_stats.parquet) for
//                     representative scoring (BUSCO_GENOME completeness + N50)
//   reuse_out_dir   — ${params.target}/_reuse_assignments/ for CSV output
//
// Output:
//   outCSV          — abinitio_reuse_assignments.csv (for FUNANNOTATE_PREDICTION)
//   outAssignments  — repr_assignments.tsv (out, asmid, species, is_representative,
//                     representative_out, ani_to_representative, reuse_eligible)
//
// Logic (mirrors nextflow/bin/species_reuse_clusters.py):
//   - group by species; only species with >=2 strains are candidates
//   - pick_representative: highest BUSCO completeness, N50 tiebreak, alphabetical
//     out final tiebreak, but ONLY among candidates with ANI coverage in ani_tsv
//     (prevents a top-BUSCO strain added after ANI ran from orphaning every other
//     strain's reuse eligibility)
//   - each non-representative strain: reuse_eligible if ANI to rep >= threshold
//   - fail-closed: no ANI row → not eligible (sparse ANI methods omit low-identity)
//   - write abinitio_reuse_assignments.csv
//   - backfill shared params if representative has predict_misc/ab_initio_parameters/
//

process PICK_REPRESENTATIVE_STRAIN {
    tag   "PICK_REPRESENTATIVE_STRAIN"
    label 'report'

    publishDir "${params.target}/_reuse_assignments", mode: 'copy', pattern: 'abinitio_reuse_assignments*.csv'
    publishDir "${params.target}/_reuse_assignments", mode: 'copy', pattern: 'repr_assignments.tsv'

    input:
        path ani_tsv
        path predict_input_tsv
        path samples_csv
        path tables_dir

    output:
        path("abinitio_reuse_assignments.csv"),  emit: outCSV
        path("repr_assignments.tsv"),            emit: outAssignments

    script:
    def threshold   = params.ani_reuse_threshold ?: 99.0
    def shared_root = params.gene_prediction_shared_abinitio ?: ''
    def target      = params.target ?: ''
    def aug_cfg     = params.augustus_config ?: ''
    """
    mkdir -p _reuse_assignments
    # duckdb is not in the python:3.12 base image (used via the 'report' label) —
    # install it first, same pattern as COMBINE_ANI_TABLE.
    pip install --quiet --target /tmp/python_packages --ignore-installed duckdb
    PYTHONPATH="/tmp/python_packages:\${PYTHONPATH}" \
    python "${projectDir}/bin/pick_representative_strain.py" \
        --ani-tsv "${ani_tsv}" \
        --predict-input "${predict_input_tsv}" \
        --samples "${samples_csv}" \
        --tables-dir "${tables_dir}" \
        --out-dir "_reuse_assignments" \
        --ani-threshold "${threshold}" \
        ${shared_root ? "--shared-root ${shared_root}" : ''} \
        ${target ? "--target ${target}" : ''} \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    mv _reuse_assignments/abinitio_reuse_assignments*.csv .
    mv _reuse_assignments/repr_assignments.tsv .
    """

    stub:
    """
    printf 'species\\tout\\tis_representative\\trepresentative_out\\tani_to_representative\\treuse_eligible\\n' > abinitio_reuse_assignments.csv
    printf '' > repr_assignments.tsv
    """
}
