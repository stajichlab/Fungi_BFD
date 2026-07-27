//
// PICK_REPRESENTATIVE_STRAIN — select one representative strain per species from
// ANI data, write abinitio_reuse_assignments.csv, and backfill the shared
// ab-initio parameter store from already-predicted representatives.
//
// Input channels:
//   ani_tsv         — merged all-pairs TSV (from CONCAT_ANI_TSVS or /dev/null)
//   predict_input   — predict_input written to TSV for strain metadata
//   samples_csv     — samples.csv for phylum/order/family column lookup
//   busco_dir       — results/genome_stats/BUSCO_genome/ for representative scoring
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

    publishDir "${params.target}/_reuse_assignments", mode: 'copy'

    input:
        path ani_tsv
        path predict_input_tsv
        path samples_csv
        path busco_dir

    output:
        path("abinitio_reuse_assignments.csv"),  emit: outCSV
        path("repr_assignments.tsv"),            emit: outAssignments

    script:
    def threshold   = params.ani_reuse_threshold ?: 99.0
    def shared_root = params.gene_prediction_shared_abinitio ?: ''
    def target      = params.target ?: ''
    def aug_cfg     = params.augustus_config ?: ''
    def out_dir     = "${params.target}/_reuse_assignments"
    """
    pick_representative_strain.py \
        --ani-tsv "${ani_tsv}" \
        --predict-input "${predict_input_tsv}" \
        --samples "${samples_csv}" \
        --busco-dir "${busco_dir}" \
        --out-dir "${out_dir}" \
        --ani-threshold "${threshold}" \
        ${shared_root ? "--shared-root ${shared_root}" : ''} \
        ${target ? "--target ${target}" : ''} \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    """

    stub:
    """
    printf 'species\\tout\\tis_representative\\trepresentative_out\\tani_to_representative\\treuse_eligible\\n' > abinitio_reuse_assignments.csv
    printf '' > repr_assignments.tsv
    """
}
