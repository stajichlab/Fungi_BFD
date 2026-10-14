// DRAFT. Gather step: merge the chunk GTFs of one genome in chunk order, renumber gene and
// transcript IDs, validate against the whole genome, and publish the genome's GeneMark result
// to genemark_store/<out>/ -- the only place the chunked route publishes. Publishing happens
// after validation, so a half-merged result never reaches the store.
// Outputs use the same names and keyed-tuple shapes as the monolithic route, so callers join
// on them unchanged. `fresh` is true when the model was just trained (publish it as the new
// .genemark.mod); false for a stored shared model (fast-reuse emits no new .mod, as before).
// The production shell is `/bin/bash -l` (no `set -e`), so strict mode is set below.
process GENEMARK_MERGE {
    tag "$out"
    publishDir path: { "${params.genemark_store}/${out}" }, mode: 'copy', overwrite: true,
        enabled: params.genemark_store ? true : false

    cpus   2
    memory '4 GB'
    time   '1h'

    input:
    tuple val(out), val(n_chunks), val(species), val(genome_fa), val(fresh), path(chunk_gtfs, stageAs: 'gtfs/*'), path(model), path(other_pending)

    output:
    tuple val(out), path("${out}.genemark.gtf"), emit: gtf
    tuple val(out), path("${out}.other.gff3"), emit: other_gff
    tuple val(out), val(species), path("${out}.genemark.mod"), optional: true, emit: mod

    script:
    """
    set -euo pipefail
    # 1. exactly the expected chunks, numbered 0..n-1 with no gap and no duplicate
    EXPECTED=${n_chunks}
    GOT=\$(ls gtfs/*.gtf | wc -l)
    [ "\$GOT" -eq "\$EXPECTED" ] || { echo "ERROR: GENEMARK_MERGE ${out}: got \$GOT chunk GTFs, expected \$EXPECTED" >&2; exit 1; }
    for i in \$(seq 0 \$((EXPECTED - 1))); do
        [ -e "\$(printf 'gtfs/${out}.chunk_%04d.gtf' \$i)" ] || { echo "ERROR: GENEMARK_MERGE ${out}: chunk \$i is missing" >&2; exit 1; }
    done
    GTFS=\$(for i in \$(seq 0 \$((EXPECTED - 1))); do printf 'gtfs/${out}.chunk_%04d.gtf ' \$i; done)

    python "${workflow.projectDir}/bin/genemark_merge_gtf.py" --out "${out}.genemark.gtf" \$GTFS
    [ -s "${out}.genemark.gtf" ] || { echo "ERROR: GENEMARK_MERGE ${out}: merged GTF is empty" >&2; exit 1; }

    # 2. every GTF contig name must be a contig of the whole genome; an empty contig list is a failure
    case "${genome_fa}" in
        *.gz) pigz -dc "${genome_fa}" > genome_for_names.fa ;;
        *)    cp "${genome_fa}" genome_for_names.fa ;;
    esac
    awk '/^>/ { print substr(\$1, 2) }' genome_for_names.fa > contig_names.txt
    [ -s contig_names.txt ] || { echo "ERROR: GENEMARK_MERGE ${out}: no contig names read from ${genome_fa}" >&2; exit 1; }
    if ! awk 'NR == FNR { ids[\$1] = 1; next }
              !/^#/ && NF >= 9 && !(\$1 in ids) { print "unknown sequence " \$1 > "/dev/stderr"; bad = 1; exit 1 }
              END { exit bad }' contig_names.txt "${out}.genemark.gtf"; then
        echo "ERROR: GENEMARK_MERGE ${out}: merged GTF names sequences that are not in the genome" >&2
        exit 1
    fi
    rm -f genome_for_names.fa contig_names.txt

    cp "${other_pending}" "${out}.other.gff3"
    if [ "${fresh}" = "true" ]; then
        [ -s "${model}" ] || { echo "ERROR: GENEMARK_MERGE ${out}: fresh model is empty" >&2; exit 1; }
        cp "${model}" "${out}.genemark.mod"
        # same predict_misc persistence as the monolithic route (see GENEMARK_TRAIN)
        LOWER_OUT=\$(printf '%s' "${out}" | tr '[:upper:]' '[:lower:]')
        AB_INITIO_DIR="${params.target}/${out}/predict_misc/ab_initio_parameters"
        mkdir -p "\$AB_INITIO_DIR"
        cp "${out}.genemark.mod" "\$AB_INITIO_DIR/\${LOWER_OUT}.genemark.mod"
    fi
    echo "[INFO] GENEMARK_MERGE complete for ${out}"
    """

    stub:
    """
    touch "${out}.genemark.gtf" "${out}.other.gff3"
    if [ "${fresh}" = "true" ]; then touch "${out}.genemark.mod"; fi
    """
}
