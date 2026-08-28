// Batched variant of GENOME_CLEAN. Receives a LIST of per-genome tuples and stages
// the FCS-GX database into /dev/shm ONCE (~30 min), then cleans every genome in the
// batch sequentially against that in-memory DB. This amortizes the expensive staging
// step over ~clean_batch_size genomes instead of paying it per genome.
//
// Outputs are written directly to ${launchDir}/input_clean_genomes/<asmid>.fa (the same
// location GENOME_CLEAN's storeDir uses) and a per-batch manifest listing every cleaned
// assembly. Genomes whose .fa already exists are skipped, so a killed/retried batch
// resumes without redoing finished assemblies.
process GENOME_CLEAN_BATCH {
    tag "clean_batch_${task.index}"

    cpus   16
    memory '450 GB'
    time   '7d'

    input:
    tuple val(items), val(taxondb)
    path samples_csv

    output:
    path "clean_batch_*.manifest.tsv", emit: manifest
    path "TO_ADD_TO_SUPRESS.csv", emit: suppress

    script:
    def batch_tsv = items.collect { row -> "${row[1]}\t${row[8]}\t${row[9]}" }.join('\n')
    """
    set -uo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load apptainer
    AAFTF_SIF=${params.aaftf_sif}
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    # AAFTF + taxonkit now come from the AAFTF SIF (replaces `module load AAFTF`
    # and `module load taxonkit`). fcs_gx_purge takes its DB via an explicit
    # --db /dev/shm/gxdb/all flag (not the image's AAFTF_DB env default), so
    # only /dev/shm needs binding for that -- no separate AAFTF_DB bind.
    # \$SCRATCH (node-local, NOT under /bigdata) also needs an explicit bind --
    # fcs_gx_purge reads/writes \$SCRATCH/\${asmid}.*.fa[sta] directly per-genome
    # in the loop below, and manually-built \$SING commands get none of
    # Nextflow's automatic task-workdir binding (same missing-bind bug class as
    # GENEMARK_RUN, confirmed 2026-08-24).
    SING_BINDS="--bind \${SCRATCH}:\${SCRATCH},${taxondb}:${taxondb},/dev/shm:/dev/shm"
    SING="apptainer exec \${SING_BINDS} \${AAFTF_SIF}"

    TAXONKIT_DB=${taxondb}
    DEST=${launchDir}/input_clean_genomes
    mkdir -p \$DEST/clean

    MANIFEST=clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST

    SUPPRESS=TO_ADD_TO_SUPRESS.csv
    : > \$SUPPRESS

    # Appends "asmid,Assembly too short <N> bp in length" to \$SUPPRESS when the
    # cleaned assembly at \$1 (gzipped or plain fasta) totals fewer than
    # params.min_assembly_len bp -- almost certainly junk/contamination-only or a
    # failed FCS-GX purge, not a usable genome.
    check_min_length() {
        local id="\$1" fa="\$2"
        local total_len
        if [[ "\$fa" == *.gz ]]; then
            total_len=\$(pigz -dc "\$fa" | grep -v '^>' | tr -d '\\n\\r' | wc -c)
        else
            total_len=\$(grep -v '^>' "\$fa" | tr -d '\\n\\r' | wc -c)
        fi
        if [ "\$total_len" -lt "${params.min_assembly_len}" ]; then
            echo "[WARN] \$id cleaned assembly is only \$total_len bp (< ${params.min_assembly_len}); flagging for suppression" >&2
            printf '%s,Assembly too short %s bp in length\\n' "\$id" "\$total_len" >> \$SUPPRESS
        fi
    }

    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF

    n_total=\$(grep -c . batch.tsv || true)
    echo "[INFO] batch ${task.index}: \$n_total genomes to consider"

    # Stage the FCS-GX DB into /dev/shm ONCE for the whole batch (~30 min). Keep it
    # across the loop; remove the RAM copy ourselves when the batch finishes.
    export FCS_GX_KEEP_SHM=1
    # Persist per-task FCS-GX /dev/shm sync timing outside work/ (cleanup=true).
    export FCS_GX_TIMING_LOG="${launchDir}/logs/nextflow/fcs_gx_shm_timing.tsv"
    source ${projectDir}/bin/setup_fcs_shm.sh
    trap 'rm -rf /dev/shm/gxdb 2>/dev/null || true' EXIT
    if [ ! -f /dev/shm/gxdb/all.gxi ]; then
        echo "[ERROR] FCS-GX DB not staged into /dev/shm/gxdb; aborting batch" >&2
        exit 1
    fi

    i=0
    while IFS=\$'\\t' read -r asmid gz taxonid; do
        [ -z "\$asmid" ] && continue
        i=\$((i+1))
        target=\$DEST/\${asmid}.fa.gz
        # Accept a prior uncompressed .fa as already-cleaned too (back-compat).
        if [ -s "\$target" ]; then
            echo "[\$i/\$n_total][SKIP] \$asmid already cleaned"
            check_min_length "\$asmid" "\$target"
            printf '%s\\t%s\\n' "\$asmid" "\$target" >> \$MANIFEST
            continue
        elif [ -s "\$DEST/\${asmid}.fa" ]; then
            echo "[\$i/\$n_total][SKIP] \$asmid already cleaned (uncompressed)"
            check_min_length "\$asmid" "\$DEST/\${asmid}.fa"
            printf '%s\\t%s\\n' "\$asmid" "\$DEST/\${asmid}.fa" >> \$MANIFEST
            continue
        fi
        if [ ! -f "\$gz" ]; then
            echo "[\$i/\$n_total][WARN] missing genome for \$asmid: \$gz" >&2
            continue
        fi

        # samples.csv already carries a curated PHYLUM column (ASMID=col1, PHYLUM=col7);
        # prefer it over a live taxonkit lineage lookup, which silently comes back empty
        # for brand-new taxids that haven't propagated into the local NCBI taxdump yet
        # (e.g. provisional taxids for freshly-deposited 'sp.' assemblies). PHYLUM is a
        # NAME (e.g. "Ascomycota"), not a taxid, so it still has to go through
        # taxonkit name2taxid -- AAFTF fcs_gx_purge -t requires a numeric NCBI taxid.
        phylum_name=\$(awk -F',' -v id="\$asmid" '\$1==id{print \$7; exit}' ${samples_csv})
        phylum=""
        if [ -n "\$phylum_name" ]; then
            phylum=\$(echo "\$phylum_name" | \$SING taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
            echo "[\$i/\$n_total][INFO] \$asmid taxonid=\$taxonid phylum_name=\$phylum_name -> phylum_taxid=\$phylum (from samples.csv PHYLUM)"
        fi
        if [ -z "\$phylum" ]; then
            phylum=\$(echo \$taxonid | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage | \$SING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | \$SING taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
            if [ -z "\$phylum" ]; then
                phylum=\$(echo \$taxonid | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage | \$SING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | \$SING taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
            fi
            echo "[\$i/\$n_total][INFO] \$asmid taxonid=\$taxonid phylum_taxid=\$phylum (from taxonkit lineage, samples.csv PHYLUM was blank or unresolvable)"
        fi
        if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
            echo "[\$i/\$n_total][WARN] could not resolve a numeric phylum taxid for \$asmid (taxonid=\$taxonid, samples.csv PHYLUM=\$phylum_name); got '\$phylum'. fcs_gx_purge will likely fail" >&2
            echo "[\$i/\$n_total][WARN] if this taxonid is recent/newly-deposited, TAXONKIT_DB=\$TAXONKIT_DB may be stale; try refreshing it (rm -rf \$TAXONKIT_DB && rerun SETUP_TAXONDB) before assuming the lookup itself is broken" >&2
        fi

        pigz -dc "\$gz" > \$SCRATCH/\${asmid}.raw.fa
        if \$SING AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
            -i \$SCRATCH/\${asmid}.raw.fa --cpus ${task.cpus} \
            -o \$SCRATCH/\${asmid}.purge.fasta \
            -t "\$phylum" -w \$SCRATCH/\${asmid}.fcs_report ; then
            cat \$SCRATCH/\${asmid}.purge.fasta | ${params.clean_script} --len ${params.min_contig_len} > \$SCRATCH/\${asmid}.clean.fa \\
                && pigz -c \$SCRATCH/\${asmid}.clean.fa > \${target}.tmp \\
                && mv \${target}.tmp \$target
            rm -f \$SCRATCH/\${asmid}.clean.fa
            echo "[\$i/\$n_total][OK] \$asmid -> \$target (\$(du -sh \$target | cut -f1))"
            check_min_length "\$asmid" "\$target"
            pigz -f \$SCRATCH/\${asmid}.purge.fasta
            [ -f \$SCRATCH/\${asmid}.purge.fcs_gx-taxonomy.tsv ] && pigz -f \$SCRATCH/\${asmid}.purge.fcs_gx-taxonomy.tsv
            mv \$SCRATCH/\${asmid}.purge.fasta.gz \$DEST/clean/ 2>/dev/null || true
            [ -f \$SCRATCH/\${asmid}.purge.fcs_gx-taxonomy.tsv.gz ] && mv \$SCRATCH/\${asmid}.purge.fcs_gx-taxonomy.tsv.gz \$DEST/clean/
            printf '%s\\t%s\\n' "\$asmid" "\$target" >> \$MANIFEST
        else
            echo "[\$i/\$n_total][FAIL] fcs_gx_purge failed for \$asmid" >&2
        fi
        rm -f \$SCRATCH/\${asmid}.raw.fa \$SCRATCH/\${asmid}.purge.fasta
    done < batch.tsv

    echo "[INFO] batch ${task.index} complete: \$(grep -c . \$MANIFEST || echo 0) cleaned genomes in manifest"
    """

    stub:
    def batch_tsv = items.collect { row -> "${row[1]}\t${row[8]}\t${row[9]}" }.join('\n')
    """
    DEST=${launchDir}/input_clean_genomes
    mkdir -p \$DEST/clean
    MANIFEST=clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    : > TO_ADD_TO_SUPRESS.csv
    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF
    while IFS=\$'\\t' read -r asmid gz taxonid; do
        [ -z "\$asmid" ] && continue
        echo ">stub_\${asmid}" | pigz -c > \$DEST/\${asmid}.fa.gz
        touch \$DEST/clean/\${asmid}.purge.fasta \$DEST/clean/\${asmid}.purge.fcs_gx-taxonomy.tsv
        printf '%s\\t%s\\n' "\$asmid" "\$DEST/\${asmid}.fa.gz" >> \$MANIFEST
    done < batch.tsv
    """
}
