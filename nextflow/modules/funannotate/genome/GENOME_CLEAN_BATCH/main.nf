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

    output:
    path "clean_batch_*.manifest.tsv", emit: manifest

    script:
    def batch_tsv = items.collect { row -> "${row[1]}\t${row[8]}\t${row[9]}" }.join('\n')
    """
    set -uo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"

    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    TAXONKIT_DB=${taxondb}
    DEST=${launchDir}/input_clean_genomes
    mkdir -p \$DEST/clean

    MANIFEST=clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST

    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF

    n_total=\$(grep -c . batch.tsv || true)
    echo "[INFO] batch ${task.index}: \$n_total genomes to consider"

    # Stage the FCS-GX DB into /dev/shm ONCE for the whole batch (~30 min). Keep it
    # across the loop; remove the RAM copy ourselves when the batch finishes.
    export FCS_GX_KEEP_SHM=1
    source ${launchDir}/scripts/setup_fcs_shm.sh
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
            printf '%s\\t%s\\n' "\$asmid" "\$target" >> \$MANIFEST
            continue
        elif [ -s "\$DEST/\${asmid}.fa" ]; then
            echo "[\$i/\$n_total][SKIP] \$asmid already cleaned (uncompressed)"
            printf '%s\\t%s\\n' "\$asmid" "\$DEST/\${asmid}.fa" >> \$MANIFEST
            continue
        fi
        if [ ! -f "\$gz" ]; then
            echo "[\$i/\$n_total][WARN] missing genome for \$asmid: \$gz" >&2
            continue
        fi

        module load taxonkit
        phylum=\$(echo \$taxonid | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
        if [ -z "\$phylum" ]; then
            phylum=\$(echo \$taxonid | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
        fi
        module unload taxonkit
        echo "[\$i/\$n_total][INFO] \$asmid taxonid=\$taxonid phylum=\$phylum"

        module load AAFTF
        pigz -dc "\$gz" > \$SCRATCH/\${asmid}.raw.fa
        if AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
            -i \$SCRATCH/\${asmid}.raw.fa --cpus ${task.cpus} \
            -o \$SCRATCH/\${asmid}.purge.fasta \
            -t "\$phylum" -w \$SCRATCH/\${asmid}.fcs_report ; then
            cat \$SCRATCH/\${asmid}.purge.fasta | ${params.clean_script} --len ${params.min_contig_len} > \$SCRATCH/\${asmid}.clean.fa \\
                && pigz -c \$SCRATCH/\${asmid}.clean.fa > \${target}.tmp \\
                && mv \${target}.tmp \$target
            rm -f \$SCRATCH/\${asmid}.clean.fa
            echo "[\$i/\$n_total][OK] \$asmid -> \$target (\$(du -sh \$target | cut -f1))"
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
