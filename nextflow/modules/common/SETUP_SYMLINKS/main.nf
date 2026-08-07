process SETUP_SYMLINKS {
    label 'setup'

    input:
        path(rows_file)

    output:
        path("symlink_manifest.txt"), emit: manifest

    script:
    def maxjobs = task.ext.symlink_jobs ?: 32
    """
    mkdir -p "${params.pep_dir}" "${params.cds_dir}" "${params.gff_dir}" \\
             "${params.genome_dir}" "${params.trna_dir}"

    genome_annotation_real=\$(realpath "${params.genome_annotation}")

    manifest_dir=\$(mktemp -d ./manifest_parts.XXXXXX)

    # NFS metadata round trips dominate this step, not CPU, so fan the
    # per-row work out across many concurrent background shells rather
    # than doing one `stat`/`ln` at a time in a serial loop. Each worker
    # writes its own manifest fragment (no shared-file locking over NFS,
    # which is unreliable on some mounts) and fragments are concatenated
    # once all workers finish.
    process_row() {
        local basename=\$1 locustag=\$2
        local src="\${genome_annotation_real}/\${basename}/predict_results"
        local misc="\${genome_annotation_real}/\${basename}/predict_misc"

        if [ ! -d "\$src" ]; then
            echo "[WARN] predict_results not found for \${basename}: \$src" >&2
            return 0
        fi

        local pair target linkname
        for pair in \\
            "\$src/\${basename}.proteins.fa|${params.pep_dir}/\${basename}.proteins.fa" \\
            "\$src/\${basename}.cds-transcripts.fa|${params.cds_dir}/\${basename}.cds-transcripts.fa" \\
            "\$src/\${basename}.gff3|${params.gff_dir}/\${basename}.gff3" \\
            "\$src/\${basename}.scaffolds.fa|${params.genome_dir}/\${basename}.scaffolds.fa"
        do
            target=\${pair%%|*}
            linkname=\${pair##*|}
            if [ -e "\$target" ]; then
                ln -sfn "\$target" "\$linkname"
            else
                echo "[WARN] source not found, skipping: \$target" >&2
            fi
        done

        if [ -f "\$misc/trnascan.no-overlaps.gff3" ]; then
            ln -sfn "\$misc/trnascan.no-overlaps.gff3" "${params.trna_dir}/\${basename}.trna.gff3"
        fi

        printf '%s\\t%s\\n' "\$locustag" "\$basename" > "\${manifest_dir}/\${basename}"
    }

    njobs=0
    while IFS=\$'\\t' read -r basename locustag; do
        process_row "\$basename" "\$locustag" &
        njobs=\$((njobs + 1))
        if (( njobs >= ${maxjobs} )); then
            wait -n
            njobs=\$((njobs - 1))
        fi
    done < ${rows_file}
    wait

    cat "\${manifest_dir}"/* > symlink_manifest.txt 2>/dev/null || : > symlink_manifest.txt
    rm -rf "\${manifest_dir}"

    echo "[INFO] SETUP_SYMLINKS: \$(wc -l < symlink_manifest.txt | tr -d ' ') species linked"
    """

    stub:
    """
    awk '{print \$1"\\t"\$2}' ${rows_file} > symlink_manifest.txt
    echo "[STUB] SETUP_SYMLINKS: \$(wc -l < symlink_manifest.txt | tr -d ' ') species"
    """
}
