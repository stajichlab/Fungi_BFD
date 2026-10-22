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
                # The result may have been uncompressed again: drop a stale .gz link
                [ -L "\$linkname.gz" ] && rm -f "\$linkname.gz"
            elif [ -e "\$target.gz" ]; then
                # Compressed result (.gff3.gz, .cds-transcripts.fa.gz): link under the .gz name
                ln -sfn "\$target.gz" "\$linkname.gz"
                # The plain link left from before compression now dangles: remove it
                [ -L "\$linkname" ] && rm -f "\$linkname"
            else
                # Neither form exists: remove dangling links of either form
                [ -L "\$linkname" ] && [ ! -e "\$linkname" ] && rm -f "\$linkname"
                [ -L "\$linkname.gz" ] && [ ! -e "\$linkname.gz" ] && rm -f "\$linkname.gz"
                echo "[WARN] source not found, skipping: \$target" >&2
            fi
        done

        # tRNA GFF3: link the plain file, or <name>.gz when only that exists (same
        # rule as the pairs above). Silent when neither exists: not every genome has one.
        local trna_target="\$misc/trnascan.no-overlaps.gff3"
        local trna_link="${params.trna_dir}/\${basename}.trna.gff3"
        if [ -e "\$trna_target" ]; then
            ln -sfn "\$trna_target" "\$trna_link"
            [ -L "\$trna_link.gz" ] && rm -f "\$trna_link.gz"
        elif [ -e "\$trna_target.gz" ]; then
            ln -sfn "\$trna_target.gz" "\$trna_link.gz"
            [ -L "\$trna_link" ] && rm -f "\$trna_link"
        else
            [ -L "\$trna_link" ] && [ ! -e "\$trna_link" ] && rm -f "\$trna_link"
            [ -L "\$trna_link.gz" ] && [ ! -e "\$trna_link.gz" ] && rm -f "\$trna_link.gz"
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
