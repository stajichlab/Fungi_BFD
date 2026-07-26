process SETUP_SYMLINKS {
    label 'setup'

    input:
        path(rows_file)

    output:
        path("symlink_manifest.txt"), emit: manifest

    script:
    """
    mkdir -p "${params.pep_dir}" "${params.cds_dir}" "${params.gff_dir}" \\
             "${params.genome_dir}" "${params.trna_dir}"

    make_link() {
        local target=\$1 linkname=\$2
        if [ ! -e "\$target" ]; then
            echo "[WARN] source not found, skipping: \$target" >&2
            return 0
        fi
        if [[ ! -L "\$linkname" || ! -e "\$linkname" ]]; then
            ln -sfn "\$target" "\$linkname"
            echo "[INFO] linked \$linkname -> \$target"
        else
            echo "[INFO] symlink already valid, skipping: \$linkname"
        fi
    }

    : > symlink_manifest.txt
    while IFS=\$'\\t' read -r basename locustag; do
        src="${params.genome_annotation}/\${basename}/predict_results"
        misc="${params.genome_annotation}/\${basename}/predict_misc"

        if [ ! -d "\$src" ]; then
            echo "[WARN] predict_results not found for \${basename}: \$src" >&2
            continue
        fi
        make_link "\$src/\${basename}.proteins.fa"        "${params.pep_dir}/\${basename}.proteins.fa"
        make_link "\$src/\${basename}.cds-transcripts.fa" "${params.cds_dir}/\${basename}.cds-transcripts.fa"
        make_link "\$src/\${basename}.gff3"               "${params.gff_dir}/\${basename}.gff3"
        make_link "\$src/\${basename}.scaffolds.fa"       "${params.genome_dir}/\${basename}.scaffolds.fa"

        if [ -f "\$misc/trnascan.no-overlaps.gff3" ]; then
            make_link "\$misc/trnascan.no-overlaps.gff3" "${params.trna_dir}/\${basename}.trna.gff3"
        else
            echo "[INFO] no trnascan GFF3 for \${basename}, skipping trna symlink"
        fi
        echo "\${locustag}\t\${basename}" >> symlink_manifest.txt
    done < ${rows_file}
    echo "[INFO] SETUP_SYMLINKS: \$(wc -l < symlink_manifest.txt | tr -d ' ') species linked"
    """

    stub:
    """
    awk '{print \$1"\\t"\$2}' ${rows_file} > symlink_manifest.txt
    echo "[STUB] SETUP_SYMLINKS: \$(wc -l < symlink_manifest.txt | tr -d ' ') species"
    """
}
