process STAGE_FILES {
    label    'comparative_setup'
    tag      params.project

    publishDir "${params.outdir}/${params.project}", mode: 'copy', pattern: '*.manifest.csv'

    input:
    path manifest
    val  pep_dir
    val  cds_dir

    output:
    val  "${params.outdir}/${params.project}/proteins", emit: proteins_dir
    val  "${params.outdir}/${params.project}/cds",      emit: cds_dir
    path manifest,                                      emit: manifest

    script:
    def proteins_out = "${params.outdir}/${params.project}/proteins"
    def cds_out      = "${params.outdir}/${params.project}/cds"
    """
    mkdir -p "${proteins_out}" "${cds_out}"
    tail -n +2 ${manifest} | while IFS=',' read -r locustag group basename; do
        [ -z "\$locustag" ] && continue
        [ -z "\$basename" ] && continue
        src="${pep_dir}/\${basename}.proteins.fa"
        dst="${proteins_out}/\${locustag}.faa"
        if [ -f "\$src" ]; then
            [ ! -e "\$dst" ] && ln -sf "\$src" "\$dst"
        else
            echo "WARN: missing protein file for \${locustag}: \$src" >&2
        fi
        src="${cds_dir}/\${basename}.cds-transcripts.fa"
        dst="${cds_out}/\${locustag}.cds.fa"
        [ -f "\$src" ] && [ ! -e "\$dst" ] && ln -sf "\$src" "\$dst"
    done
    n_staged=\$(find "${proteins_out}" -maxdepth 1 -name '*.faa' | wc -l)
    echo "Staged \${n_staged} protein file(s) into ${proteins_out}"
    if [ "\${n_staged}" -eq 0 ]; then
        echo "ERROR: no protein files staged — check pep_dir (${pep_dir}) and BASENAME mapping" >&2
        exit 1
    fi
    """
}
