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
    # A plain file is linked. A compressed one (<name>.gz, see SETUP_SYMLINKS) is
    # inflated to the plain destination name, since the comparative tools read
    # ${proteins_out}/*.faa and ${cds_out}/*.cds.fa as plain FASTA.
    stage_one() {
        local src=\$1 dst=\$2
        [ -e "\$dst" ] && return 0
        if [ -f "\$src" ]; then
            ln -sf "\$src" "\$dst"
        elif [ -f "\$src.gz" ]; then
            gzip -dc "\$src.gz" > "\$dst.tmp.\$\$" && mv "\$dst.tmp.\$\$" "\$dst" || { rm -f "\$dst.tmp.\$\$"; return 1; }
        else
            return 2
        fi
    }
    tail -n +2 ${manifest} | while IFS=',' read -r locustag group basename; do
        [ -z "\$locustag" ] && continue
        [ -z "\$basename" ] && continue
        src="${pep_dir}/\${basename}.proteins.fa"
        stage_one "\$src" "${proteins_out}/\${locustag}.faa" || echo "WARN: missing protein file for \${locustag}: \$src[.gz]" >&2
        stage_one "${cds_dir}/\${basename}.cds-transcripts.fa" "${cds_out}/\${locustag}.cds.fa" || true
    done
    n_staged=\$(find "${proteins_out}" -maxdepth 1 -name '*.faa' | wc -l)
    echo "Staged \${n_staged} protein file(s) into ${proteins_out}"
    if [ "\${n_staged}" -eq 0 ]; then
        echo "ERROR: no protein files staged — check pep_dir (${pep_dir}) and BASENAME mapping" >&2
        exit 1
    fi
    """
}
