process EARLGREY_BUILD_LIB {
    tag    "${species} (${rep_asmid})"
    label  'earlgrey_build'

    storeDir { "${params.outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}" }

    input:
        tuple val(species), val(rep_asmid), path(genome)

    output:
        tuple val(species), val(rep_asmid), path("${rep_asmid}.families.fa"), emit: lib
        tuple val(rep_asmid), path("${rep_asmid}.masked.fasta.gz"), emit: masked
        path("*_RepeatLandscape"), emit: landscape, optional: true
        path("*_summaryFiles"),    emit: summary,   optional: true

    script:
    def sp_safe = species.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def mem_mb = task.memory ? (task.memory.toMega() * 0.9) as long : 3200
    def egdir = "${params.earlgrey_workdir}/${sp_safe}"
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load earlgrey/${params.earlgrey_version}

    GENOME_IN="${genome}"
    case "${genome}" in
        *.gz) echo "[INFO] Inflating compressed genome ${genome}"; gzip -dc "${genome}" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
    esac

    mkdir -p ${egdir}

    earlGrey \
        -g "\$GENOME_IN" \
        -s ${sp_safe} \
        -o ${egdir} \
        -r ${params.repeat_taxon} \
        -d yes \
        -c yes \
        -v yes \
        -q yes \
        -M ${mem_mb} \
        -t ${task.cpus}

    lib=\$(find ${egdir} -name '*-families.fa' | head -n1)
    soft=\$(find ${egdir} -name '*.softmasked.fasta' | head -n1)

    if [ -z "\$lib" ] || [ -z "\$soft" ]; then
        echo "ERROR: EarlGrey outputs not found (lib=\$lib soft=\$soft)" >&2
        find ${egdir} -maxdepth 3 -type f >&2
        exit 1
    fi

    cp "\$lib"  ${rep_asmid}.families.fa
    gzip -c "\$soft" > ${rep_asmid}.masked.fasta.gz

    for d in \$(find ${egdir} -maxdepth 3 -type d \\( -name '*_RepeatLandscape' -o -name '*_summaryFiles' \\)); do
        cp -r "\$d" ./
    done

    gz=\$(command -v pigz >/dev/null 2>&1 && echo "pigz -p ${task.cpus}" || echo "gzip")
    for d in ./*_summaryFiles; do
        [ -d "\$d" ] || continue
        find "\$d" -type f ! -name '*.gz' ! -name '*.bgz' ! -name '*.zip' -print0 \
            | xargs -0 -r \$gz -f
    done

    rm -rf ${egdir}
    """

    stub:
    """
    printf '>stub_family_1#Unknown\\nACGTACGTACGT\\n' > ${rep_asmid}.families.fa
    printf '>stub_${rep_asmid}\\nacgtACGTacgt\\n' | gzip -c > ${rep_asmid}.masked.fasta.gz

    sp_safe=\$(echo '${species}' | sed 's/[^A-Za-z0-9._-]\\+/_/g')
    mkdir -p "\${sp_safe}_RepeatLandscape" "\${sp_safe}_summaryFiles"
    printf 'stub landscape\\n' > "\${sp_safe}_RepeatLandscape/stub.txt"
    printf 'stub summary\\n'   > "\${sp_safe}_summaryFiles/stub.txt"
    """
}
