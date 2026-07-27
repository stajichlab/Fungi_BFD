process REPEATMASK_STRAIN {
    tag    "${asmid} (${species})"
    label  'earlgrey_repeatmask'

    storeDir { "${params.outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}/strains" }

    input:
        tuple val(asmid), val(species), path(genome), path(library)

    output:
        tuple val(asmid), path("${asmid}.masked.fasta.gz")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load RepeatMasker/${params.repeatmasker_version}

    GENOME_IN="${genome}"
    case "${genome}" in
        *.gz) echo "[INFO] Inflating compressed genome ${genome}"; gzip -dc "${genome}" > genome_input.fa; GENOME_IN=genome_input.fa ;;
    esac

    mkdir -p rmask_out
    RepeatMasker \
        -lib ${library} \
        -xsmall \
        -pa ${task.cpus} \
        -dir rmask_out \
        "\$GENOME_IN"

    MASKED="rmask_out/\$(basename "\$GENOME_IN").masked"
    if [ -f "\$MASKED" ]; then
        gzip -c "\$MASKED" > ${asmid}.masked.fasta.gz
    else
        echo "WARN: no repeats masked for ${asmid}; using unmasked genome" >&2
        gzip -c "\$GENOME_IN" > ${asmid}.masked.fasta.gz
    fi
    """

    stub:
    """
    printf '>stub_${asmid}\\nacgtACGTacgt\\n' | gzip -c > ${asmid}.masked.fasta.gz
    """
}
