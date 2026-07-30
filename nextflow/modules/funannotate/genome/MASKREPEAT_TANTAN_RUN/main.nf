// Soft-mask each assembly using funannotate mask with tantan.
// storeDir caches the masked FASTA alongside the clean genome.
process MASKREPEAT_TANTAN_RUN {
    tag "$asmid"

    storeDir "${launchDir}/input_clean_genomes"

    cpus   8
    memory '16 GB'
    time   '2h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(taxonid)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.masked.fasta.gz"), val(taxonid), emit: masked

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate
    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac
    funannotate mask -i "\$GENOME_IN" -o ${asmid}.masked.fasta -m tantan --cpus ${task.cpus}
    pigz -f ${asmid}.masked.fasta
    """

    stub:
    """
    echo ">stub_${asmid}_masked" | pigz -c > ${asmid}.masked.fasta.gz
    """
}
