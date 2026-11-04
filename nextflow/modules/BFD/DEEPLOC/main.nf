include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

// DeepLoc 2.1 (DTU, academic license): multi-label subcellular localization,
// sorting signals (signal peptide, mitochondrial/chloroplast transit peptide,
// NLS/NES, PTS) and membrane type, for every protein. Replaces TargetP 2.0 in
// -profile BFD (TargetP only gave SP/mTP/noTP).
//
// No Lmod module is usable from a container-agnostic task, so this runs the
// lab-built SIF (~/projects/containers/deeploc_container, deeploc-2.1.def).
// Both model backbones (Fast = ESM1b, Accurate = ProtT5) are baked into the
// image and HF_HUB_OFFLINE=1 is set, so tasks download nothing.
//
// Model notes (DeepLoc2/deeploc2.py): proteins are run one at a time (no
// batching). Fast clips proteins > 1022 aa to the first + last 511 aa;
// Accurate clips > 4000 aa to the first + last 2000 aa.
//
// Output: DeepLoc's own results CSV, renamed to a fixed name and gzipped. The
// model is part of the file name so Fast and Accurate results never mix in
// storeDir or in MERGE_DEEPLOC's glob.
process RUN_DEEPLOC {
    tag        "${meta.locustag}"
    label      'deeploc'
    storeDir   { "${params.outdir}/deeploc/${hashBucketForType('deeploc', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.deeploc_${params.deeploc_model.toLowerCase()}.csv.gz"), emit: results

    script:
    def proteins_plain = plainName(proteins)
    def model  = params.deeploc_model.toLowerCase() == 'accurate' ? 'Accurate' : 'Fast'
    def device = params.deeploc_gpu ? 'cuda' : 'cpu'
    def out    = "${meta.locustag}.deeploc_${model.toLowerCase()}.csv.gz"
    """
    ${gzInflate(proteins)}
    module load apptainer
    if [ "${params.deeploc_gpu}" = "true" ]; then nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader || true; fi
    # The protein file may be a symlink to a shared genome_annotation dir;
    # bind its real parent too, or apptainer cannot follow it.
    PROT_REAL_DIR=\$(dirname \$(readlink -f ${proteins_plain}))
    export OMP_NUM_THREADS=${task.cpus}
    rm -rf deeploc_out
    apptainer exec ${params.deeploc_gpu ? '--nv' : ''} \\
        --bind \${PWD}:\${PWD},\${PROT_REAL_DIR}:\${PROT_REAL_DIR} \\
        ${params.deeploc_sif} \\
        deeploc2 -f \${PWD}/${proteins_plain} -o \${PWD}/deeploc_out -m ${model} -d ${device}

    shopt -s nullglob
    RES=(deeploc_out/results_*.csv)
    if [ \${#RES[@]} -ne 1 ]; then
        echo "[ERROR] expected 1 DeepLoc results CSV, found \${#RES[@]}" >&2; exit 1
    fi
    # DeepLoc keys results by FASTA id, so duplicate ids would collapse.
    # Require one row per input sequence.
    NSEQ=\$(grep -c '^>' ${proteins_plain} || true)
    NROW=\$(( \$(wc -l < "\${RES[0]}") - 1 ))
    if [ "\$NSEQ" -ne "\$NROW" ]; then
        echo "[ERROR] ${meta.locustag}: \$NSEQ input sequences but \$NROW DeepLoc rows" >&2; exit 1
    fi
    pigz -c "\${RES[0]}" > ${out}
    rm -rf deeploc_out
    """

    stub:
    def out = "${meta.locustag}.deeploc_${params.deeploc_model.toLowerCase()}.csv.gz"
    """
    printf 'Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome,Peripheral,Transmembrane,Lipid anchor,Soluble\\n' | gzip > ${out}
    """
}
