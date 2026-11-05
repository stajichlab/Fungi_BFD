// WGD_DMD_BUNDLE — wgd dmd for params.wgd_dmd_bundle genomes per task (about 1 h per task).
//
// Same command as WGD_DMD (all-vs-all DIAMOND on the translated CDS + MCL families) but one task
// runs many genomes one after the other: wgd dmd takes seconds to 2 min per genome, far too short
// for one SLURM job each. Output goes flat into wgd_dmd/ (exact file names, so a storeDir hit only
// happens when every genome of the bundle is present). A genome that fails is listed in
// dmd_<first locustag>.failed.txt in the same folder; the task fails only when every genome failed.
// Use with --wgd_stage dmd to build and cache the families for a whole sample sheet first.
include { plainName } from '../../../common/utils.nf'

process WGD_DMD_BUNDLE {
    tag        "${metas.size()} genomes: ${metas[0].id}"
    label      'comparative_wgd'

    storeDir   { "${params.outdir}/wgd_dmd" }

    input:
    tuple val(metas), path(cds, stageAs: 'cds_?/*')

    output:
    // only the families files are outputs: an input path in the output tuple would put the staged cds_N/ folders into the storeDir
    tuple val(metas), path({ metas.collect { "${it.id}.cds-transcripts.fa.tsv" } }), emit: families, optional: true
    path("dmd_${metas[0].locustag}.failed.txt"), emit: failed, optional: true

    script:
    def files = (cds instanceof List ? cds : [cds]).collect { it.toString() }.join(' ')
    def ids   = metas.collect { it.id }.join(' ')
    """
    module load apptainer squashfuse fuse
    export TMPDIR=\${SCRATCH:-/tmp}
    PWD_REAL=\$(realpath \${PWD})
    SING_BINDS="--bind \${PWD_REAL}:\${PWD_REAL},${projectDir}:${projectDir},/bigdata:/bigdata,/rhome:/rhome,\$TMPDIR:\$TMPDIR"
    FILES=(${files}); IDS=(${ids})
    FAILED=dmd_${metas[0].locustag}.failed.txt
    nok=0
    for i in \${!FILES[@]}; do
        f=\${FILES[\$i]}; id=\${IDS[\$i]}; b=\$(basename \$f); plain=\${b%.gz}
        if [ "\$b" != "\$plain" ]; then pigz -dc \$f > \$plain || { echo "\$id: cannot inflate" >> \$FAILED; continue; }; else cp \$f \$plain; fi
        rm -rf dmd_out
        if apptainer exec \${SING_BINDS} --env PYTHONNOUSERSITE=1 ${params.wgd_sif} wgd dmd \\
              --outdir dmd_out --eval ${params.wgd_mcl_evalue} --inflation ${params.wgd_mcl_inflation} \\
              \$plain < /dev/null && [ -s dmd_out/\$plain.tsv ]; then
            mv dmd_out/\$plain.tsv .; nok=\$((nok+1))
        else
            echo "\$id: wgd dmd failed" >> \$FAILED
        fi
        rm -rf dmd_out \$plain
    done
    echo "genomes ok: \$nok of \${#FILES[@]}" >&2
    [ \$nok -gt 0 ] || { echo "ERROR: every genome in the bundle failed" >&2; exit 1; }
    """

    stub:
    """
    for f in ${(cds instanceof List ? cds : [cds]).collect { plainName(it) }.join(' ')}; do printf 'family\\tgene_a\\tgene_b\\n' > \${f}.tsv; done
    """
}
