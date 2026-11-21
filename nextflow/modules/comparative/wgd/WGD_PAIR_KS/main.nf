// WGD_PAIR_KS — pair-level dS for wgd dmd paralog families, as a replacement for WGD_KSD.
//
// For every gene pair in each family (bin/wgd_pair_ks.py): global protein alignment, back-
// translation, protein identity + NG86 dS, a cheap pre-filter (identity >= wgd_prefilter_identity,
// >= wgd_prefilter_codons aligned codons, optional NG86 dS cap), then one pairwise codeml per
// surviving pair. No family-wide MSA. Writes the same ks.tsv columns that MERGE_WGD_KSD reads,
// plus a `node` column (neighbour-joining tree on protein identity) and `status`.
// Evaluation and limits: eval/sim_controls/report.md (2026-11-01 sections).
//
// Bundling: one task runs params.wgd_pair_bundle genomes one after the other (default 20), to
// reach about 1-1.5 h per task. Output goes flat into wgd_ksd_pp/ (or wgd_ksd_pp_full/ when
// params.wgd_speedups = false: pre-filter off). A genome that fails is listed in
// bundle_<locustag>.failed.txt (copied into the table folder by the script); the merge step stops with an error when any
// such file exists. The task itself fails only when every genome in the bundle failed.
// Genetic code: meta.transl_table (samples.csv TRANSL_TABLE) -> --table; codeml icode 8 for
// table 12; other non-standard tables are NOT validated and the task fails for them unless
// params.wgd_allow_untested_table = true.
// dN/dS estimator: params.wgd_dnds_method = yn00 (kaks, Yang & Nielsen 2000; about 30x faster per pair) or codeml. kaks is in the
// container params.wgd_pair_sif (wgd2_complete with subopt-kaks).
process WGD_PAIR_KS {
    tag        "${metas.size()} genomes: ${metas[0].id}"
    label      'comparative_wgd_pair'

    // results of different estimators never share a folder: wgd_ksd_pp (codeml), wgd_ksd_pp_yn00; _full = wgd_speedups false
    storeDir   { "${params.outdir}/${params.wgd_speedups.toString().toBoolean() ? 'wgd_ksd_pp' : 'wgd_ksd_pp_full'}${params.wgd_dnds_method == 'yn00' ? '_yn00' : ''}" }

    input:
    tuple val(metas), path(families, stageAs: 'fam_?/*'), path(cds, stageAs: 'cds_?/*')

    output:
    path({ metas.collect { "${it.id}.cds-transcripts.fa.tsv.ks.tsv.gz" } }), emit: ks, optional: true
    path({ metas.collect { "${it.id}.cds-transcripts.fa.tsv.ks.tsv.gz.nodes.tsv.gz" } }), emit: nodes, optional: true
    // Required marker, named after the exact set of genomes. All other outputs are optional (a genome may fail),
    // and with only optional outputs Nextflow treats a storeDir as "already done" as soon as one file of the bundle
    // exists, which silently skipped new genomes of a bundle and every re-run of a failed one. Failure lists are
    // copied straight into the storeDir by the script (not declared), so they do not count as outputs either.
    path("bundle_${metas.collect { it.id }.sort().join(',').md5()}.done"), emit: done

    script:
    def fast    = params.wgd_speedups.toString().toBoolean()
    def ident   = fast ? params.wgd_prefilter_identity : 0
    def ngcod   = fast ? params.wgd_prefilter_codons : 0
    def maxngds = (fast && params.wgd_prefilter_ngds) ? "--max-ngds ${params.wgd_prefilter_ngds}" : ''
    // identical-CDS collapse is a speed-up: off when wgd_speedups=false or wgd_collapse_identical=false
    def nocoll  = (fast && params.wgd_collapse_identical.toString().toBoolean()) ? '' : '--no-collapse'
    // node sampling: K pairs per node go to codeml (0 = all pairs). A speed-up: off when wgd_speedups=false.
    def ns      = fast ? (params.wgd_node_sample as int) : 0
    // a node draws pairs in rounds until it has ns pairs inside the gate (dS cut, minimum codons) or wgd_node_max_draws were run
    def refine  = (params.wgd_dnds_method == 'yn00' && (params.wgd_refine_codeml_max as double) > 0) ? "--refine-codeml-min ${params.wgd_refine_codeml_min} --refine-codeml-max ${params.wgd_refine_codeml_max}" : ''
    def nsflag  = ns > 0 ? "--node-sample ${ns} --node-max-draws ${params.wgd_node_max_draws} --gate-max-ds ${params.wgd_max_ds ?: 0} --gate-min-codons ${params.wgd_min_pair_codons ?: 0}" : ''
    def allow   = params.wgd_allow_untested_table.toString().toBoolean() ? '1' : '0'
    def ids     = metas.collect { it.id }.join(' ')
    def tabs    = metas.collect { (it.transl_table ?: '1').toString() }.join(' ')
    def fams    = (families instanceof List ? families : [families]).collect { it.toString() }.join(' ')
    def cdss    = (cds instanceof List ? cds : [cds]).collect { it.toString() }.join(' ')
    """
    module load apptainer squashfuse fuse
    export TMPDIR=\${SCRATCH:-/tmp}
    PWD_REAL=\$(realpath \${PWD})
    SING_BINDS="--bind \${PWD_REAL}:\${PWD_REAL},${projectDir}:${projectDir},/bigdata:/bigdata,/rhome:/rhome,\$TMPDIR:\$TMPDIR"
    IDS=(${ids}); TABS=(${tabs}); FAMS=(${fams}); CDSS=(${cdss})
    FAILED=bundle_${metas[0].locustag}.failed.txt
    STORE=${params.outdir}/${params.wgd_speedups.toString().toBoolean() ? 'wgd_ksd_pp' : 'wgd_ksd_pp_full'}${params.wgd_dnds_method == 'yn00' ? '_yn00' : ''}
    nok=0
    for i in \${!IDS[@]}; do
        id=\${IDS[\$i]}; tab=\${TABS[\$i]}; fam=\${FAMS[\$i]}; cds=\${CDSS[\$i]}
        out=\${id}.cds-transcripts.fa.tsv.ks.tsv.gz
        # already stored by an earlier run (another bundle): do not redo it
        if [ -s "\$STORE/\$out" ]; then echo "\$id: already stored, skipped" >&2; nok=\$((nok+1)); continue; fi
        icode=0; table=1
        case \$tab in
            1) ;;
            12) table=12; icode=8 ;;
            *) if [ "${allow}" = 1 ]; then table=\$tab; icode=0; else
                 echo "\$id: genetic code \$tab not validated; skipped" >> \$FAILED; continue; fi ;;
        esac
        if apptainer exec \${SING_BINDS} --env PYTHONNOUSERSITE=1 ${params.wgd_pair_sif} \\
            python3 ${projectDir}/bin/wgd_pair_ks.py --cds \$cds --families \$fam --out \$out \\
            --threads ${task.cpus} --min-identity ${ident} --min-ng-codons ${ngcod} ${maxngds} \\
            --dnds-method ${params.wgd_dnds_method} ${refine} ${nsflag} ${nocoll} --table \$table --icode \$icode --tmpdir \$TMPDIR < /dev/null \\
            && [ -s \$out ]; then
            nok=\$((nok+1)); cp \$out.status.txt \${id}.pair_status.txt 2>/dev/null || true
        else
            echo "\$id: wgd_pair_ks.py failed" >> \$FAILED; rm -f \$out
        fi
    done
    echo "genomes ok: \$nok of \${#IDS[@]}" >&2
    if [ -s \$FAILED ]; then mkdir -p "\$STORE" && cp \$FAILED "\$STORE/"; fi
    [ \$nok -gt 0 ] || { echo "ERROR: every genome in the bundle failed" >&2; exit 1; }
    touch bundle_${metas.collect { it.id }.sort().join(',').md5()}.done
    """

    stub:
    """
    for id in ${metas.collect { it.id }.join(' ')}; do printf 'pair\\tdS\\n' | gzip > \${id}.cds-transcripts.fa.tsv.ks.tsv.gz; done
    touch bundle_${metas.collect { it.id }.sort().join(',').md5()}.done
    """
}
