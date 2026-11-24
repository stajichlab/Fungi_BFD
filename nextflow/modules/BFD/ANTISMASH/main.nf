include { hashBucketForType } from '../../common/utils.nf'

// Secondary-metabolite biosynthetic gene cluster (BGC) prediction. Runs
// directly off BFD's existing gff3 + genome-fasta inputs via
// --genefinding-gff3 (no GBK conversion step, unlike
// ~/projects/nf/nf_funannotate1's ANTISMASH_RUN, which starts from a
// funannotate-produced .gbk). Reuses the already-populated shared DB at
// params.antismash_dbdir (/bigdata/stajichlab/shared/lib/antismash_db,
// 9.6 GB, confirmed populated 2026-08-20) and the prebuilt SIF in the shared
// singularity cache -- no new download needed.
//
// BATCHED: one task handles params.antismash_batch_size genomes (default 10) and
// runs params.antismash_batch_parallel of them at a time (default 2), each with
// task.cpus / antismash_batch_parallel CPUs. A genome takes 13-19 min in lite mode at
// 8 CPUs (6 genomes measured 2026-11-21), so 10 genomes at 2 in parallel is
// 5 rounds, about 65-95 min, which fits the 2 h limit of the 'short' queue.
// Per-genome time at other CPU counts and the slowest genomes are NOT measured.
//
// PER-GENOME STORE, no storeDir: each genome is checked, run and written to
// params.outdir/antismash/<bucket>/ by the task itself, as soon as it finishes
// (function one() below). So a failed batch keeps the genomes that finished, and a
// retried task (more memory, see conf) skips them. The write order makes the
// regions table the completion marker: run.json, cds, known, known_pairs (and the
// raw JSON in full mode) first, each through a temp file and mv; regions last.
// A genome whose regions table exists is complete. One failed genome still makes
// the task exit non-zero after the other genomes in the batch have finished.
// The caller (BFD_FUNCTIONAL) also drops genomes whose regions table exists
// before batching, so a rerun builds new batches from the missing genomes only.
// Output `stored` is a text file listing the stored paths of every genome in the
// batch (also the ones skipped as already done).
//
// Two modes (params.antismash_mode, default 'lite'):
//   lite -- cluster detection + MIBiG known-cluster comparison only. Drops
//           --fullhmmer, --clusterhmmer (Pfam annotation of genes; they do not
//           detect clusters) and --pfam2go (GO, unused). Measured on 6 fungal
//           genomes 2026-11-21: 13-19 min/genome at 8 CPUs, against 109 min for
//           the full flag set on Af293, with identical regions on Af293 (44 = 44,
//           same coordinates and classes). The raw antiSMASH output (~240 MB per
//           genome) is NOT kept. Only small gzip TSV tables are:
//             <LOCUSTAG>.antismash.regions.tsv.gz  location, class, contig edge, core
//             <LOCUSTAG>.antismash.cds.tsv.gz      genes inside regions (locus tag,
//                                            coordinates, strand, kind, function)
//             <LOCUSTAG>.antismash.known.tsv.gz    top-5 MIBiG hits per region (with entry size, synteny)
//             <LOCUSTAG>.antismash.known_pairs.tsv.gz  our gene -> MIBiG gene hits (identity, coverage)
//             <LOCUSTAG>.antismash.run.json        antiSMASH and MIBiG versions, mode
//   full -- the original flag set (--fullhmmer --clusterhmmer --cb-general
//           --pfam2go) and the raw JSON kept as <LOCUSTAG>.antismash.json.gz, plus
//           the same tables.
// Outputs are hash-bucketed by LOCUSTAG under params.outdir/antismash/<bucket>/ (same
// scheme as cazy and pfam_hmmscan; width 2).
// The reduction is bin/extract_antismash_subset.py (--minimal). A copy with
// more tables (domains, modules, A-domains) lives in
// DeltaGain_Fungi/analysis/bgc_subset/; keep the two in step.

process RUN_ANTISMASH {
    tag        "${metas instanceof List ? metas[0].locustag + '+' + (metas.size() - 1) : metas.locustag}"
    label      'antismash'
    // No storeDir: the script writes to the store itself (see above).

    input:
        tuple val(metas), path(gff3s), path(genomes)

    output:
        path("stored.txt"), emit: stored   // absolute paths of the stored files, one per line

    script:
    def metaList = metas instanceof List ? metas : [metas]
    def gffList  = gff3s instanceof List ? gff3s : [gff3s]
    def genList  = genomes instanceof List ? genomes : [genomes]
    def full = (params.antismash_mode ?: 'lite') == 'full'
    def par  = Math.max(1, (params.antismash_batch_parallel ?: 2) as int)
    def cpus_each = Math.max(1, task.cpus.intdiv(par))
    def store = "${params.outdir}/antismash"
    // One line per genome: LOCUSTAG GFF3 GENOME TRANSL_TABLE BUCKET (xargs -n 5).
    def rows = [metaList, gffList, genList].transpose().collect { m, g, d ->
        "${m.locustag} ${g.name} ${d.name} ${(m.transl_table ?: '1') as String} ${hashBucketForType('antismash', m.locustag as String)}"
    }.join('\n')
    // Every file a complete genome has in the store.
    def suffixes = ['regions.tsv.gz', 'cds.tsv.gz', 'known.tsv.gz', 'known_pairs.tsv.gz', 'run.json'] + (full ? ['json.gz'] : [])
    def stored = metaList.collect { m ->
        def b = hashBucketForType('antismash', m.locustag as String)
        suffixes.collect { s -> "${store}/${b}/${m.locustag}.antismash.${s}" }
    }.flatten().join('\n')
    // antiSMASH 8.0.4 has no genetic-code option and reads the GFF3 as code 1 (CUG = Leu).
    // For TRANSL_TABLE != 1 (samples.csv; 1,164 genomes use 12), tag every CDS line with
    // transl_table so translations match BFD. Measured on D. hansenii CBS767 (table 12):
    // 3,923 of 5,951 CDS differ without the tag, 0 with it. Table 26 is untested.
    // lite: no --cb-general. Measured 2026-11-22 on 3 genomes: 3.8x-6.9x faster, regions,
    // known and known_pairs identical. Only the unstored blast_top table needs it.
    // --no-*: smaller raw output in $SCRATCH (206 -> 102 MB on A1163); no time saving.
    def flags = full ? '--fullhmmer --clusterhmmer --cb-general --cb-knownclusters --pfam2go'
                     : '--cb-knownclusters --no-zip-output --no-summary-gbk --no-region-gbks'
    """
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    export MPLCONFIGDIR=\$TMPDIR/mpl; mkdir -p \$MPLCONFIGDIR
    export AS_DB='${params.antismash_dbdir}'
    export AS_SIF='${params.antismash_standalone_sif}'
    export AS_TAXON='${params.antismash_taxon}'
    export AS_CPUS=${cpus_each}
    export AS_MODE=${full ? 'full' : 'lite'}
    export AS_EXTRACT='${projectDir}/bin/extract_antismash_subset.py'
    export AS_STORE='${store}'
    export AS_ROOT=\$(mktemp -d -p \$TMPDIR)

    # Process one genome: 1. skip if already stored, 2. run antiSMASH in node-local
    # scratch, 3. reduce, then write to the store (temp file + mv; regions table last).
    one() {
        set -euo pipefail
        local lt=\$1 gff=\$2 gen=\$3 tt=\$4 bucket=\$5
        local dest=\$AS_STORE/\$bucket
        if [ -s "\$dest/\$lt.antismash.regions.tsv.gz" ]; then
            echo "[skip] \$lt already stored"
            return 0
        fi
        local w=\$AS_ROOT/\$lt
        mkdir -p \$w/in \$w/out
        local gffp=\${gff%.gz} genp=\${gen%.gz}
        # This tool may not read gzip, so inflate into scratch (FASTA/GFF3 may be .gz).
        if [[ \$gff == *.gz ]]; then pigz -dc "\$gff" > \$w/in/\$gffp; else cp -L "\$gff" \$w/in/\$gffp; fi
        if [[ \$gen == *.gz ]]; then pigz -dc "\$gen" > \$w/in/\$genp; else cp -L "\$gen" \$w/in/\$genp; fi
        if [ "\$tt" != "1" ]; then
            awk -F'\\t' -v OFS='\\t' -v t="\$tt" '\$0 !~ /^#/ && \$3=="CDS" {\$9=\$9";transl_table="t} {print}' \$w/in/\$gffp > \$w/in/tt.gff3
            mv \$w/in/tt.gff3 \$w/in/\$gffp
        fi
        apptainer exec --bind "\$AS_DB:\$AS_DB,\$TMPDIR:\$TMPDIR" "\$AS_SIF" \\
            antismash --taxon "\$AS_TAXON" \\
                --databases "\$AS_DB" \\
                --output-dir \$w/out/r \\
                --genefinding-gff3 \$w/in/\$gffp \\
                ${flags} \\
                -c \$AS_CPUS \\
                \$w/in/\$genp
        local J
        J=\$(ls \$w/out/r/*.json | head -1)
        # Fail loudly if antiSMASH wrote no JSON: an empty result must not look like "no clusters".
        test -s "\$J"
        python3 "\$AS_EXTRACT" "\$J" \$w/sub \$lt --minimal --mode \$AS_MODE
        mkdir -p "\$dest"
        ${full ? 'pigz -c "\$J" > \$w/sub/\$lt.json.gz' : 'true'}
        # Each file: copy to a temp name in the store dir, then rename (atomic on one filesystem).
        put() { cp "\$1" "\$dest/.\$2.tmp.\$\$" && mv "\$dest/.\$2.tmp.\$\$" "\$dest/\$2"; }
        put \$w/sub/\$lt.run.json  \$lt.antismash.run.json
        ${full ? 'put \$w/sub/\$lt.json.gz \$lt.antismash.json.gz' : 'true'}
        for t in cds known known_pairs; do put \$w/sub/\$lt.\$t.tsv.gz \$lt.antismash.\$t.tsv.gz; done
        put \$w/sub/\$lt.regions.tsv.gz \$lt.antismash.regions.tsv.gz   # last: marks the genome complete
        rm -rf \$w
    }
    export -f one

    cat > genomes.list <<'LIST_EOF'
${rows}
LIST_EOF

    # xargs keeps running the other genomes when one fails and exits 123 at the end.
    xargs -a genomes.list -P ${par} -n 5 bash -c 'one "\$@"' _
    rm -rf \$AS_ROOT

    # Every genome of the batch must now be in the store.
    cat > stored.txt <<'STORED_EOF'
${stored}
STORED_EOF
    while read -r f; do test -s "\$f" || { echo "[ERROR] missing \$f" >&2; exit 1; }; done < stored.txt
    """

    stub:
    // Stub writes under stub_store/ in the task dir, never into the real store.
    def metaList = metas instanceof List ? metas : [metas]
    """
    : > stored.txt
    ${ metaList.collect { m ->
        def b = hashBucketForType('antismash', m.locustag as String)
        "mkdir -p stub_store/${b} && for t in regions cds known known_pairs; do printf '' | gzip > stub_store/${b}/${m.locustag}.antismash.\$t.tsv.gz; echo \"\$PWD/stub_store/${b}/${m.locustag}.antismash.\$t.tsv.gz\" >> stored.txt; done"
    }.join('\n') }
    """
}
