#!/usr/bin/env nextflow

/*
 * compare_ANI — Average Nucleotide Identity comparison pipeline
 *
 * Groups genomes from samples.csv by a taxonomic rank (--compare GENUS by default)
 * and computes all-vs-all identity within each group, then writes a clustering +
 * outlier report per group.
 *
 * Identity method (--ani_method):
 *   skani    (default)  fast, accurate ANI; sketch-once + sparse triangle.
 *                       Best for clade-level (CLASS/PHYLUM) where fastANI is too slow.
 *   mash                fastest, roughest; MinHash distance -> ANI = 100*(1-dist).
 *   sourmash            MinHash containment-ANI matrix (sourmash compare --ani).
 *   fastani             classic fragment-mapping ANI. Optionally pre-filtered with
 *                       mash (--fastani_prefilter) so fastANI only runs all-vs-all
 *                       *within* mash connected components, not across the whole clade.
 *
 * Sketch caching: skani/mash/sourmash sketch each genome ONCE into
 *   ${params.sketch_cache}/<method>/<params>/  (storeDir).
 * The same sketch is reused across groups and across re-runs / different --compare
 * levels, so a genome is never re-sketched.
 *
 * Usage (from project root):
 *   nextflow run nextflow/compare_ANI.nf \
 *       -c nextflow/nextflow.config -profile ani -resume
 *
 *   # clade level with the default (skani)
 *   nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --taxon PHYLUM:Ascomycota --compare CLASS -resume
 *
 *   # pick a method
 *   nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --ani_method mash -resume
 *
 *   # fastANI with mash prefilter cascade
 *   nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       --ani_method fastani --fastani_prefilter true -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
 *       -stub-run --n_test 3
 */

// ════════════════════════════════════════════════════════════════════════════
// PARAMETER DEFAULTS
// (override via --param value on CLI or in a params-file)
// ════════════════════════════════════════════════════════════════════════════

params.genome_dir             = "${launchDir}/input/dna"
params.genome_suffix          = '.scaffolds.fa'   // alt: '.fa.gz' / '.masked.fasta.gz' with input_clean_genomes/
// Controls how the genome filename stem is derived from samples.csv:
//   'species' (default) → ${SPECIES}_${STRAIN}${genome_suffix}  e.g. Fusarium_oxysporum_CBS123.scaffolds.fa
//   'asmid'             → ${ASMID}${genome_suffix}               e.g. GCF_010015735.1_Aaoar1.fa.gz
// NOTE: input_clean_genomes/ files are gzipped (.fa.gz / .masked.fasta.gz). skani,
// mash and sourmash read gzipped FASTA natively, so no decompression is needed.
// (fastANI does NOT read .gz — use uncompressed input if --ani_method fastani.)
params.genome_name_style      = 'species'
params.compare                = 'GENUS'
params.ani_cluster_threshold  = 95.0
params.ani_outlier_threshold  = 90.0
params.min_group_size         = 2
params.outdir                 = "${launchDir}/results/ANI"

// ── Method selection ─────────────────────────────────────────────────────────
params.ani_method             = 'skani'   // skani | mash | sourmash | fastani

// ── Sketch cache (shared across groups, runs, and --compare levels) ──────────
params.sketch_cache           = "${launchDir}/work/ANI/sketch_cache"

// ── skani ────────────────────────────────────────────────────────────────────
params.skani_preset           = 'medium'  // fast | medium | slow (must match sketch+triangle)
params.skani_min_af           = 15        // minimum aligned fraction (%) to report a pair
params.skani_compression      = 0         // -c override; 0 = use preset default
params.skani_sketch_chunk     = 50        // genomes sketched per skani job (>=1)

// ── mash ──────────────────────────────────────────────────────────────────────
params.mash_kmer              = 21
params.mash_sketch_size       = 10000

// ── sourmash ──────────────────────────────────────────────────────────────────
params.sourmash_kmer          = 21
params.sourmash_scaled        = 1000

// ── fastANI ───────────────────────────────────────────────────────────────────
params.fastani_fraglen        = 3000
params.fastani_kmer           = 16
params.ani_batch_size         = 0         // 0 = single job per group; >0 = upper-triangle batches
params.fastani_prefilter      = false     // true = mash pre-cluster, then fastANI within components
params.prefilter_ani          = 80.0      // mash-ANI%% floor for the prefilter cascade

// ════════════════════════════════════════════════════════════════════════════
// HELPERS
// ════════════════════════════════════════════════════════════════════════════

// skani preset → CLI flag (medium is the tool default → no flag)
def skaniPresetFlag(String p) {
    def flags = [fast: '--fast', medium: '', slow: '--slow']
    def key   = (p ?: 'medium').toLowerCase()
    if (!flags.containsKey(key)) {
        error "--skani_preset must be one of: fast, medium, slow"
    }
    flags[key]
}

// ════════════════════════════════════════════════════════════════════════════
// SKETCH PROCESSES (one genome → one cached sketch; storeDir = sketch cache)
// ════════════════════════════════════════════════════════════════════════════

process MASH_SKETCH {
    tag   "${genome.name}"
    label 'mash'
    cpus   2
    memory '4 GB'

    storeDir "${params.sketch_cache}/mash/k${params.mash_kmer}_s${params.mash_sketch_size}"

    input:
        tuple val(group_name), path(genome)

    output:
        tuple val(group_name), path("${genome.baseName}.msh")

    script:
    """
    mash sketch -k ${params.mash_kmer} -s ${params.mash_sketch_size} \\
        -o ${genome.baseName} ${genome}
    """

    stub:
    """
    touch ${genome.baseName}.msh
    """
}

// skani sketching is fast, so genomes are batched --skani_sketch_chunk per job
// (instead of one SLURM job per genome) to cut scheduler churn. Output filenames
// are declared deterministically (one <genome>.sketch per input) so storeDir
// still caches each chunk and reuses it across runs / --compare levels.
process SKANI_SKETCH {
    tag   "${group_name} [${genomes.size()} genomes]"
    label 'skani'
    cpus   4
    memory '8 GB'

    storeDir "${params.sketch_cache}/skani/${params.skani_preset}_af${params.skani_min_af}"

    input:
        tuple val(group_name), path(genomes), val(sketch_names)

    output:
        tuple val(group_name), path(sketch_names)

    script:
    def preset = skaniPresetFlag(params.skani_preset)
    def cflag  = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    """
    # skani sketch errors out if its -o dir already exists; clear any partial
    # output left behind by an interrupted/retried run before re-sketching.
    rm -rf sk_out
    printf '%s\\n' ${genomes} > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sk_out
    mv sk_out/*.sketch .
    """

    stub:
    """
    for s in ${sketch_names.join(' ')}; do touch \$s; done
    """
}

process SOURMASH_SKETCH {
    tag   "${genome.name}"
    label 'sourmash'
    cpus   2
    memory '4 GB'

    storeDir "${params.sketch_cache}/sourmash/k${params.sourmash_kmer}_scaled${params.sourmash_scaled}"

    input:
        tuple val(group_name), path(genome)

    output:
        tuple val(group_name), path("${genome.name}.sig.zip")

    script:
    """
    sourmash sketch dna -p k=${params.sourmash_kmer},scaled=${params.sourmash_scaled} \\
        --name "${genome.name}" ${genome} -o ${genome.name}.sig.zip
    """

    stub:
    """
    touch ${genome.name}.sig.zip
    """
}

// ════════════════════════════════════════════════════════════════════════════
// COMPARE PROCESSES — each emits normalized  query<TAB>reference<TAB>ANI
// ════════════════════════════════════════════════════════════════════════════

// ── skani triangle (sparse) over cached sketches ─────────────────────────────
process SKANI_COMPARE {
    tag   "${group_name} n=${sketches.size()}"
    label 'skani'

    cpus   { sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8 }
    memory { sketches.size() > 500 ? '64 GB' : sketches.size() > 200 ? '32 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(sketches)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    ls *.sketch > sketch_list.txt
    skani triangle --sparse -l sketch_list.txt \\
        --min-af ${params.skani_min_af} -t ${task.cpus} -o skani_raw.tsv

    # skani cols: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query [names]
    awk -F'\\t' 'NR>1 {
        nr=split(\$1,b,"/"); r=b[nr]; sub(/\\.sketch\$/,"",r);
        nq=split(\$2,a,"/"); q=a[nq]; sub(/\\.sketch\$/,"",q);
        if (q!=r) print q"\\t"r"\\t"\$3
    }' skani_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}

// ── mash dist over cached sketches ───────────────────────────────────────────
process MASH_COMPARE {
    tag   "${group_name} n=${sketches.size()}"
    label 'mash'

    cpus   { sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8 }
    memory { sketches.size() > 500 ? '32 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(sketches)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    ls *.msh > sketch_list.txt
    mash paste combined -l sketch_list.txt
    mash dist -p ${task.cpus} combined.msh combined.msh > mash_raw.tsv

    # mash cols: ref  query  distance  p-value  shared-hashes  ->  ANI = 100*(1-dist)
    awk -F'\\t' '{
        nr=split(\$1,b,"/"); r=b[nr];
        nq=split(\$2,a,"/"); q=a[nq];
        if (q!=r) printf "%s\\t%s\\t%.4f\\n", q, r, (1-\$3)*100
    }' mash_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}

// ── sourmash compare (ANI matrix) over cached signatures ─────────────────────
process SOURMASH_COMPARE {
    tag   "${group_name} n=${sigs.size()}"
    label 'sourmash'

    cpus   { sigs.size() > 500 ? 32 : sigs.size() > 200 ? 16 : 8 }
    memory { sigs.size() > 500 ? '128 GB' : sigs.size() > 200 ? '64 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(sigs)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    ls *.sig.zip > sig_list.txt
    sourmash compare --ani --containment --quiet \\
        --from-file sig_list.txt --csv cmp.csv -k ${params.sourmash_kmer}

    python3 ${projectDir}/bin/sourmash_matrix_to_long.py --input cmp.csv --output ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}

// ── fastANI (query list × ref list) ──────────────────────────────────────────
// Used both for full-group single jobs and per-component prefiltered jobs.
// query/ref staged into separate subdirs so identical files don't collide.
// group_size drives CPU/memory scaling.
process FASTANI_COMPARE {
    tag   "${group_name} [${batch_tag}] n=${group_size}"
    label 'fastani'

    cpus   { group_size > 500 ? 64 : group_size > 200 ? 24 : 8 }
    memory { group_size > 500 ? '128 GB' : group_size > 200 ? '48 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name),
              path(query_genomes, stageAs: 'query/*'),
              path(ref_genomes,   stageAs: 'ref/*'),
              val(batch_tag), val(group_size)

    output:
        tuple val(group_name), path("${group_name}.${batch_tag}.ani.tsv")

    script:
    """
    ls query/* > query_list.txt
    ls ref/*   > ref_list.txt
    fastANI \\
        --ql query_list.txt \\
        --rl ref_list.txt \\
        -o ${group_name}.${batch_tag}.ani.tsv \\
        --fragLen ${params.fastani_fraglen} \\
        -k ${params.fastani_kmer} \\
        -t ${task.cpus}
    """

    stub:
    """
    printf 'query\\treference\\t99.0\\t100\\t100\\n' > ${group_name}.${batch_tag}.ani.tsv
    """
}

// ── mash prefilter: paste sketches, all-vs-all dist (cheap) ───────────────────
process MASH_PREFILTER {
    tag   "${group_name} n=${sketches.size()}"
    label 'mash'

    cpus   { sketches.size() > 500 ? 32 : sketches.size() > 200 ? 16 : 8 }
    memory { sketches.size() > 500 ? '32 GB' : '16 GB' }

    input:
        tuple val(group_name), path(sketches)

    output:
        tuple val(group_name), path("${group_name}.mash_dist.tsv")

    script:
    """
    ls *.msh > sketch_list.txt
    mash paste combined -l sketch_list.txt
    mash dist -p ${task.cpus} combined.msh combined.msh > ${group_name}.mash_dist.tsv
    """

    stub:
    """
    touch ${group_name}.mash_dist.tsv
    """
}

// ── connected components from mash dist (python, in the report env) ──────────
process MASH_COMPONENTS {
    tag   "${group_name}"
    label 'report'

    input:
        tuple val(group_name), path(mash_dist)

    output:
        tuple val(group_name), path("${group_name}.components.tsv")

    script:
    """
    python3 ${projectDir}/bin/mash_components.py --input ${mash_dist} \\
        --ani ${params.prefilter_ani} --min-size ${params.min_group_size} \\
        > ${group_name}.components.tsv
    """

    stub:
    """
    touch ${group_name}.components.tsv
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MERGE + REPORT (shared tail for all methods)
// ════════════════════════════════════════════════════════════════════════════

// Concatenate one-or-many partial TSVs into the group's full ANI table.
process MERGE_ANI {
    tag   "${group_name}"
    label 'report'

    storeDir "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}"

    input:
        tuple val(group_name), path(part_tsvs)

    output:
        tuple val(group_name), path("${group_name}.ani.tsv")

    script:
    """
    cat ${part_tsvs} > ${group_name}.ani.tsv
    """

    stub:
    """
    cat ${part_tsvs} > ${group_name}.ani.tsv
    """
}

process REPORT_ANI {
    tag   "${group_name}"
    label 'report'

    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}" }, mode: 'copy'

    input:
        tuple val(group_name), path(ani_tsv), path(names_tsv)

    output:
        path("${group_name}_ANI_report.txt")
        path("${group_name}_genome_names.tsv")

    script:
    """
    cp ${names_tsv} ${group_name}_genome_names.tsv
    python3 ${projectDir}/bin/report_ani.py \\
        --input    ${ani_tsv} \\
        --names    ${group_name}_genome_names.tsv \\
        --group    "${group_name}" \\
        --level    "${params.compare}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --outlier-threshold ${params.ani_outlier_threshold} \\
        --output   ${group_name}_ANI_report.txt
    """

    stub:
    """
    printf '=== ANI Report: ${group_name} (stub) ===\\n' > ${group_name}_ANI_report.txt
    cp ${names_tsv} ${group_name}_genome_names.tsv
    """
}

// Merge every group's pair table + names lookup into one labeled, queryable
// table (CSV + SQLite) so pairwise identities can be looked up across the
// whole run without re-parsing per-group files.
process COMBINE_ANI_TABLE {
    label 'report'

    publishDir "${params.outdir}/${params.ani_method}/${params.compare}", mode: 'copy'

    input:
        path ani_tsvs,   stageAs: 'ani_in/*'
        path names_tsvs, stageAs: 'names_in/*'

    output:
        path("all_pairs.csv")
        path("ani.db")

    script:
    """
    python3 ${projectDir}/bin/combine_ani_table.py \\
        --ani-dir   ani_in \\
        --names-dir names_in \\
        --compare-level "${params.compare}" \\
        --csv-output all_pairs.csv \\
        --db-output  ani.db
    """

    stub:
    """
    touch all_pairs.csv ani.db
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Validate params ───────────────────────────────────────────────────────
    def validRanks = ['PHYLUM','SUBPHYLUM','CLASS','SUBCLASS','ORDER','FAMILY','GENUS']
    def compareRank = (params.compare as String).toUpperCase()
    if (!(compareRank in validRanks)) {
        error "--compare must be one of: ${validRanks.join(', ')}"
    }

    def method = (params.ani_method as String).toLowerCase()
    if (!(method in ['skani','mash','sourmash','fastani'])) {
        error "--ani_method must be one of: skani, mash, sourmash, fastani"
    }

    def nameStyle = (params.genome_name_style as String).toLowerCase()
    if (!(nameStyle in ['species', 'asmid'])) {
        error "--genome_name_style must be 'species' or 'asmid'"
    }
    log.info "ANI method: ${method}"
    log.info "Genome filename style: ${nameStyle} (suffix: ${params.genome_suffix})"
    if (method == 'fastani' && (params.fastani_prefilter as boolean)) {
        log.info "fastANI mash-prefilter cascade ENABLED (prefilter_ani=${params.prefilter_ani}%)"
    }

    // ── Taxonomy filter ───────────────────────────────────────────────────────
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be RANK:VALUE, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { _row -> true }
    }

    // ── Build sample channel ──────────────────────────────────────────────────
    // Emits: tuple(group_key, locustag, species, genome_path, genus, strain, asmid)
    def sample_ch = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def locustag  = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def species   = row.SPECIES?.trim() ?: ''
            def genus     = row.GENUS?.trim() ?: ''
            def strain    = row.STRAIN?.trim() ?: ''
            def asmid     = row.ASMID?.trim() ?: ''
            def groupKey  = row[compareRank]?.trim() ?: ''
            def stem      = nameStyle == 'asmid'
                                ? row.ASMID?.trim()
                                : SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def genome    = file("${params.genome_dir}/${stem}${params.genome_suffix}", glob: false)

            if (!groupKey) {
                log.warn "Skipping ${locustag}: empty ${compareRank} field"
                return null
            }
            if (!genome.exists()) {
                log.warn "Skipping ${locustag}: genome not found at ${genome}"
                return null
            }
            tuple(groupKey, locustag, species, genome, genus, strain, asmid)
        }
        .filter { item -> item != null }

    // ── Group by taxonomic rank ───────────────────────────────────────────────
    // n_test limits *groups* (applied after groupTuple so --n_test 3 = 3 groups).
    def grouped_ch = sample_ch
        .map { groupKey, locustag, species, genome, genus, strain, asmid ->
            tuple(groupKey, tuple(locustag, species, genome, genus, strain, asmid))
        }
        .groupTuple()
        .filter { _gname, members -> members.size() >= params.min_group_size as int }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // ── Per-group genome list + names TSV ─────────────────────────────────────
    // names_tsv (filename\tspecies) is the authoritative genome universe for the
    // report — genomes with no surviving pair still show up there as outliers.
    // Only genomes whose file actually exists on disk are listed: the names file
    // (and the genome list fed to the compare step) must never reference a genome
    // that is absent from the input folder.
    def prepared_ch = grouped_ch
        .map { group_name, members ->
            def present = members.findAll { m -> m[2].exists() }
            def missing = members.size() - present.size()
            if (missing > 0) {
                log.warn "${group_name}: omitting ${missing} genome(s) missing from ${params.genome_dir}"
            }
            tuple(group_name, present)
        }
        .filter { _gname, present -> present.size() >= params.min_group_size as int }
        .map { group_name, present ->
            def genomes   = present.collect { m -> m[2] }
            // names TSV columns: filename, asmid, genus, species, strain
            // (m = locustag, species, genome, genus, strain, asmid)
            def namesText = present.collect { m ->
                "${m[2].getName()}\t${m[5]}\t${m[3]}\t${m[1]}\t${m[4]}"
            }.join('\n')
            def nameFile  = file("${workflow.workDir}/names_${group_name}.tsv")
            nameFile.text = "filename\tasmid\tgenus\tspecies\tstrain\n" + namesText + "\n"
            tuple(group_name, genomes, nameFile)
        }

    def prepared_split = prepared_ch.multiMap { group_name, genomes, nameFile ->
        ani:   tuple(group_name, genomes)
        names: tuple(group_name, nameFile)
    }
    def genome_ch = prepared_split.ani         // tuple(group_name, [genome files])
    def names_map = prepared_split.names        // tuple(group_name, names_file)

    // Flattened (group, genome) for per-genome sketching.
    def genome_flat = genome_ch.flatMap { group_name, genomes ->
        genomes.collect { g -> tuple(group_name, g) }
    }

    // ── Dispatch by method → produce ani_tsv_ch: tuple(group_name, full_ani_tsv)
    def ani_tsv_ch

    if (method == 'skani') {
        // Batch genomes into chunks of skani_sketch_chunk per sketch job.
        def chunk = Math.max(1, params.skani_sketch_chunk as int)
        def skani_sketch_in = genome_ch.flatMap { gn, genomes ->
            genomes.collate(chunk).collect { sub ->
                tuple(gn, sub, sub.collect { g -> "${g.name}.sketch" })
            }
        }
        ani_tsv_ch = SKANI_SKETCH(skani_sketch_in)
            .groupTuple()
            .map { gn, sketch_lists -> tuple(gn, sketch_lists.flatten()) }
            | SKANI_COMPARE

    } else if (method == 'mash') {
        ani_tsv_ch = MASH_SKETCH(genome_flat)
            .groupTuple()
            | MASH_COMPARE

    } else if (method == 'sourmash') {
        ani_tsv_ch = SOURMASH_SKETCH(genome_flat)
            .groupTuple()
            | SOURMASH_COMPARE

    } else { // fastani
        def part_ch

        if (params.fastani_prefilter as boolean) {
            // (b) mash pre-cluster, then fastANI all-vs-all within each component.
            def grouped_msk = MASH_SKETCH(genome_flat).groupTuple()
            def comps_ch    = MASH_PREFILTER(grouped_msk) | MASH_COMPONENTS

            // Explode components.tsv (comp_id <TAB> genome_filename) into rows
            // keyed by (group, filename) so we can re-attach the genome file path.
            def comp_rows = comps_ch
                .splitCsv(elem: 1, sep: '\t')
                .map { gn, row -> tuple(tuple(gn, row[1]), tuple(gn, row[0])) }

            def genome_lut = genome_ch.flatMap { gn, genomes ->
                genomes.collect { g -> tuple(tuple(gn, g.name), g) }
            }

            def per_comp = comp_rows
                .combine(genome_lut, by: 0)                       // (key, (gn,comp), gpath)
                .map { _key, gncomp, gpath -> tuple(gncomp[0], gncomp[1], gpath) }
                .groupTuple(by: [0, 1])                           // (gn, comp, [gpaths])
                .filter { _gn, _comp, gs -> gs.size() >= params.min_group_size as int }
                .map { gn, comp, gs -> tuple(gn, gs, gs, "c${comp}", gs.size()) }

            part_ch = FASTANI_COMPARE(per_comp)

        } else {
            def batchSize = params.ani_batch_size as int

            if (batchSize > 0) {
                // Upper-triangle batch pairs; carry total group size for CPU scaling.
                def batch_pairs_ch = genome_ch.flatMap { group_name, genomes ->
                    def N      = genomes.size()
                    def nBatch = (N + batchSize - 1).intdiv(batchSize)
                    def batches = (0..<nBatch).collect { i ->
                        def end = (i + 1) * batchSize < N ? (i + 1) * batchSize : N
                        genomes.subList(i * batchSize, end)
                    }
                    (0..<nBatch).collectMany { i ->
                        (i..<nBatch).collect { j ->
                            tuple(group_name, batches[i], batches[j], "b${i}_b${j}", N)
                        }
                    }
                }
                part_ch = FASTANI_COMPARE(batch_pairs_ch)

            } else {
                // Single job per group; query == ref == all genomes.
                def single_ch = genome_ch.map { group_name, genomes ->
                    tuple(group_name, genomes, genomes, "full", genomes.size())
                }
                part_ch = FASTANI_COMPARE(single_ch)
            }
        }

        // Merge per-group partials (1 for single job, many for batch/prefilter).
        ani_tsv_ch = part_ch
            .groupTuple()
            .map { gn, tsv_list -> tuple(gn, tsv_list.flatten()) }
            | MERGE_ANI
    }

    // ── Generate reports ──────────────────────────────────────────────────────
    REPORT_ANI(ani_tsv_ch.join(names_map, by: 0))

    // ── Combine every group's pairs + labels into one queryable table ────────
    COMBINE_ANI_TABLE(
        ani_tsv_ch.map { _gn, tsv -> tsv }.collect(),
        names_map.map  { _gn, tsv -> tsv }.collect()
    )
}
