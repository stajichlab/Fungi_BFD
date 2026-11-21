include { plainName; gzInflate; hashBucketForType } from '../../common/utils.nf'

// CAZY stores one file per genome per output type, hash-bucketed by LOCUSTAG,
// the same layout every other function/* type uses (todo/genome_stats_storage_reorg.md
// §A: "LOCUSTAG for ... all of results/function/*", width 2).
//
// It was the last type on the old one-subdirectory-per-genome layout keyed by
// meta.id (Genus_species_strain). That key is NOT stable: it is rebuilt from
// SPECIES and STRAIN, so any taxonomic rename moved a genome to a fresh
// directory and orphaned the old one, which the merge glob then kept reading.
// By 2026-11-21 that had produced 24 genomes present under two names at once
// -- a silent double count in cazy.overview.parquet on any merge_all build --
// plus 8 directories for genomes no longer in samples.csv. LOCUSTAG does not
// change when taxonomy does, so the rewrite now lands on the same path and the
// failure cannot recur. Migrated by scripts/one-off/migrate_cazy_to_hash_buckets.py.
//
// storeDir is the cache here, not just a published copy: Nextflow skips this
// process when the declared output files already exist under it. The migration
// and this module therefore have to be applied together -- either one alone
// leaves storeDir looking at an empty path and re-runs dbCAN for every genome.
process RUN_CAZY {
    tag        "${meta.locustag}"
    label      'cazy'
    storeDir   { "${params.outdir}/cazy/${hashBucketForType('cazy', meta.locustag as String)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.overview.tsv.gz"),   emit: overview
        path("${meta.locustag}.cazymes.tsv.gz"),    emit: cazymes
        path("${meta.locustag}.substrates.tsv.gz"), emit: substrates

    script:
    def proteins_plain = plainName(proteins)
    """
    ${gzInflate(proteins)}
    module load dbcanlight
    mkdir -p ${meta.locustag}
    dbcanlight search -i ${proteins_plain} -m cazyme -o ${meta.locustag} -t ${task.cpus}
    dbcanlight search -i ${proteins_plain} -m sub    -o ${meta.locustag} -t ${task.cpus}
    dbcanlight conclude ${meta.locustag}
    pigz -f ${meta.locustag}/cazymes.tsv ${meta.locustag}/substrates.tsv ${meta.locustag}/overview.tsv
    mv ${meta.locustag}/overview.tsv.gz   ${meta.locustag}.overview.tsv.gz
    mv ${meta.locustag}/cazymes.tsv.gz    ${meta.locustag}.cazymes.tsv.gz
    mv ${meta.locustag}/substrates.tsv.gz ${meta.locustag}.substrates.tsv.gz
    """

    stub:
    """
    printf 'Gene_ID\\tEC\\tcazyme_fam\\tsub_fam\\tdiamond_fam\\tSubstrate\\t#ofTools\\n' | gzip > ${meta.locustag}.overview.tsv.gz
    printf 'HMM_Profile\\tProfile_Length\\tGene_ID\\tGene_Length\\tEvalue\\tProfile_Start\\tProfile_End\\tGene_Start\\tGene_End\\tCoverage\\n' | gzip > ${meta.locustag}.cazymes.tsv.gz
    printf 'Gene_ID\\tSubstrate\\n' | gzip > ${meta.locustag}.substrates.tsv.gz
    """
}
