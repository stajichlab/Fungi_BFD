include { hashBucketForType } from '../../common/utils.nf'

// tRNA gene calling, run standalone in BFD instead of reusing funannotate's
// cached tRNAscan-SE output -- keeps BFD's functional-annotation stack
// self-contained and independently re-runnable. Mirrors funannotate predict's
// own tRNA handling (funannotate/library.py runtRNAscan()+validate_tRNA()):
//   1. tRNAscan-SE on the genome FASTA (containerized, not `module load`, per
//      the "less HPCC-dependent" direction).
//   2. NCBI length rule (50-150 bp) + drop Pseudo/Sup/Undet calls, done by
//      bin/trnascan2gff3.pl (ported verbatim from funannotate's aux_scripts/,
//      so the exact same filtering logic is reused, not reimplemented).
//   3. bedtools intersect -v against the existing protein-coding gene GFF3 --
//      drops any tRNA gene model that overlaps a predicted gene, same as
//      funannotate's validate_tRNA(). No assembly-gap track here (BFD's
//      cleaned assemblies are the same ones already gap-filtered upstream by
//      GENOME_CLEAN), so only the gene-overlap intersection is applied.
process RUN_TRNASCAN {
    tag        "${meta.id}"
    label      'trnascan'
    storeDir   { "${params.outdir}/trnascan/${hashBucketForType('trnascan', meta.id)}" }

    input:
        tuple val(meta), path(gff3), path(genome)

    output:
        path("${meta.id}.trnascan.no-overlaps.gff3.gz"), emit: gff3

    script:
    """
    module load singularity
    singularity exec ${params.trnascan_sif} tRNAscan-SE -o ${meta.id}.tRNAscan.out --thread ${task.cpus} ${genome}

    awk -F'\\t' 'NR<=3 {print; next} {len=(\$3>\$4)?\$3-\$4:\$4-\$3; if (len>=50 && len<=150) print}' \\
        ${meta.id}.tRNAscan.out > ${meta.id}.tRNAscan.len-filtered.out

    singularity exec ${params.trnascan_sif} perl ${projectDir}/bin/trnascan2gff3.pl \\
        --input ${meta.id}.tRNAscan.len-filtered.out > ${meta.id}.trnascan.raw.gff3

    module load bedtools/2.30.0
    bedtools sort -i ${meta.id}.trnascan.raw.gff3   > ${meta.id}.trnascan.sorted.gff3
    bedtools sort -i ${gff3}                        > ${meta.id}.genes.sorted.gff3
    bedtools intersect -sorted -v \\
        -a ${meta.id}.trnascan.sorted.gff3 \\
        -b ${meta.id}.genes.sorted.gff3 > ${meta.id}.trnascan.no-overlaps.gff3

    pigz ${meta.id}.trnascan.no-overlaps.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.id}.trnascan.no-overlaps.gff3.gz
    """
}
