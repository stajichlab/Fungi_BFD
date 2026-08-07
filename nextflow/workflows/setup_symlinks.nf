//
// SETUP_SYMLINKS_ONLY — run just the input/ symlinking step from samples.csv,
// without triggering any of BFD.nf's downstream functional-annotation or
// genome-stats work.
//
// Shares the same rows_ch construction (taxon/asmid/suppress filters,
// meta.id/locustag mapping) and params (samples, genome_annotation, pep_dir,
// etc.) as BFD.nf, so `-profile BFD` config applies unchanged.
//
// Usage:
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD --pipeline setup_symlinks -resume
//

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { INPUT_SETUP } from '../subworkflows/local/INPUT_SETUP.nf'
include { taxonRowFilter; asmidRowFilter; loadSuppressSet; suppressRowFilter } from '../modules/common/utils.nf'
include { cleanStrain; makeSampleTag } from '../modules/common/utils.nf'

workflow SETUP_SYMLINKS_ONLY {

    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    def taxonFilter    = taxonRowFilter()
    def asmidFilter    = asmidRowFilter()
    def suppressFilter = suppressRowFilter(loadSuppressSet())

    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .filter(suppressFilter)
        .map { row ->
            [
                id      : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: ''),
                locustag: row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim(),
                species : row.SPECIES?.trim() ?: '',
                strain  : cleanStrain(row.STRAIN?.trim() ?: ''),
            ]
        }
        .take((params.n_test as int) > 0 ? (params.n_test as int) : -1)

    INPUT_SETUP(rows_ch)
}
