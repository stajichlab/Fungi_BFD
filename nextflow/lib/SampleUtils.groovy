/**
 * SampleUtils — shared helpers for constructing filesystem-safe sample identifiers
 * from samples.csv SPECIES and STRAIN columns.
 *
 * Placed in nextflow/lib/ so it is automatically on the Groovy classpath for all
 * workflow scripts in nextflow/*.nf.
 *
 * The Python equivalent is collect_busco_stats.py::build_basename_map() and
 * fix_low_trinity.py::species_key_from_row().  Keep them in sync.
 */
class SampleUtils {

    /**
     * Build a filesystem-safe "{species}_{strain}" tag from raw samples.csv values.
     *
     * Canonicalises:
     *   - strips leading/trailing whitespace and quote characters (', ") from both fields
     *   - takes only the first semicolon-delimited token of strain (some rows list
     *     multiple synonymous strains separated by ';')
     *   - replaces colons with spaces in strain (colons appear as ' colon ' separators)
     *   - collapses runs of whitespace, /, #, [, ], *, ?, {, } into single underscores
     *
     * Examples:
     *   makeSampleTag("Saccharomyces cerevisiae", "CBS 1171")    → "Saccharomyces_cerevisiae_CBS_1171"
     *   makeSampleTag("Aspergillus fumigatus", "Af293; CBS 101")  → "Aspergillus_fumigatus_Af293"
     *   makeSampleTag("Fusarium oxysporum", "")                   → "Fusarium_oxysporum"
     */
    static String makeSampleTag(String rawSpecies, String rawStrain) {
        def sp = (rawSpecies ?: '').trim().replaceAll(/['"]/, '')
        def st = (rawStrain  ?: '').trim()
                    .replaceAll(/['"]/, '')
                    .split(';')[0]
                    .trim()
                    .replace(':', ' ')
        return [sp, st].findAll { it }
                       .join('_')
                       .replaceAll(/[\s\/\#\[\]\*\?\{\}]+/, '_')
    }
}
