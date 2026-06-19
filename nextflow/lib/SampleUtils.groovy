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
     * Canonicalise a raw strain value into a human-readable strain name.
     *
     *   - strips leading/trailing whitespace and quote characters (', ")
     *   - takes only the first semicolon-delimited token (some rows list
     *     multiple synonymous strains separated by ';')
     *   - replaces colons with spaces (colons appear as ' colon ' separators)
     *   - handles asterisks: a '*' at the start or end of the strain is removed
     *     entirely; a '*' between two other words is replaced with '-'
     *
     * Examples:
     *   cleanStrain("Af293; CBS 101") → "Af293"
     *   cleanStrain("T-34*")          → "T-34"
     *   cleanStrain("ABC*DEF")        → "ABC-DEF"
     *   cleanStrain("ARSEF * 2860")   → "ARSEF-2860"
     */
    static String cleanStrain(String rawStrain) {
        return (rawStrain ?: '').trim()
                    .replaceAll(/['"]/, '')
                    .split(';')[0]
                    .trim()
                    .replace(':', ' ')
                    .replaceAll(/^\s*\*+/, '')      // '*' at the start of the strain -> removed
                    .replaceAll(/\*+\s*$/, '')      // '*' at the end of the strain   -> removed
                    .replaceAll(/\s*\*+\s*/, '-')   // '*' between two words          -> '-'
                    .trim()
    }

    /**
     * Build a filesystem-safe "{species}_{strain}" tag from raw samples.csv values.
     *
     * Canonicalises:
     *   - strips leading/trailing whitespace and quote characters (', ") from both fields
     *   - cleans the strain via {@link #cleanStrain} (first ';' token, colon→space,
     *     asterisk handling)
     *   - collapses runs of whitespace, /, #, [, ], ?, {, } into single underscores
     *
     * Examples:
     *   makeSampleTag("Saccharomyces cerevisiae", "CBS 1171")    → "Saccharomyces_cerevisiae_CBS_1171"
     *   makeSampleTag("Aspergillus fumigatus", "Af293; CBS 101")  → "Aspergillus_fumigatus_Af293"
     *   makeSampleTag("Fusarium oxysporum", "")                   → "Fusarium_oxysporum"
     *   makeSampleTag("Beauveria bassiana", "ARSEF 2860*")        → "Beauveria_bassiana_ARSEF_2860"
     */
    static String makeSampleTag(String rawSpecies, String rawStrain) {
        def sp = (rawSpecies ?: '').trim().replaceAll(/['"]/, '')
        def st = cleanStrain(rawStrain)
        return [sp, st].findAll { it }
                       .join('_')
                       .replaceAll(/[\s\/\#\[\]\?\{\}]+/, '_')
    }
}
