#!/usr/bin/env python3
"""Shared sanitizers for building samples.csv from NCBI accession + taxonomy tables.

These functions are the single source of truth for cleaning SPECIES / STRAIN /
ASMID values.  clean_strain() is intentionally kept in lock-step with the Groovy
helper nextflow/lib/SampleUtils.groovy::cleanStrain so that the strain content in
samples.csv and the filesystem tags built by makeSampleTag() agree.

Aggressive strain normalisation (per project decision):
  - first synonym token only (split on ';' and '=')
  - ':'  -> ' '
  - '*'  -> removed at start/end, '-' between two words
  - '#'  -> '-'
  - ','  -> '-'   (also avoids needing CSV quoting)
  - '+'  -> removed (rare mating-type marker, e.g. 'S mat+')
  - surrounding quotes stripped, whitespace / dashes collapsed
Spaces and the safe set [A-Za-z0-9 ._-] are preserved so strain stays readable.
"""

import re

# species epithets are frequently quoted by NCBI, e.g.  sp. 'richmondensis'
_QUOTES = re.compile(r"""['"]""")
_NOM_INVAL = re.compile(r"\s*\(nom\.\s*inval\.\)", re.IGNORECASE)
_WS = re.compile(r"\s+")
_MULTIDASH = re.compile(r"-{2,}")


def clean_species(raw):
    """Canonicalise a species name (binomial or with embedded strain).

    Strips surrounding/embedded quote characters and '(nom. inval.)' notes and
    collapses whitespace.  Does NOT remove the strain that some NCBI species
    names carry (that distinction is SPECIES vs SPECIES_IN upstream)."""
    s = (raw or "").strip()
    s = _NOM_INVAL.sub("", s)
    s = _QUOTES.sub("", s)
    # upstream add_taxonomy converts NCBI brackets [Candida] -> _Candida_ ; undo
    # only when the underscores wrap a whole token (genus), so legitimate inner
    # underscores like 'GC_Crypt_4' or 'PMI_857' are preserved.
    s = re.sub(r"(^|(?<=\s))_([^_\s]+)_(?=\s|$)", r"\1\2", s)   # _Candida_ -> Candida
    s = re.sub(r"(?<=\w)#(?=\w)", "-", s)     # Mi166#008 -> Mi166-008
    s = s.replace("#", "")                     # 'str. #12' -> 'str. 12'
    s = _WS.sub(" ", s).strip()
    return s


def clean_strain(raw):
    """Aggressively normalise a strain label to clean, near-filesystem-safe text.

    Mirrors nextflow/lib/SampleUtils.groovy::cleanStrain for the shared rules
    (';' first token, ':'->space, '*' handling) and adds the project's
    aggressive normalisation for the remaining special characters."""
    s = (raw or "").strip()
    s = _QUOTES.sub("", s)
    s = s.split(";")[0]            # synonymous strains separated by ';'
    s = s.split("=")[0]            # synonymous strains separated by '='
    s = s.strip()
    s = s.replace(":", " ")
    # asterisk handling (kept identical in spirit to SampleUtils.cleanStrain)
    s = re.sub(r"^\s*\*+", "", s)        # leading '*'  -> removed
    s = re.sub(r"\*+\s*$", "", s)        # trailing '*' -> removed
    s = re.sub(r"\s*\*+\s*", "-", s)     # internal '*' -> '-'
    # remaining aggressive normalisation
    s = s.replace("#", "-")
    s = s.replace(",", "-")
    s = s.replace("+", "")
    s = _WS.sub(" ", s).strip()
    s = _MULTIDASH.sub("-", s)
    s = s.strip(" -.")
    return s


# trailing genome-file extension artefacts seen in ASM_NAME, e.g.
#   GA10_6-scaf.final.scaffolds.fasta   ->  GA10_6-scaf
#   A6_trimmed_clc_assembly_flt.fa      ->  A6_trimmed_clc_assembly_flt
# Legitimate version tokens such as _v1.0 are preserved.
_ASM_EXT = re.compile(
    r"(?:\.(?:final|scaffolds|renamed|scaf|genomic))*"
    r"\.(?:fa|fasta|fna)$",
    re.IGNORECASE,
)


def clean_asmid(raw):
    """Strip genome file-extension tails from an assembly identifier."""
    s = (raw or "").strip()
    s = _ASM_EXT.sub("", s)
    return s


def backfill_strain(strain, species_in, species):
    """If strain is empty but the verbatim species name carries trailing strain
    tokens (SPECIES_IN beyond the binomial SPECIES), recover them.

    e.g. species_in='Thermochaetoides thermophila DSM 1495',
         species='Thermochaetoides thermophila'  ->  'DSM 1495'
    """
    if strain:
        return strain
    si = (species_in or "").strip()
    sp = (species or "").strip()
    if not (sp and si.startswith(sp) and len(si) > len(sp)):
        return strain
    remainder = si[len(sp):].strip(" -")
    # Only accept remainders that look like a culture-collection strain token:
    # 1-3 whitespace tokens AND containing at least one digit.  This rejects
    # hybrid/taxonomic remainders such as '-Leptosphaeria biglobosa ... group'.
    if remainder and any(c.isdigit() for c in remainder) and len(remainder.split()) <= 3:
        return remainder
    return strain


if __name__ == "__main__":
    # tiny self-test / demonstration
    for fn, val in [
        (clean_species, "Acidomyces sp. 'richmondensis'"),
        (clean_species, "Foo bar (nom. inval.)"),
        (clean_strain, "IIF2*SW-F2"),
        (clean_strain, "SR23;CBS 7157"),
        (clean_strain, "EXF-153=EXF-2781"),
        (clean_strain, "okayama7#130"),
        (clean_strain, "BRIP:72097"),
        (clean_strain, "S mat+"),
        (clean_strain, "CRUB 1588,7"),
        (clean_asmid, "GCA_051397355.1_GA10_6-scaf.final.scaffolds.fasta"),
        (clean_asmid, "GCA_000401675.1_A6_trimmed_clc_assembly_flt.fa"),
        (clean_asmid, "GCA_001592465.1_Acidomyces_richmondensis_BFW_v1.0"),
    ]:
        print(f"{fn.__name__}({val!r}) = {fn(val)!r}")
