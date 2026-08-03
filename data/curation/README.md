# samples.csv — reproducible production & curation

`samples.csv` is built from the NCBI fungal genome download by
`scripts/create_samples_file.py`, using the shared sanitizer
`scripts/sample_sanitize.py` and the curation data in this folder.

## End-to-end chain

```
../../1KFG/2026/NCBI_fungi/   (pixi run make)
  datasets summary genome taxon fungi
    -> assembly_json_process.py -> ncbi_accessions.csv
    -> add_taxonomy.py (taxonkit) -> ncbi_accessions_taxonomy.csv
                       |
this repo:             v
  scripts/create_samples_file.py
    + scripts/sample_sanitize.py        (clean SPECIES / STRAIN / ASMID)
    + data/curation/exclude_asmids.txt  (hard removals)
    + data/curation/keep_dupes.csv      (protected multi-assembly isolates)
    + data/curation/overrides.csv       (optional per-ASMID field fixes)
    -> samples.csv
```

Run it:

```bash
python3 scripts/create_samples_file.py            # writes ./samples.csv
python3 scripts/audit_strain_changes.py           # writes the SPECIES/STRAIN audit
```

## Sanitization rules (scripts/sample_sanitize.py)

- **SPECIES**: strip quotes (`sp. 'x'` → `sp. x`), drop `(nom. inval.)`, undo NCBI
  bracket artifacts (`_Candida_` → `Candida`, from upstream `[Candida]`), strip `#`.
- **STRAIN** (aggressive): first `;`/`=` synonym token; `:`→space; `*` removed at
  ends / `-` between words; `#`→`-`; `,`→`-`; `+` dropped; collapse dashes/spaces.
  Kept byte-compatible with `nextflow/lib/SampleUtils.groovy::cleanStrain`.
- **STRAIN backfill**: if empty, recover a culture-collection token from the
  verbatim species name (only short, digit-bearing remainders).
- **ASMID**: strip genome-file extension tails (`.fa`, `.final.scaffolds.fasta`, …);
  preserve version tokens like `_v1.0`.

## Curation files

| File | Meaning |
|------|---------|
| `exclude_asmids.txt` | raw ASMIDs to drop (curated removals + raw-filename dedup). `<ASMID>\t<reason>`; `#` comments. |
| `preferred_asmids.txt` | ASMIDs that win their species+strain dedup group, overriding the default tie-breaker. `<ASMID>\t<reason>`; `#` comments. |
| `keep_dupes.csv` | species+strain isolates intentionally kept as >1 assembly (protected from the auto dedup tie-breaker). |
| `overrides.csv` | `ASMID,FIELD,VALUE` — force-correct a field value for a specific assembly (after sanitization; matched on cleaned ASMID). Changes content, not which assembly is kept. |
| `strain_species_audit.csv` | full `ASMID, ORIG/MODIFIED SPECIES, ORIG/MODIFIED STRAIN` record. |
| `strain_species_audit.changed.csv` | the subset that changed. |

**Pick the right file**: drop a genome → `exclude_asmids.txt`; keep a *specific*
assembly of a duplicated isolate → `preferred_asmids.txt`; keep *several* →
`keep_dupes.csv`; fix a *field value* on a kept row → `overrides.csv`.

### Dedup policy

Known duplicate assemblies are curated in `exclude_asmids.txt` (these were
case-by-case calls favouring the best-annotated/reference genome — neither
"prefer GCF" nor "newest" reproduces them, so they are recorded explicitly).
Any **new** species+strain collision not covered by curation is auto-resolved by
a default tie-breaker and logged. Winner precedence:

1. an ASMID listed in `preferred_asmids.txt` (force a specific assembly), else
2. RefSeq `GCF`, else
3. newest `GCA` accession.

Add the isolate to `keep_dupes.csv` to keep more than one.

## IMPORTANT: regenerate stale taxonomy before a production run

As of this writing the upstream `ncbi_accessions_taxonomy.csv` is **stale**:
it covers only ~8,061 of 22,412 accessions and every row is duplicated ~2.8×.
`create_samples_file.py` handles duplicates and missing rows defensively (logs
counts, falls back to the verbatim species name), but assemblies without
taxonomy get empty lineage and a `fungi` BUSCO default. Re-run the upstream
pipeline first:

```bash
cd ../../1KFG/2026/NCBI_fungi && pixi run make lib/ncbi_accessions_taxonomy.csv
```
