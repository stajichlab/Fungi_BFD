# Reference files for the DB build

Loaded by `sql/bfd_traits.sql` into `ref_fungaltraits`, `ref_growth_form_override` and `ref_class_subphylum`.

| File | Source | Notes |
|---|---|---|
| `fungaltraits_1.2.utf8.csv` | FungalTraits 1.2, `DeltaGain_Fungi/refs/lifestyle/FungalTraits_1.2.csv` (latin-1), converted to UTF-8 with the same text-mode conversion as `bin/atlas/build_labels.py` | Cite as Polme et al. 2020, Fungal Diversity 105:1-16 (verify). One row per genus except a few duplicates. |
| `growth_form_overrides.tsv` | Copy of `atlas_data/labels/growth_form_overrides.tsv` (2026-11-18) | Genera that FungalTraits lacks; `source` is `class_rule`. |
| `class_subphylum.tsv` | `taxonomy_py` rows from `bin/summary/taxonomy.py`; `atlas` rows from the class lists in `bin/atlas/build_labels.py`; `proposed` rows added 2026-11-18 and 2026-11-19; `ncbi` rows (8 classes) copied from the NCBI `species.SUBPHYLUM` values of 2026-11-18 | `proposed` rows are not yet checked against a cited classification. Classes in phyla without subphyla, and classes that could not be placed, are left out, so `SUBPHYLUM` stays NULL for them. |

| `phylum_subphylum.tsv` | Added 2026-11-18 at the user's request | Fallback when the CLASS is missing or not in `class_subphylum.tsv`: Kickxellomycota, Mortierellomycota, Glomeromycota, Basidiobolomycota, Blastocladiomycota, Microsporidia, Rozellomycota, Aphelidiomycota (the last five added 2026-11-19). Phyla with several subphyla (Mucoromycota, Zoopagomycota) have no phylum entry. |

`Pneumocystomycetes` is the spelling in the genome data. `taxonomy.py` has `Pneumocystidomycetes`, which matched no genome.
