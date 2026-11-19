-- REF_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs
-- bfd_traits.sql: FungalTraits ecology columns and SUBPHYLUM for the local BFD.duckdb.
-- Run after bfd_compat_views.sql (needs the `species` view) and BEFORE bfd_extra_views.sql
-- (v_genome_labels joins v_genome_traits). Line 1 is replaced with the real config/refs path by
-- workflows/build_bfd_db.sh. Idempotent. The three ref_* tables are small (about 11,000 rows in total) and are stored
-- in the database so a dumped copy does not depend on files elsewhere.
-- Rules are those of bin/atlas/build_labels.py (growth_form, growth_form_source, in_fungaltraits, yeast_origin).
-- FungalTraits is a GENUS-level label, not a genome trait. Cite: Polme et al. 2020, Fungal Diversity 105:1-16 (verify).

CREATE OR REPLACE TABLE ref_fungaltraits AS
SELECT GENUS,
       Growth_form_template AS growth_form_ft,
       primary_lifestyle,
       Secondary_lifestyle AS secondary_lifestyle,
       Plant_pathogenic_capacity_template AS plant_pathogenic_capacity,
       Decay_substrate_template AS decay_substrate,
       Animal_biotrophic_capacity_template AS animal_biotrophic_capacity
FROM (SELECT *, row_number() OVER (PARTITION BY GENUS ORDER BY jrk_template) AS rn
      FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs/fungaltraits_1.2.utf8.csv', header=true, all_varchar=true))
WHERE rn = 1 AND GENUS IS NOT NULL;
COMMENT ON TABLE ref_fungaltraits IS 'FungalTraits 1.2, one row per GENUS (first row by jrk_template when a genus is listed twice). Source: config/refs/fungaltraits_1.2.utf8.csv.';

CREATE OR REPLACE TABLE ref_growth_form_override AS
SELECT GENUS, max(nullif(growth_form, '')) AS growth_form_ovr
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs/growth_form_overrides.tsv', header=true, delim='\t', all_varchar=true)
WHERE GENUS IS NOT NULL
GROUP BY GENUS;
COMMENT ON TABLE ref_growth_form_override IS 'Growth form for genera that FungalTraits lacks (config/refs/growth_form_overrides.tsv). An empty growth_form is ignored.';

CREATE OR REPLACE TABLE ref_class_subphylum AS
SELECT CLASS, min(SUBPHYLUM) AS SUBPHYLUM, min(source) AS source, min(note) AS note
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs/class_subphylum.tsv', header=true, delim='\t', all_varchar=true)
WHERE CLASS IS NOT NULL
GROUP BY CLASS;
COMMENT ON TABLE ref_class_subphylum IS 'CLASS to SUBPHYLUM (config/refs/class_subphylum.tsv). source: taxonomy_py, atlas, or proposed (not yet checked). Classes not listed have no subphylum here.';

CREATE OR REPLACE TABLE ref_phylum_subphylum AS
SELECT PHYLUM, min(SUBPHYLUM) AS SUBPHYLUM, min(source) AS source, min(note) AS note
FROM read_csv('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/refs/phylum_subphylum.tsv', header=true, delim='\t', all_varchar=true)
WHERE PHYLUM IS NOT NULL
GROUP BY PHYLUM;
COMMENT ON TABLE ref_phylum_subphylum IS 'PHYLUM to SUBPHYLUM fallback for genomes whose CLASS is missing or not in ref_class_subphylum (config/refs/phylum_subphylum.tsv).';

CREATE OR REPLACE VIEW v_genome_traits AS
SELECT t.LOCUSTAG, t.growth_form, t.growth_form_source, t.primary_lifestyle, t.secondary_lifestyle,
       t.plant_pathogenic_capacity, t.decay_substrate, t.animal_biotrophic_capacity, t.in_fungaltraits,
       CASE WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes') THEN 'Saccharomycotina'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Schizosaccharomycetes','Taphrinomycetes','Pneumocystomycetes','Neolectomycetes') THEN 'Taphrinomycotina'
            WHEN t.growth_form = 'yeast' AND t.CLASS = 'Tremellomycetes' THEN 'Tremellomycetes'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Microbotryomycetes','Cystobasidiomycetes','Agaricostilbomycetes','Mixiomycetes','Spiculogloeomycetes','Classiculomycetes','Tritirachiomycetes') THEN 'Pucciniomycotina-yeast'
            WHEN t.growth_form = 'yeast' AND t.CLASS IN ('Malasseziomycetes','Ustilaginomycetes','Exobasidiomycetes','Moniliellomycetes') THEN 'Ustilaginomycotina-yeast'
            WHEN t.growth_form = 'yeast' THEN 'yeast-other'
            WHEN t.growth_form = 'dimorphic_yeast' AND t.CLASS IN ('Eurotiomycetes','Dothideomycetes','Sordariomycetes','Leotiomycetes','Lecanoromycetes','Pezizomycetes','Orbiliomycetes','Arthoniomycetes','Xylonomycetes','Geoglossomycetes','Lichinomycetes','Laboulbeniomycetes','Coniocybomycetes','Candelariomycetes') THEN 'Pezizomycotina-dimorphic'
            WHEN t.growth_form = 'dimorphic_yeast' AND t.CLASS IN ('Malasseziomycetes','Ustilaginomycetes','Exobasidiomycetes','Moniliellomycetes') THEN 'Ustilaginomycotina-dimorphic'
            WHEN t.growth_form = 'dimorphic_yeast' THEN 'dimorphic-other'
            WHEN t.growth_form = 'filamentous_mycelium' THEN 'filamentous'
            ELSE 'other' END AS yeast_origin,
       t.SUBPHYLUM, t.subphylum_source
FROM (
  SELECT s.LOCUSTAG, s.CLASS,
         coalesce(o.growth_form_ovr, f.growth_form_ft) AS growth_form,
         CASE WHEN o.growth_form_ovr IS NOT NULL THEN 'override'
              WHEN f.growth_form_ft IS NOT NULL THEN 'fungaltraits' END AS growth_form_source,
         f.primary_lifestyle, f.secondary_lifestyle, f.plant_pathogenic_capacity, f.decay_substrate,
         f.animal_biotrophic_capacity,
         (f.GENUS IS NOT NULL) AS in_fungaltraits,
         coalesce(s.SUBPHYLUM, c.SUBPHYLUM, ph.SUBPHYLUM) AS SUBPHYLUM,
         CASE WHEN s.SUBPHYLUM IS NOT NULL THEN 'species'
              WHEN c.SUBPHYLUM IS NOT NULL THEN 'class_table'
              WHEN ph.SUBPHYLUM IS NOT NULL THEN 'phylum_table' END AS subphylum_source
  FROM species s
  LEFT JOIN ref_fungaltraits f ON s.GENUS = f.GENUS
  LEFT JOIN ref_growth_form_override o ON s.GENUS = o.GENUS
  LEFT JOIN ref_class_subphylum c ON s.CLASS = c.CLASS
  LEFT JOIN ref_phylum_subphylum ph ON s.PHYLUM = ph.PHYLUM
) t;
COMMENT ON VIEW v_genome_traits IS 'One row per genome (key LOCUSTAG): FungalTraits growth form and lifestyle by GENUS (in_fungaltraits false = genus not in FungalTraits 1.2), yeast_origin, and SUBPHYLUM (species.SUBPHYLUM if filled, else the class table, else the phylum table; NULL when none). Same names and rules as atlas genome_labels.parquet.';
