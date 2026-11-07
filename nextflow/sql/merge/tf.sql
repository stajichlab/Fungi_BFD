-- MERGE_TF: Pfam <LOCUSTAG>.domtblout.gz files -> tf_inventory.parquet
-- Same rows, columns and types as bin/merge_tf_domains.py followed by
-- read_csv_auto(sample_size=-1). Needs the TF domain list:
--   --param tf_domains=<path to assets/fungal_tf_pfam_domains.csv>
-- (columns pfam_acc,pfam_name,tf_family). A domain hit is kept when its Pfam
-- accession without the version (PF00172.25 -> PF00172) is in the list.
-- The list is loaded into a MAP variable (not a JOIN) so the output keeps the
-- input row order. If an accession is listed twice the last entry wins, as in
-- the Python dict. The list is read as CSV (quoted fields allowed), like
-- csv.DictReader on this file. tf_family values that contain a comma are quoted
-- (PF11955 WOPR; before 2026-11-06 it was unquoted and was cut at the comma).
SET VARIABLE tf_map = (
  WITH d AS (
    SELECT pfam_acc AS acc, pfam_name AS nm, tf_family AS fam, row_number() OVER () AS rn
    FROM read_csv('{tf_domains}', header = true, quote = '"', auto_detect = false,
                  columns = {'pfam_acc': 'VARCHAR', 'pfam_name': 'VARCHAR', 'tf_family': 'VARCHAR'})
    WHERE pfam_acc IS NOT NULL AND pfam_acc <> ''
  ), u AS (
    SELECT acc, arg_max(nm, rn) AS nm, arg_max(fam, rn) AS fam FROM d GROUP BY acc
  )
  SELECT map(list(acc), list({'name': nm, 'fam': fam})) FROM u
);
-- Cheap pre-filter: keep only lines that contain a listed accession as a
-- whole word followed by '.' or whitespace. This is a superset of the hits
-- the exact field-5 test below keeps; it avoids splitting every line.
SET VARIABLE tf_re = (
  SELECT '(^|\s)(' || string_agg(DISTINCT regexp_escape(k), '|') || ')(\.|\s|$)'
  FROM (SELECT unnest(map_keys(getvariable('tf_map'))) AS k)
);
COPY (
  WITH l AS (
    SELECT regexp_split_to_array(regexp_replace(line, '^\s+|\s+$', '', 'g'), '\s+') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = E'\x01',
                  quote = '', escape = '', auto_detect = false)
    WHERE line IS NOT NULL AND line NOT LIKE '#%'
      AND regexp_matches(line, getvariable('tf_re'))
  ), h AS (
    SELECT f, split_part(f[5], '.', 1) AS acc,
           getvariable('tf_map')[split_part(f[5], '.', 1)][1] AS tf
    FROM l
    WHERE len(f) >= 21
  )
  SELECT NULLIF(split_part(f[1], '_', 1), '') AS species_prefix,
         NULLIF(f[1], '') AS protein_id,
         NULLIF(acc, '') AS pfam_acc,
         NULLIF(tf.name, '') AS pfam_name,
         NULLIF(tf.fam, '') AS tf_family,
         f[13]::DOUBLE AS domain_i_evalue,
         f[14]::DOUBLE AS domain_score
  FROM h
  WHERE tf IS NOT NULL
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
