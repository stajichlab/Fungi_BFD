-- BUILD_SWISSPROT_ANNOT: swissprot_subject.parquet, one row per SwissProt entry.
-- sp_id is the integer key the hit table stores instead of accession/entry/name
-- strings (those took 38.9 GB of the 146 GB swissprot.parquet). sp_id is assigned
-- by accession order and is stable only within one SwissProt release; the hit
-- tables are rebuilt whenever the release (swissprot_dat) changes.
-- Placeholder {annot}: the swissprot_annot.csv written by parse_uniprot_sprot.py.
COPY (
  SELECT row_number() OVER (ORDER BY accession)::INTEGER AS sp_id,
         accession, entry_name,
         split_part(entry_name, '_', 2) AS org_code,
         length::INTEGER AS s_len,
         protein_name, gene_name, gene_synonyms, organism, ncbi_taxid::BIGINT AS ncbi_taxid,
         go_terms, ec_numbers, interpro, pfam
  FROM read_csv({annot}, header = true, auto_detect = true, sample_size = -1)
  ORDER BY sp_id
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
