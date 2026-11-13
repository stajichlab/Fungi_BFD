-- MERGE_SWISSPROT (one genome bucket per task): DIAMOND/blastp -outfmt 6 tables
-- (18 columns, see merge_swissprot.py) -> swissprot_hits/bucket=XX/data.parquet.
-- Design: analysis/swissprot_design/REVIEW.md (2026-11-06).
--   * subject strings are replaced by sp_id (join swissprot_subject)
--   * rk = rank of the hit within its protein (bitscore desc, evalue, sp_id)
--   * pident / coverages as x10 integers; qcov = DIAMOND qcovhsp; scov =
--     (send - sstart + 1) / slen. The old query_cov / hit_cov divided the
--     alignment length (gaps included) by the sequence length and exceeded 1.0
--     in ~10-17% of rows.
--   * evalue as -log10 (REAL); evalue = 0 (underflow) -> 400.0
--   * nident / positive / mismatch / gapopen stay in the per-genome raw files
-- Placeholders: {files}, {out}, {subject} (swissprot_subject.parquet path, quoted).
-- Fails if a subject accession is missing from swissprot_subject.
CREATE TEMP TABLE raw AS
  SELECT * FROM read_csv({files}, delim = '\t', header = false, quote = '', escape = '',
    auto_detect = false,
    columns = {'qseqid':'VARCHAR','sseqid':'VARCHAR','pident':'DOUBLE','positive':'INTEGER',
               'nident':'INTEGER','length':'INTEGER','mismatch':'INTEGER','gapopen':'INTEGER',
               'qstart':'INTEGER','qend':'INTEGER','sstart':'INTEGER','send':'INTEGER',
               'evalue':'DOUBLE','bitscore':'DOUBLE','qcovhsp':'DOUBLE','qlen':'INTEGER',
               'slen':'INTEGER','stitle':'VARCHAR'});
CREATE TEMP TABLE j AS
  SELECT r.*, s.sp_id
  FROM raw r LEFT JOIN read_parquet({subject}) s
    ON s.accession = CASE WHEN contains(r.sseqid, '|') THEN split_part(r.sseqid, '|', 2) ELSE r.sseqid END;
SELECT CASE WHEN count(*) > 0 THEN error('swissprot_hits: ' || count(*) || ' hits have no swissprot_subject entry, e.g. ' || any_value(sseqid)) END
FROM j WHERE sp_id IS NULL;
COPY (
  SELECT qseqid AS protein_id,
         split_part(qseqid, '_', 1) AS species_prefix,
         sp_id,
         row_number() OVER (PARTITION BY qseqid ORDER BY bitscore DESC, evalue, sp_id)::UTINYINT AS rk,
         round(pident * 10)::USMALLINT AS pident_x10,
         round(qcovhsp * 10)::USMALLINT AS qcov_x10,
         round((send - sstart + 1) * 1000.0 / slen)::USMALLINT AS scov_x10,
         length AS aln_length,
         qstart AS q_start, qend AS q_end, sstart AS s_start, send AS s_end,
         qlen AS q_len,
         bitscore::REAL AS bitscore,
         (CASE WHEN evalue = 0 THEN 400.0 ELSE -log10(evalue) END)::REAL AS neg_log10_evalue
  FROM j
  ORDER BY protein_id, rk
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);
