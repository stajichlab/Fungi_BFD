#!/usr/bin/env python3
"""Export the merged antiSMASH parquet tables as the gzip TSV tables sent to the EBC-KT agents.

Inputs : --inputs input_files.txt (the file list given to the merge); --pq DIR with bgc_cluster, bgc_cluster_gene, bgc_mibig_hit parquet (merge_antismash.py output);
         --samples samples.csv; --runjson GLOB of <LOCUSTAG>.antismash.run.json;
         --linked-missing TSV (output of check_antismash_complete.py --out: genomes with linked inputs and no regions table; column LOCUSTAG).
Run status: ok (in the merge), no_run_yet (linked-missing list, or finished after the merge), not_annotated (every other genome without a run).
Outputs (in --out): regions_v5.tsv.gz, genomes_v5.tsv.gz, mibig_hits_v5.tsv.gz,
         flagged_clusters_v5.tsv.gz, flagged_summary_v5.tsv, export_counts.txt
Column choices follow regions_with_mibig_v3 / regions_v4 where the data allow. Differences:
  genome = LOCUSTAG (v3 used a short strain label). best_hit_rank = position by similarity
  (1 = best); best_antismash_rank = antiSMASH rank of the same hit. No pubchem/formula/inchikey columns.
"""
import argparse, os, subprocess, sys

SQL = r"""
SET threads = {threads};
CREATE VIEW cl AS SELECT * FROM read_parquet('{pq}/bgc_cluster.parquet');
CREATE VIEW cg AS SELECT * FROM read_parquet('{pq}/bgc_cluster_gene.parquet');
CREATE VIEW hh AS SELECT * FROM read_parquet('{pq}/bgc_mibig_hit.parquet');
CREATE TABLE sm AS SELECT * FROM read_csv('{samples}', header=true, all_varchar=true);
CREATE TABLE snap AS SELECT DISTINCT regexp_extract(f, '([^/]+)\.antismash\.regions\.tsv\.gz$', 1) AS genome
  FROM read_csv('{inputs}', header=false, columns={{f:'VARCHAR'}}) WHERE f LIKE '%.antismash.regions.tsv.gz';
CREATE TABLE rj_all AS SELECT * FROM read_json('{runjson}', format='newline_delimited', union_by_name=true,
   columns={{genome:'VARCHAR', antismash_version:'VARCHAR', mibig_data_version:'VARCHAR', mode:'VARCHAR', n_records:'BIGINT', n_regions:'BIGINT'}});
CREATE TABLE rj AS SELECT * FROM rj_all WHERE genome IN (SELECT genome FROM snap);
CREATE TABLE ms AS {missing_sql};

-- hits ordered by similarity (ties keep antiSMASH rank order); sim_rank 1 = best
CREATE TABLE hr AS
  SELECT *, row_number() OVER (PARTITION BY locustag, record, region ORDER BY similarity DESC, rank ASC) AS sim_rank,
            count(*) FILTER (WHERE similarity = max_sim) OVER (PARTITION BY locustag, record, region) AS n_tied_best
  FROM (SELECT *, max(similarity) OVER (PARTITION BY locustag, record, region) AS max_sim FROM hh);

CREATE TABLE bg AS
  SELECT locustag, record, region, string_agg(gene_id, ',' ORDER BY start, gene_id) AS biosynthetic_genes
  FROM cg WHERE gene_kind = 'biosynthetic' GROUP BY ALL;

CREATE TABLE reg AS
SELECT c.locustag AS genome, s.SPECIES AS species, s.STRAIN AS strain, s.NCBI_TAXONID AS ncbi_taxid, s.ASMID AS assembly,
       c.record, c.region, c.start, c."end", c.products, c.contig_edge, c.n_cds,
       CASE WHEN b.similarity >= 75 THEN 'known' WHEN b.similarity >= 1 THEN 'weak' ELSE 'none' END AS status,
       CASE WHEN b.similarity >= 75 AND b.n_entry_proteins >= 3 AND b.gene_hits::DOUBLE / b.n_entry_proteins >= 0.5 AND b.core_gene_hits >= 1 THEN 'known_multigene'
            WHEN b.similarity >= 75 AND b.n_entry_proteins < 3 THEN 'known_single'
            WHEN b.similarity >= 75 THEN 'known_other'
            WHEN b.similarity >= 1 AND b.n_entry_proteins >= 3 AND b.gene_hits::DOUBLE / b.n_entry_proteins >= 0.5 AND b.core_gene_hits >= 1 THEN 'partial_multigene'
            WHEN b.similarity >= 1 THEN 'weak' ELSE 'none' END AS evidence_class,
       b.mibig_accession, b.compounds AS mibig_compounds, b.hit_cluster_type AS mibig_class,
       CASE WHEN b.similarity >= 1 THEN 1 END AS best_hit_rank, b.rank AS best_antismash_rank, b.n_tied_best,
       t.top3_hits, b.similarity, b.core_gene_hits, b.gene_hits, b.n_entry_proteins,
       round(b.gene_hits::DOUBLE / nullif(b.n_entry_proteins, 0), 3) AS entry_fraction_hit, b.synteny_score,
       c.length, c.core_start, c.core_end, c.candidate_kinds, c.categories, c.n_biosynthetic, c.n_additional,
       g.biosynthetic_genes
FROM cl c
LEFT JOIN sm s ON s.LOCUSTAG = c.locustag
LEFT JOIN hr b ON b.locustag=c.locustag AND b.record=c.record AND b.region=c.region AND b.sim_rank=1 AND b.similarity >= 1
LEFT JOIN (SELECT locustag, record, region, string_agg(mibig_accession||'|'||coalesce(compounds,'')||'|'||similarity, ';' ORDER BY sim_rank) AS top3_hits
           FROM hr WHERE sim_rank <= 3 GROUP BY ALL) t ON t.locustag=c.locustag AND t.record=c.record AND t.region=c.region
LEFT JOIN bg g ON g.locustag=c.locustag AND g.record=c.record AND g.region=c.region;

CREATE TABLE rc AS SELECT locustag, count(*) AS n_regions_table FROM cl GROUP BY 1;
CREATE TABLE gen AS
SELECT s.LOCUSTAG AS genome, s.ASMID AS assembly, s.NCBI_TAXONID AS ncbi_taxid, s.SPECIES_IN AS species_in, s.SPECIES AS species, s.STRAIN AS strain,
       s.BIOPROJECT AS bioproject, s.PHYLUM AS phylum, s.SUBPHYLUM AS subphylum, s.CLASS AS class, s.SUBCLASS AS subclass, s."ORDER" AS "order",
       s.FAMILY AS family, s.GENUS AS genus, s.TRANSL_TABLE AS transl_table,
       CASE WHEN r.genome IS NOT NULL THEN 'ok'
            WHEN ra.genome IS NOT NULL OR m.LOCUSTAG IS NOT NULL THEN 'no_run_yet'
            ELSE 'not_annotated' END AS run_status,
       CASE WHEN r.genome IS NOT NULL THEN NULL
            WHEN ra.genome IS NOT NULL THEN 'run finished after the merge snapshot; not in these tables'
            WHEN m.LOCUSTAG IS NOT NULL THEN 'inputs linked, no antiSMASH regions table (failed or not yet run); BFD contig_count=' || coalesce(m.contig_count, 'NA')
            ELSE 'no gene models / assembly input linked in BFD; antiSMASH cannot run' END AS status_note,
       r.antismash_version, r.mibig_data_version, r.mode AS antismash_mode, r.n_records, r.n_regions AS n_regions_runjson,
       coalesce(rc.n_regions_table, 0) AS n_regions_table
FROM sm s LEFT JOIN rj r ON r.genome = s.LOCUSTAG LEFT JOIN rj_all ra ON ra.genome = s.LOCUSTAG LEFT JOIN rc ON rc.locustag = s.LOCUSTAG LEFT JOIN ms m ON m.LOCUSTAG = s.LOCUSTAG
ORDER BY s.LOCUSTAG;

COPY (SELECT * FROM reg ORDER BY genome, record, region) TO '{out}/regions_v5.tsv.gz' (FORMAT csv, DELIMITER '\t', HEADER, COMPRESSION gzip);
COPY (SELECT * FROM gen) TO '{out}/genomes_v5.tsv.gz' (FORMAT csv, DELIMITER '\t', HEADER, COMPRESSION gzip);
COPY (SELECT h.locustag AS genome, h.record, h.region, h.sim_rank, h.rank AS antismash_rank, h.mibig_accession, h.compounds AS mibig_compounds,
             h.hit_cluster_type AS mibig_class, h.similarity, h.core_gene_hits, h.gene_hits, h.n_entry_proteins,
             round(h.gene_hits::DOUBLE / nullif(h.n_entry_proteins,0),3) AS entry_fraction_hit, h.synteny_score
      FROM hr h ORDER BY 1,2,3,4) TO '{out}/mibig_hits_v5.tsv.gz' (FORMAT csv, DELIMITER '\t', HEADER, COMPRESSION gzip);

CREATE TABLE flag(acc VARCHAR, label VARCHAR);
INSERT INTO flag VALUES ('BGC0000356','fumitremorgin/verruculogen'),('BGC0000420','roquefortine'),('BGC0000129','pyripyropene'),('BGC0000686','helvolic acid'),('BGC0001067','fumagillin');
CREATE TABLE fh AS
  SELECT f.acc AS flagged_accession, f.label, h.locustag AS genome, h.record, h.region, h.mibig_accession, h.sim_rank, h.similarity, h.n_entry_proteins, h.gene_hits, h.core_gene_hits
  FROM hr h JOIN flag f ON split_part(h.mibig_accession,'.',1) = f.acc;
COPY (SELECT x.*, s.SPECIES AS species, s.ASMID AS assembly FROM fh x LEFT JOIN sm s ON s.LOCUSTAG = x.genome ORDER BY 1,3,4,5) TO '{out}/flagged_clusters_v5.tsv.gz' (FORMAT csv, DELIMITER '\t', HEADER, COMPRESSION gzip);
-- denominators: genomes per species with a run (ok), genomes with a region whose top-5 hit is the flagged entry (any similarity / >=75)
COPY (
  SELECT f.acc AS flagged_accession, f.label, g.species, count(*) FILTER (WHERE g.run_status='ok') AS n_genomes_ok, count(*) AS n_genomes_all,
         count(DISTINCT g.genome) FILTER (WHERE x.n_any > 0) AS n_genomes_any_hit, count(DISTINCT g.genome) FILTER (WHERE x.n_ge75 > 0) AS n_genomes_hit_ge75,
         round(100.0 * count(DISTINCT g.genome) FILTER (WHERE x.n_ge75 > 0) / nullif(count(*) FILTER (WHERE g.run_status='ok'),0), 1) AS pct_ok_genomes_hit_ge75,
         round(100.0 * count(DISTINCT g.genome) FILTER (WHERE x.n_any > 0) / nullif(count(*) FILTER (WHERE g.run_status='ok'),0), 1) AS pct_ok_genomes_any_hit
  FROM flag f CROSS JOIN gen g
  LEFT JOIN (SELECT flagged_accession, genome, count(*) AS n_any, count(*) FILTER (WHERE similarity >= 75) AS n_ge75 FROM fh GROUP BY 1,2) x
         ON x.flagged_accession = f.acc AND x.genome = g.genome
  GROUP BY f.acc, f.label, g.species
  HAVING bool_or(g.genus IN ('Aspergillus','Penicillium','Talaromyces','Fusarium','Neosartorya')) OR count(DISTINCT g.genome) FILTER (WHERE x.n_any > 0) > 0
  ORDER BY 1, n_genomes_hit_ge75 DESC, 3) TO '{out}/flagged_summary_v5.tsv' (FORMAT csv, DELIMITER '\t', HEADER);

COPY (SELECT 'regions' AS tbl, count(*) AS n FROM reg UNION ALL SELECT 'genomes', count(*) FROM gen UNION ALL SELECT 'genomes_ok', count(*) FROM gen WHERE run_status='ok'
      UNION ALL SELECT 'genomes_without_run', count(*) FROM gen WHERE run_status<>'ok' UNION ALL SELECT 'genomes_no_run_yet', count(*) FROM gen WHERE run_status='no_run_yet' UNION ALL SELECT 'genomes_not_annotated', count(*) FROM gen WHERE run_status='not_annotated' UNION ALL SELECT 'hits', count(*) FROM hh
      UNION ALL SELECT 'regions_without_samples_row', count(*) FROM reg WHERE species IS NULL
      UNION ALL SELECT 'genomes_runjson_vs_table_regions_mismatch', count(*) FROM gen WHERE run_status='ok' AND n_regions_runjson <> n_regions_table
      UNION ALL SELECT 'runjson_without_samples_row', count(*) FROM rj WHERE genome NOT IN (SELECT LOCUSTAG FROM sm)) TO '{out}/export_counts.txt' (FORMAT csv, DELIMITER '\t', HEADER);
"""

def q(s): return s.replace("'", "''")

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pq", required=True); ap.add_argument("--samples", required=True)
    ap.add_argument("--runjson", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--inputs", required=True, help="list of merged input files (input_files.txt of the merge); genomes not in it are not counted as ok");
    ap.add_argument("--linked-missing", required=True); ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    miss = f"SELECT * FROM read_csv('{q(a.linked_missing)}', delim='\\t', header=true, all_varchar=true)"
    sql = SQL.format(pq=q(a.pq), samples=q(a.samples), runjson=q(a.runjson), inputs=q(a.inputs), out=q(a.out), threads=a.threads, missing_sql=miss)
    r = subprocess.run(["duckdb", "-bail"], input=sql, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("ERROR: duckdb failed:\n" + r.stderr[:3000])
    print(open(os.path.join(a.out, "export_counts.txt")).read())

if __name__ == "__main__":
    main()
