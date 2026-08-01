#!/usr/bin/env python3
"""Extract a taxon-scoped subset of the master BFD DuckDB into a standalone DB file.

Per T-014 §D.2: the master BFD DuckDB (built by build_BFD_duckDB.sh) is always
full/unscoped now -- taxon-specific slices are no longer built at merge time,
they're extracted post-hoc from the master DB here. `species` is the
authoritative key table: filtered first by --taxon, then every other table is
restricted by INNER JOIN against that filtered key set (LOCUSTAG, or the
derived species_prefix for LOCUSTAG-prefixed protein/domain-level tables) --
never by independently re-evaluating the taxon predicate per table, which
would risk orphaned rows if a table's taxonomy columns ever drifted from
`species`'s. This is purely a performance/convenience slice, not an access
boundary (confirmed 2026-08-01) -- treat it as a correctness-focused tool, not
a security-focused one.

Usage
-----
    python3 scripts/extract_bfd_taxonomic_subset.py \\
        --taxon GENUS:Malassezia --master db/BFD.duckdb -o db/BFD.Malassezia.duckdb

    # List valid ranks / preview a match count without writing anything:
    python3 scripts/extract_bfd_taxonomic_subset.py \\
        --taxon GENUS:Malassezia --master db/BFD.duckdb -o /dev/null --dry-run
"""

import argparse
import os
import sys

import duckdb

# species/asm_stats/gene_info carry LOCUSTAG directly. Every other table is
# keyed by species_prefix (derived at DB-build time via
# string_split(<id_col>, '_')[1] -- see build_BFD_duckDB.sh).
LOCUSTAG_KEYED_TABLES = ["asm_stats", "gene_info"]
SPECIES_PREFIX_KEYED_TABLES = [
    "gene_proteins", "gene_transcripts", "gene_exons", "gene_CDS",
    "gene_introns", "gene_trna", "gene_intergenic_distances",
    "pfam", "cazy", "cazy_overview", "merops", "signalp", "targetp",
    "wolfpsort", "predgpi", "tmhmm", "idp", "idp_summary",
    "aa_frequency", "codon_frequency",
]

# Index DDL mirrored from build_BFD_duckDB.sh -- kept in sync manually since
# DuckDB has no "describe indexes on read-only attached DB, replay them" API
# that would let this be derived automatically without also re-deriving
# constraints (UNIQUE) that CREATE TABLE AS SELECT doesn't preserve.
INDEX_DDL = {
    "species": [
        "CREATE UNIQUE INDEX idx_species_locustag ON species(LOCUSTAG)",
        "CREATE UNIQUE INDEX idx_species_asm ON species(ASMID)",
        'CREATE INDEX idx_species_taxonomy ON species(PHYLUM, CLASS, "ORDER", FAMILY, GENUS)',
    ],
    "asm_stats": [
        "CREATE UNIQUE INDEX idx_asm_locustag ON asm_stats(LOCUSTAG)",
        "CREATE UNIQUE INDEX idx_asm_asmid ON asm_stats(ASMID)",
    ],
    "gene_info": [
        "CREATE UNIQUE INDEX idx_gene_info_gene_id ON gene_info(gene_id)",
        "CREATE INDEX idx_gene_info_locustag ON gene_info(LOCUSTAG)",
    ],
    "gene_proteins": [
        "CREATE UNIQUE INDEX idx_gene_prot_protein_id ON gene_proteins(protein_id)",
        "CREATE INDEX idx_gene_prot_locustag ON gene_proteins(species_prefix)",
    ],
    "gene_transcripts": [
        "CREATE UNIQUE INDEX idx_tx_transcript_id ON gene_transcripts(transcript_id)",
        "CREATE INDEX idx_tx_locustag ON gene_transcripts(species_prefix)",
    ],
    "gene_exons": [
        "CREATE UNIQUE INDEX idx_exon_exon_id ON gene_exons(exon_id)",
        "CREATE INDEX idx_exon_locustag ON gene_exons(species_prefix)",
    ],
    "gene_CDS": [
        "CREATE UNIQUE INDEX idx_cds_cds_id ON gene_CDS(cds_id)",
        "CREATE INDEX idx_cds_locustag ON gene_CDS(species_prefix)",
    ],
    "gene_introns": [
        "CREATE UNIQUE INDEX idx_intron_intron_id ON gene_introns(intron_id)",
        "CREATE INDEX idx_intron_locustag ON gene_introns(species_prefix)",
    ],
    "gene_trna": [
        "CREATE UNIQUE INDEX idx_trna_gene_id ON gene_trna(gene_id)",
        "CREATE INDEX idx_trna_locustag ON gene_trna(species_prefix)",
    ],
    "gene_intergenic_distances": [
        "CREATE INDEX idx_intergenic_locustag ON gene_intergenic_distances(species_prefix)",
    ],
    "pfam": [
        "CREATE UNIQUE INDEX idx_pfam_domain_hit ON pfam(protein_id, pfam_id, domain_num)",
        "CREATE INDEX idx_pfam_locustag ON pfam(species_prefix)",
        "CREATE INDEX idx_pfam_pfam_id ON pfam(pfam_id)",
    ],
    "cazy": [
        "CREATE UNIQUE INDEX idx_cazy_hit ON cazy(protein_id, HMM_id, s_start)",
        "CREATE INDEX idx_cazy_locustag ON cazy(species_prefix)",
        "CREATE INDEX idx_cazy_hmm_id ON cazy(HMM_id)",
    ],
    "cazy_overview": [
        "CREATE UNIQUE INDEX idx_cazyov_protein_id ON cazy_overview(protein_id)",
        "CREATE INDEX idx_cazyov_locustag ON cazy_overview(species_prefix)",
    ],
    "merops": [
        "CREATE INDEX idx_merops_locustag ON merops(species_prefix)",
        "CREATE INDEX idx_merops_merops_id ON merops(merops_id)",
    ],
    "signalp": [
        "CREATE INDEX idx_signalp_locustag ON signalp(species_prefix)",
        "CREATE INDEX idx_signalp_protein ON signalp(protein_id)",
    ],
    "targetp": [
        "CREATE INDEX idx_targetp_locustag ON targetp(species_prefix)",
        "CREATE INDEX idx_targetp_protein ON targetp(protein_id)",
    ],
    "wolfpsort": [
        "CREATE INDEX idx_wolfpsort_locustag ON wolfpsort(species_prefix)",
        "CREATE INDEX idx_wolfpsort_protein ON wolfpsort(protein_id)",
    ],
    "predgpi": [
        "CREATE INDEX idx_predgpi_locustag ON predgpi(species_prefix)",
        "CREATE INDEX idx_predgpi_protein ON predgpi(protein_id)",
    ],
    "tmhmm": [
        "CREATE INDEX idx_tmhmm_locustag ON tmhmm(species_prefix)",
        "CREATE INDEX idx_tmhmm_protein ON tmhmm(protein_id)",
    ],
    "idp_summary": [
        "CREATE UNIQUE INDEX idx_idpsum_protein ON idp_summary(protein_id)",
        "CREATE INDEX idx_idpsum_locustag ON idp_summary(species_prefix)",
    ],
    "idp": [
        "CREATE UNIQUE INDEX idx_idp_region ON idp(protein_id, IDP_start)",
        "CREATE INDEX idx_idp_locustag ON idp(species_prefix)",
    ],
    "aa_frequency": [
        "CREATE INDEX idx_aa_locustag ON aa_frequency(species_prefix)",
    ],
    "codon_frequency": [
        "CREATE INDEX idx_codon_locustag ON codon_frequency(species_prefix)",
    ],
}

VALID_RANKS = {
    "PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER", "FAMILY",
    "GENUS", "SPECIES",
}


def parse_taxon(taxon: str):
    if ":" not in taxon:
        sys.exit(f"ERROR: --taxon must be RANK:VALUE, e.g. --taxon GENUS:Malassezia (got {taxon!r})")
    rank, _, value = taxon.partition(":")
    rank = rank.strip().upper()
    value = value.strip()
    if not rank or not value:
        sys.exit(f"ERROR: --taxon must be RANK:VALUE, e.g. --taxon GENUS:Malassezia (got {taxon!r})")
    if rank not in VALID_RANKS:
        sys.exit(f"ERROR: unknown rank '{rank}' -- must be one of {sorted(VALID_RANKS)}")
    return rank, value


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxon", required=True, help="RANK:VALUE, e.g. GENUS:Malassezia (RANK matches a species table column)")
    p.add_argument("--master", default="db/BFD.duckdb", help="path to the master BFD DuckDB [db/BFD.duckdb]")
    p.add_argument("-o", "--outfile", required=True, help="path for the extracted subset DuckDB (overwritten if it exists)")
    p.add_argument("--dry-run", action="store_true", help="report the match count and exit without writing an output DB")
    args = p.parse_args()

    rank, value = parse_taxon(args.taxon)

    if not os.path.exists(args.master):
        sys.exit(f"ERROR: master DB not found: {args.master}")

    con = duckdb.connect(":memory:")
    con.execute(f"ATTACH '{args.master}' AS master (READ_ONLY)")

    match_count = con.execute(
        f'SELECT count(*) FROM master.species WHERE "{rank}" = ?', [value]
    ).fetchone()[0]
    print(f"Taxon filter {rank}={value!r}: {match_count} genomes match in master.species")
    if match_count == 0:
        sys.exit(f"ERROR: no genomes matched {rank}={value!r} -- check spelling/case (exact match, same as --taxon elsewhere in this pipeline)")

    if args.dry_run:
        print("(--dry-run: not writing an output DB)")
        return

    if os.path.exists(args.outfile):
        os.remove(args.outfile)
    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    con.execute(f"ATTACH '{args.outfile}' AS out")

    con.execute(f'CREATE TABLE out.species AS SELECT * FROM master.species WHERE "{rank}" = ?', [value])
    for stmt in INDEX_DDL["species"]:
        con.execute(f"{stmt.replace('ON species', 'ON out.species', 1)}")

    for table in LOCUSTAG_KEYED_TABLES:
        con.execute(
            f"CREATE TABLE out.{table} AS "
            f"SELECT m.* FROM master.{table} m "
            f"INNER JOIN out.species s ON m.LOCUSTAG = s.LOCUSTAG"
        )
        for stmt in INDEX_DDL[table]:
            con.execute(stmt.replace(f"ON {table}", f"ON out.{table}", 1))

    for table in SPECIES_PREFIX_KEYED_TABLES:
        con.execute(
            f"CREATE TABLE out.{table} AS "
            f"SELECT m.* FROM master.{table} m "
            f"INNER JOIN out.species s ON m.species_prefix = s.LOCUSTAG"
        )
        for stmt in INDEX_DDL[table]:
            con.execute(stmt.replace(f"ON {table}", f"ON out.{table}", 1))

    # Recreate the two analytical views verbatim from master, rather than
    # duplicating their SQL here -- avoids drift if build_BFD_duckDB.sh's view
    # definitions ever change. Table names inside the SQL are unqualified, so
    # they resolve against out's own just-built tables once executed with out
    # as the connection's default catalog.
    view_rows = con.execute(
        "SELECT view_name, sql FROM duckdb_views() WHERE internal = false AND database_name = 'master'"
    ).fetchall()
    con.execute("USE out")
    for view_name, view_sql in view_rows:
        con.execute(view_sql)
    con.execute("USE master")

    # Post-extract row-count assertion (robust-analysis convention: fail
    # loudly on unexpected data, don't trust the JOIN silently). Every
    # LOCUSTAG-keyed and species_prefix-keyed table must have zero rows whose
    # key is absent from the just-built species table -- trivially true given
    # INNER JOIN, asserted anyway as a guard against a future refactor
    # accidentally loosening a JOIN to LEFT JOIN.
    errors = []
    for table in LOCUSTAG_KEYED_TABLES:
        orphans = con.execute(
            f"SELECT count(*) FROM out.{table} t "
            f"WHERE t.LOCUSTAG NOT IN (SELECT LOCUSTAG FROM out.species)"
        ).fetchone()[0]
        if orphans:
            errors.append(f"{table}: {orphans} rows with LOCUSTAG not in species")
    for table in SPECIES_PREFIX_KEYED_TABLES:
        orphans = con.execute(
            f"SELECT count(*) FROM out.{table} t "
            f"WHERE t.species_prefix NOT IN (SELECT LOCUSTAG FROM out.species)"
        ).fetchone()[0]
        if orphans:
            errors.append(f"{table}: {orphans} rows with species_prefix not in species")
    if errors:
        con.close()
        os.remove(args.outfile)
        sys.exit("ERROR: orphaned rows found after extraction (output DB removed):\n  " + "\n  ".join(errors))

    print(f"\nExtracted {rank}={value!r} to {args.outfile}:")
    for table in ["species"] + LOCUSTAG_KEYED_TABLES + SPECIES_PREFIX_KEYED_TABLES:
        n = con.execute(f"SELECT count(*) FROM out.{table}").fetchone()[0]
        print(f"  {table}: {n} rows")

    con.close()


if __name__ == "__main__":
    main()
