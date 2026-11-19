#!/usr/bin/env python3
"""Compare the trait columns of v_genome_labels in a built BFD.duckdb with the atlas genome_labels.parquet.

    python3 bin/check_traits_vs_atlas.py --db PATH.duckdb --atlas PATH.parquet [--ft-csv PATH]

Compares growth_form and primary_lifestyle per LOCUSTAG. A genus listed twice in FungalTraits can legitimately
differ (the atlas keeps an arbitrary row, the database keeps the first by jrk_template): those genomes are
reported but do not fail the check. Also prints the genomes with NULL SUBPHYLUM, by CLASS.
Genomes present on only one side are counted and reported but do not fail the check (the atlas is a snapshot).
Exit 1 when any other genome differs in growth_form or primary_lifestyle.
"""
import argparse
import sys
from pathlib import Path

import duckdb

DEFAULT_FT = Path(__file__).resolve().parents[1] / "config" / "refs" / "fungaltraits_1.2.utf8.csv"
COLS = ["growth_form", "primary_lifestyle"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--ft-csv", default=str(DEFAULT_FT))
    a = ap.parse_args(argv)

    con = duckdb.connect(a.db, read_only=True)
    con.execute(f"CREATE TEMP VIEW atlas AS SELECT * FROM read_parquet('{a.atlas}')")
    con.execute(f"""CREATE TEMP VIEW dup_genera AS
                    SELECT GENUS FROM read_csv('{a.ft_csv}', header=true, all_varchar=true)
                    GROUP BY GENUS HAVING count(*) > 1""")
    n_db = con.execute("SELECT count(*) FROM v_genome_labels").fetchone()[0]
    n_atlas = con.execute("SELECT count(*) FROM atlas").fetchone()[0]
    only_db = con.execute("SELECT count(*) FROM v_genome_labels l ANTI JOIN atlas a USING (LOCUSTAG)").fetchone()[0]
    only_atlas = con.execute("SELECT count(*) FROM atlas a ANTI JOIN v_genome_labels l USING (LOCUSTAG)").fetchone()[0]
    print(f"genomes: database {n_db:,}, atlas {n_atlas:,}; only in database {only_db:,}, only in atlas {only_atlas:,}")

    bad = False   # genome-set differences are reported only: the atlas is a snapshot and the live tables change
    for c in COLS:
        rows = con.execute(f"""
            SELECT l.LOCUSTAG, l.GENUS, l.{c}, a.{c}, (l.GENUS IN (SELECT GENUS FROM dup_genera)) AS dup
            FROM v_genome_labels l JOIN atlas a USING (LOCUSTAG)
            WHERE l.{c} IS DISTINCT FROM a.{c} ORDER BY 1""").fetchall()
        hard = [r for r in rows if not r[4]]
        dup = [r for r in rows if r[4]]
        print(f"{c}: {len(hard):,} differ; {len(dup):,} more differ in genera listed twice in FungalTraits (duplicate rows)")
        for r in hard[:20]:
            print(f"  {r[0]}\t{r[1]}\tdatabase={r[2]}\tatlas={r[3]}")
        bad = bad or bool(hard)

    nul = con.execute("SELECT count(*) FROM v_genome_labels WHERE SUBPHYLUM IS NULL").fetchone()[0]
    print(f"SUBPHYLUM NULL: {nul:,} of {n_db:,} genomes. By class (top 25):")
    for phy, cls, n in con.execute("""SELECT PHYLUM, CLASS, count(*) FROM v_genome_labels WHERE SUBPHYLUM IS NULL
                                      GROUP BY ALL ORDER BY 3 DESC, 2 LIMIT 25""").fetchall():
        print(f"  {n:>6,}\t{phy}\t{cls}")
    print("FAIL" if bad else "OK")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
