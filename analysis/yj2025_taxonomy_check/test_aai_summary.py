"""Fixture test for aai_summary.py (synthetic DIAMOND hits, runs in well under a second).

Run: python3 -m pytest analysis/yj2025_taxonomy_check/test_aai_summary.py
"""
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aai_summary as A  # noqa: E402


def test_aai_from_hits_filters_and_dedups():
    lines = [
        'p1\ts1\t80.0\t90\t100\t100\t1e-30\t200\n',   # counted
        'p1\ts9\t10.0\t90\t100\t100\t1e-30\t50\n',    # second row of p1: ignored
        'p2\ts2\t60.0\t40\t100\t100\t1e-30\t100\n',   # qcov 0.4 < 0.5: not counted
        'p3\ts3\t70.0\t60\t100\t100\t1e-3\t40\n',     # evalue too weak: not counted
        'p4\ts4\t50.0\t50\t100\t100\t1e-10\t90\n',    # qcov exactly 0.5: counted
    ]
    n, aai = A.aai_from_hits(lines)
    assert n == 2
    assert aai == 65.0


def test_cli_ranks_and_skips_self(tmp_path):
    panel = tmp_path / 'panel.tsv'
    cols = ['LOCUSTAG', 'role', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES', 'STRAIN', 'busco_complete', 'pep']
    rows = [['Q1', 'query', 'Ascomycota', '', '', '', '', 'Q sp.', 's', '90', 'x'],
            ['R1', 'rhodotorula', 'Basidiomycota', 'Microbotryomycetes', 'Sporidiobolales', 'Sporidiobolaceae',
             'Rhodotorula', 'Rhodotorula a', 's', '90', 'x'],
            ['S1', 'control', 'Ascomycota', 'Saccharomycetes', 'Saccharomycetales', 'Saccharomycetaceae',
             'Saccharomyces', 'Saccharomyces b', 's', '90', 'x']]
    with open(panel, 'w') as o:
        o.write('\t'.join(cols) + '\n')
        for r in rows:
            o.write('\t'.join(r) + '\n')
    hits = tmp_path / 'hits'
    hits.mkdir()
    (hits / 'Q1__R1.tsv').write_text('a\tx\t75\t100\t100\t100\t1e-50\t300\nb\tx\t85\t100\t100\t100\t1e-50\t300\n')
    (hits / 'Q1__S1.tsv').write_text('a\tx\t40\t100\t100\t100\t1e-20\t100\n')
    (hits / 'Q1__Q1.tsv').write_text('a\tx\t100\t100\t100\t100\t0\t500\n')
    nprot = tmp_path / 'nprot.tsv'
    nprot.write_text('Q1\t4\n')
    out = tmp_path / 'aai.tsv'
    subprocess.run([sys.executable, str(HERE / 'aai_summary.py'), '--panel', str(panel), '--hits', str(hits),
                    '--nprot', str(nprot), '--out', str(out)], check=True)
    got = list(csv.DictReader(open(out), delimiter='\t'))
    assert [r['ref'] for r in got] == ['R1', 'S1']          # self excluded, sorted by AAI
    assert got[0]['rank'] == '1' and got[0]['aai'] == '80.0' and got[0]['frac_hit'] == '0.5'
    assert got[1]['ref_class'] == 'Saccharomycetes'
