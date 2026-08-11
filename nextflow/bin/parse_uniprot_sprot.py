#!/usr/bin/env python3
"""Parse the UniProtKB/Swiss-Prot flatfile (uniprot_sprot.dat.gz) into a
per-accession annotation CSV for the SwissProt BLAST/diamond results table.

Output columns (one row per primary accession):
    accession, entry_name, protein_name, gene_name, gene_synonyms,
    organism, go_terms, ec_numbers, interpro, pfam

Multi-valued fields are joined with '; '. GO / InterPro / Pfam entries carry
both the accession and the term name, e.g. "GO:0004674 (Protein
serine/threonine kinase activity)". This is dependency-free (pure stdlib)
so it runs under any python3, matching merge_swissprot.py.
"""

import gzip
import re
import csv
import sys
import argparse

_EC_RE = re.compile(r'EC=([0-9]+(?:\.[0-9]+)+)')

HEADER = [
    'accession', 'entry_name', 'protein_name', 'gene_name', 'gene_synonyms',
    'organism', 'go_terms', 'ec_numbers', 'interpro', 'pfam',
]


def uniq(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def parse_entry(lines):
    """lines: raw lines of one //-terminated entry. Returns a row dict."""
    accs = []
    entry_name = ''
    recname = []
    gene_name = ''
    gene_synonyms = []
    organism = []
    go = []
    interpro = []
    pfam = []
    ec = []

    for raw in lines:
        ln = raw.rstrip('\n')
        if len(ln) < 2:
            continue
        code = ln[:2]
        body = ln[2:].strip()

        if code == 'ID':
            entry_name = body.split()[0] if body.split() else ''
        elif code == 'AC':
            accs.extend(t.strip().rstrip(';') for t in body.split(';') if t.strip())
        elif code == 'DE':
            m = re.search(r'RecName: Full=([^;{]+)', body)
            if m:
                recname.append(m.group(1).strip())
        elif code == 'GN':
            m = re.search(r'Name=([^;{]+)', body)
            if m and not gene_name:
                gene_name = m.group(1).strip()
            for gm in re.finditer(r'Synonyms=([^;{]+)', body):
                gene_synonyms.extend(s.strip() for s in gm.group(1).split(',') if s.strip())
        elif code == 'OS':
            organism.append(body)
        elif code == 'DR':
            fields = [f.strip() for f in body.split(';')]
            if len(fields) >= 3:
                db, dbid = fields[0], fields[1]
                name = fields[2].rstrip('.')
                if db == 'GO' and dbid.startswith('GO:'):
                    go.append(f"{dbid} ({name})" if name else dbid)
                elif db == 'InterPro' and dbid.startswith('IPR'):
                    interpro.append(f"{dbid} ({name})" if name else dbid)
                elif db == 'Pfam' and dbid.startswith('PF'):
                    pfam.append(f"{dbid} ({name})" if name else dbid)
        # EC numbers also appear on DE lines (EC=2.7.11.1 {ECO:...}).
        if code in ('DE', 'CC'):
            ec.extend(m.group(1) for m in _EC_RE.finditer(body))

    return {
        'accession': accs[0] if accs else '',
        'entry_name': entry_name,
        'protein_name': ' / '.join(uniq(recname)),
        'gene_name': gene_name,
        'gene_synonyms': ', '.join(uniq(gene_synonyms)),
        'organism': ' '.join(organism).rstrip('.') if organism else '',
        'go_terms': '; '.join(uniq(go)),
        'ec_numbers': '; '.join(uniq(ec)),
        'interpro': '; '.join(uniq(interpro)),
        'pfam': '; '.join(uniq(pfam)),
    }


def main():
    parser = argparse.ArgumentParser(description="Parse Swiss-Prot flatfile to CSV")
    parser.add_argument("dat", help="uniprot_sprot.dat.gz (or plain .dat)")
    parser.add_argument("-o", "--outfile", default="swissprot_annot.csv")
    args = parser.parse_args()

    opener = gzip.open if args.dat.endswith('.gz') else open
    n = 0
    with opener(args.dat, 'rt') as fh, open(args.outfile, 'w', newline='') as of:
        w = csv.DictWriter(of, fieldnames=HEADER)
        w.writeheader()
        entry = []
        for ln in fh:
            if ln.strip() == '//':
                row = parse_entry(entry)
                w.writerow(row)
                n += 1
                entry = []
            elif not ln.strip():
                continue
            else:
                entry.append(ln)
        if entry:  # trailing unterminated entry (should not happen)
            w.writerow(parse_entry(entry))
            n += 1

    print(f"Parsed {n} entries -> {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
