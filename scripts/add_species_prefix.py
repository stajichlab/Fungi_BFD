#!/usr/bin/env python

import sys
import argparse
import csv
import re

parser = argparse.ArgumentParser(description = 'Add species prefix to gene names')

parser.add_argument('-s','--samples', type = str, default = 'samples.csv', help = 'species prefix')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

samples = {}
with open(args.samples, 'r') as fh:
    sampleinfo = csv.DictReader(fh, delimiter=",")
    for row in sampleinfo:
        samples[row['LOCUSTAG']] = row['CLASS'] + '_' + re.sub(' ','_',row['SPECIES'])
        
idmatch = re.compile(r'>([^_]+)_(\d+)')
for line in args.infile:
    line = line.rstrip()
    if line.startswith('>'):
        m = idmatch.match(line)
        if m:
            locustag = m.group(1)
            geneid = m.group(2)
            if locustag in samples:
                species = samples[locustag]
                print(f'>{species}__{locustag}_{geneid}', file=args.outfile)
            else:
                print(f'{line}', file=args.outfile)
        else:
            print(f'Error: no match for ID {line}',file=sys.stderr)
            print(line, file=args.outfile)
            continue
    else:
        print(line, file=args.outfile)
