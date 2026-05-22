#!/usr/bin/env python3

import json
import csv
import argparse
import re

def load_json(file_path):
    with open(file_path, 'r') as fh:
        data = json.load(fh)
    return data

def load_samples(file_path):
    samples = []
    with open(file_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            samples.append(row)
    return samples

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate funguilds data for bigquery",
                                    epilog='Example: build_funguild_bigquery.py')
    parser.add_argument("--funguild", default="lib/funguild.pp.json", help="funguild json")
    parser.add_argument("--samples", default="samples.csv", help="samples.csv file for fungi5k")
    parser.add_argument('-v','--debug', help='Debugging output', action='store_true')

    parser.add_argument("-o","--outfile", default='bigquery/species_funguild.csv', 
                        help="Output file")
    
    args = parser.parse_args()
    data = load_json(args.funguild)
    samples = load_samples(args.samples)    
    
    with open(args.outfile, "w",newline="") as outfh:
        csvout = csv.writer(outfh)
        csvout.writerow(["species_prefix","species","growthForm","guild","trophicMode","confidenceRanking"])
        guilds = {}
        for d in data:
            taxon = d['taxon']
            if taxon not in guilds:
                guilds[taxon] = d
            else:
                print(f"Duplicate taxon: {taxon}")            
            
            # print(d['taxon'])
        
        for sample in samples:
            genus = sample['GENUS']
            sp  = sample['SPECIES']
            if not genus:
                genus = sp.split(" ")[0]
            wrote = False
            for itype in [sp, genus]:
                if itype in guilds:
                    drow = guilds[itype]
                    drow['trophicMode'] = re.sub('^\s+','',drow['trophicMode'])
                    csvout.writerow([sample['LOCUSTAG'],itype,drow['growthForm'],drow['guild'],drow['trophicMode'],drow['confidenceRanking']])
                    wrote = True
                    break
            if not wrote:
                print(f"Missing guild {genus} {sp}")
            
            
