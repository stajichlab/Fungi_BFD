#!/usr/bin/env python3

# this script is from 1KFG/common_annotate but demonstrates how the LOCUS tag prefix was generated based on the assembly id

import os
import csv
import xml.etree.ElementTree as ET
import hashlib
from Bio import Entrez
Entrez.email = 'jason.stajich@ucr.edu'
from joblib import Memory

cachedir = os.path.join(os.getcwd(), 'cache')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)
memory = Memory(cachedir, verbose=0)


#@memory.cache
def get_locus_from_asm(assembly_id):
    print(f"asmid = '{assembly_id}'")
    mdsum = hashlib.md5(assembly_id.encode('utf-8')).hexdigest().upper()
    return 'F' + mdsum[-7:]
    
@memory.cache
def get_bioproject_prefix(BIOPROJECTID):
    LOCUSTAG = ""
    for project in BIOPROJECTID.split(';'):
        # remove umbrella projects
        print(f'Checking {project}')
        if project == "PRJNA533106" or project == "PRJEB40665" or project == "PRJEB43510":
            continue
        try:
            bioproject_handle = Entrez.efetch(db="bioproject",id = project)
            projtree = ET.parse(bioproject_handle)
            projroot = projtree.getroot()
        
            lt = projroot.iter('LocusTagPrefix')
            for locus in lt:
                LOCUSTAG = locus.text
        except Exception as e:
            print(f'Error {e}')
        if len(LOCUSTAG) > 0:
            break
    return LOCUSTAG

# Families where all members use the Alternative Yeast Nuclear Code (CUG→Ser, table 12)
TRANSL_TABLE_12_FAMILIES = {'Debaryomycetaceae'}
# Individual genera using table 12 not covered by the family-level rule above
TRANSL_TABLE_12_GENERA = {'Clavispora'}
# Genera using Pachysolen tannophilus Nuclear Code (CUG→Ala, table 26)
TRANSL_TABLE_26_GENERA = {'Pachysolen'}

def get_transl_table(family, genus):
    if family in TRANSL_TABLE_12_FAMILIES or genus in TRANSL_TABLE_12_GENERA:
        return 12
    if genus in TRANSL_TABLE_26_GENERA:
        return 26
    return 1

accessions = 'ncbi_accessions.csv'
accession_taxonomy = 'ncbi_accessions_taxonomy.csv'

outsamples = 'samples.csv'

accession_dict = {}
fields = ['ASMID','SPECIES_IN', 'STRAIN', 'BIOPROJECT', 'NCBI_TAXONID', 'BUSCO_LINEAGE']
with open(accessions, 'r', newline='') as f:
    reader = csv.reader(f)
    header = next(reader)
    hdr2dict = {header[i]: i for i in range(len(header))}
    for row in reader:
        acc = row[hdr2dict['ACCESSION']]
        species = row[hdr2dict['SPECIES']]
        strain = row[hdr2dict['STRAIN']]
        asm_name = row[hdr2dict['ASM_NAME']]
        asm_base = f'{acc}_{asm_name}'
        accession_dict[asm_base] = [species, strain, 
                                    row[hdr2dict['BIOPROJECT']], 
                                    row[hdr2dict['NCBI_TAXID']] ]
        
with open(accession_taxonomy, 'r', newline='') as f:
    reader = csv.reader(f)
    header = next(reader)
    fields.extend(header[4:])
    fields.append('TRANSL_TABLE')
    tax_hdr = {header[i]: i for i in range(len(header))}
    for row in reader:
        asm_base = row[0]

        # add the BUSCO lineage to the dictionary
        buscodb = 'fungi'
        if row[4] == "Ascomycota" or row[4] == "Basidiomycota":
            buscodb = 'dikarya'
        elif row[4] == "Mucoromycota":
            buscodb = 'mucoromycota'
        elif row[4] == 'Microsporidia':
            buscodb = 'microsporidia'
        if asm_base not in accession_dict:
            print(f'{asm_base} not found in accession_dict, skipping')
            continue
        family = row[tax_hdr['FAMILY']] if 'FAMILY' in tax_hdr else ''
        genus = row[tax_hdr['GENUS']] if 'GENUS' in tax_hdr else ''
        transl_table = get_transl_table(family, genus)
        accession_dict[asm_base].append(buscodb)
        accession_dict[asm_base].extend(row[4:])
        accession_dict[asm_base].append(transl_table)
            
            
print(f'{len(accession_dict)} assemblies processed')

# add LOCUS tag field
fields.append('LOCUSTAG')

# add LOCUSTAG to the dictionary
# seen maps base hash -> collision count; first occurrence keeps no suffix,
# subsequent collisions get A, B, C ... to guarantee uniqueness.
seen = {}  # base_hash -> count of times seen so far
for asm in accession_dict:
    locus = get_locus_from_asm(asm)
    if locus in seen:
        collision_index = seen[locus]  # 1 for second occurrence → suffix 'A'
        suffixed = locus + chr(ord('A') + collision_index - 1)
        print(f'Collision: {locus} reassigned to {suffixed} for {asm}')
        seen[locus] += 1
        locus = suffixed
    else:
        seen[locus] = 1
    accession_dict[asm].append(locus)
    
with open(outsamples, 'wt',newline="") as f:
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow(fields)
    for asm in accession_dict:
        writer.writerow([asm] + accession_dict[asm])

