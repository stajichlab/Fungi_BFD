#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Updated by Jason Stajich for purpose of downloading db only.

Based on Guilds.py from Zewei Song.

optional arguments:
  -h, --help       Show this help message and exit
  --out lib/funguild         Path to funguild saves
  --db              Database to use ('fungi' or 'nematode') [default:fungi]
  
This is an example command to run this script:
python download_funguild.py [--out lib/funguild] [--db fungi]

'''
from __future__ import print_function
from __future__ import division
#Import modules#################
import argparse
import os
import timeit
import sys
#import urllib
from operator import itemgetter
import csv

start = timeit.default_timer()
################################

#Command line parameters#####################################################################
parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", help="Stem for output file.", default='lib/funguild')
parser.add_argument("-d", "--db", choices=['fungi','nematode'], default='fungi', 
                    help="Assign a specified database to the script")
args = parser.parse_args()

outbase = args.output
outfolder = os.path.dirname(outbase)
# add protection in future
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

database_name = args.db
if database_name == 'fungi':
    url = 'http://www.stbates.org/funguild_db_2.php'
elif database_name == 'nematode':
    url = 'http://www.stbates.org/nemaguild_db.php'

import requests
import json

print('Connecting with FUNGuild database ...')
db_url = requests.get(url)
#db_url = db_url.content.decode('utf-8').split('\n')[6].strip('[').strip(']</body>').replace('} , {', '} \n {').split('\n')
db_url = db_url.content.decode('utf-8')
db_url = db_url.split('\n')[6].strip('</body>')
data = json.loads(db_url)

with open(f'{outbase}.json', 'w') as f:
    json.dump(data, f)

with open(f'{outbase}.pp.json', 'w') as f:
    json.dump(data, f, indent=4) # 'indent=4' makes the JSON output human-readable with 4-space indentation
