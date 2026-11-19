#!/usr/bin/env python3
"""Score telomere callers against a proxy truth. See TIDK_COMPARISON.md.
usage: score_ground_truth.py LIST.tsv RUNDIR [--copies 8 --win 300 --minlen 100000]
Truth: scaffold (>= minlen) end whose terminal --win bp holds >= --copies copies of the
genome's unit (or its reverse complement). Unit = denovo candidate with most ends."""
import csv, re, sys, argparse, collections
ap=argparse.ArgumentParser(); ap.add_argument('lst'); ap.add_argument('run')
ap.add_argument('--copies',type=int,default=8); ap.add_argument('--win',type=int,default=300)
ap.add_argument('--minlen',type=int,default=100000)
ap.add_argument('--find',nargs='+',default=['find_tol10:find10.tsv','find_tol30:find30.tsv'],help='name:file pairs in each genome dir')
a=ap.parse_args()
RC=str.maketrans('ACGT','TGCA')
def rc(s): return s.translate(RC)[::-1]
def read(p):
    d={};n=None
    for l in open(p):
        if l[0]=='>': n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l.strip().upper())
    return {k:''.join(v) for k,v in d.items()}
tot=collections.Counter(); per=[]
for line in open(a.lst):
    asm=line.split('\t')[0]; o=f'{a.run}/{asm}'
    try: seqs=read(f'{o}/genome.fa')
    except FileNotFoundError: continue
    units=[(int(r['n_ends']),r['unit']) for r in csv.DictReader(open(f'{o}/denovo.tsv'),delimiter='\t') if r['unit']!='.']
    if not units: per.append((asm,'no-unit')); continue
    unit=max(units)[1]
    # rotation-aware: count occurrences of any rotation by doubling
    pats=set()
    for s in (unit,rc(unit)):
        for i in range(len(s)): pats.add(s[i:]+s[:i])
    rx=re.compile('|'.join(sorted(pats)))
    truth=set()
    for s,q in seqs.items():
        if len(q)<a.minlen: continue
        if len(rx.findall(q[:a.win]))>=a.copies: truth.add((s,'5prime'))
        if len(rx.findall(q[-a.win:]))>=a.copies: truth.add((s,'3prime'))
    ends={(s,e) for s,q in seqs.items() if len(q)>=a.minlen for e in ('5prime','3prime')}
    calls={}
    for k,f in (x.split(':') for x in a.find):
        calls[k]={(r['scaffold'],r['end']) for r in csv.DictReader(open(f'{o}/{f}'),delimiter='\t')
                  if r.get('terminal')=='True' and (r['scaffold'],r['end']) in ends} \
                 if True else set()
    calls['denovo']={(r['scaffold'],r['end']) for r in csv.DictReader(open(f'{o}/denovo_ends.tsv'),delimiter='\t') if (r['scaffold'],r['end']) in ends}
    row=[asm,unit,len(ends),len(truth)]
    for k,c in calls.items():
        tp=len(c&truth); fp=len(c-truth); row+=[tp,fp]; tot[k+'_tp']+=tp; tot[k+'_fp']+=fp
    tot['truth']+=len(truth); tot['ends']+=len(ends); per.append(tuple(row))
w=csv.writer(sys.stdout,delimiter='\t'); w.writerow(['asm','unit','ends','truth']+[f'{k}_{x}' for k in calls for x in ('tp','fp')])
for r in per: w.writerow(r)
print('# TOTAL ends',tot['ends'],'truth',tot['truth'],file=sys.stderr)
for k in [x.split(':')[0] for x in a.find]+['denovo']:
    tp,fp=tot[k+'_tp'],tot[k+'_fp']
    print(f'# {k}: recall {tp}/{tot["truth"]}={tp/max(1,tot["truth"]):.3f}  false-pos(not in truth) {fp}',file=sys.stderr)
print('# genomes without unit:',sum(1 for p in per if len(p)==2),file=sys.stderr)
