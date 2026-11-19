import csv,sys,collections
genomes=sys.argv[1:]
MINLEN=int(1e5)
def fa_lens(p):
    L={};n=None
    for l in open(p):
        if l[0]=='>': n=l[1:].split()[0]; L[n]=0
        else: L[n]+=len(l.strip())
    return L
for g in genomes:
    L=fa_lens(f'{g}/genome.fa')
    big={k:v for k,v in L.items() if v>=MINLEN}
    ends=[(s,e) for s in big for e in ('5prime','3prime')]
    fnd=set(); 
    for r in csv.DictReader(open(f'{g}/find.tsv'),delimiter='\t'):
        if r['terminal']=='True': fnd.add((r['scaffold'],r['end']))
    dn=set()
    for r in csv.DictReader(open(f'{g}/denovo_ends.tsv'),delimiter='\t'):
        dn.add((r['scaffold'],r['end']))
    # tidk search windows: end if first/last window has >=TH repeats (either strand)
    TH=5
    win=collections.defaultdict(list)
    for r in csv.DictReader(open(f'{g}/tidk_TTAGGG_telomeric_repeat_windows.csv')):
        win[r['id']].append((int(r['window']),int(r['forward_repeat_number'])+int(r['reverse_repeat_number'])))
    ts=set()
    for s,w in win.items():
        if s in big:
            if w[0][1]>=TH: ts.add((s,'5prime'))
            if w[-1][1]>=TH or (len(w)>1 and w[-1][1]+w[-2][1]>=TH and False): ts.add((s,'3prime'))
    # tidk explore: any hit within 2kb of an end
    te=set()
    for r in csv.DictReader(open(f'{g}/tidk_explore_telomeric_locations.tsv'),delimiter='\t'):
        s=r['id']
        if s in big:
            if int(r['start_pos'])<2000: te.add((s,'5prime'))
            if int(r['end_pos'])>L[s]-2000: te.add((s,'3prime'))
    print(f'== {g}: {len(big)} scaffolds >=100kb, {len(ends)} ends (total scaffolds {len(L)})')
    for name,S in (('find_telomeres',fnd),('denovo',dn),('tidk search TTAGGG w1000>=5',ts),('tidk explore <2kb',te)):
        print(f'  {name:32s} {len({e for e in ends if e in S}):3d}/{len(ends)}')
    allf=fnd|dn|ts
    miss=[e for e in ends if e not in fnd]
    print('  ends not in find but in denovo/tidk:',[(e,'dn' if e in dn else '', 'ts' if e in ts else '') for e in miss if e in dn or e in ts][:10])
