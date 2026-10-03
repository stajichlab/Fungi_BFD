#!/usr/bin/env python3
"""Build per-genome ab-initio parameter stores from earlier BUSCO-trained predictions, for
genomes with no RNA-seq (Fungi_BFD_runs DECISIONS D136). FUNANNOTATE_PREDICTION passes
<store>/<out>/parameters.json to `funannotate predict -p` (utils.nf ownAbinitioParamsFor),
so predict skips BUSCO and Augustus/SNAP training.

A genome qualifies when all hold:
  - it is in the run's samples CSV and has no GBK in <run>/genome_annotation yet;
  - its species has no RNA-seq reads (<run>/rnaseq_reads/<species>_norm_R1.fastq.gz absent
    or empty) -- RNA-seq genomes always retrain;
  - <source>/<out>/predict_results/<out>.parameters.json says Augustus AND SNAP were trained
    from BUSCO (funannotate writes "BUCSCO <lineage>");
  - the source copies exist: predict_misc/ab_initio_parameters/augustus/species/<out>/ and
    <out>.snap.hmm;
  - logfiles/funannotate-predict.log reports >= --min-busco valid BUSCO training models.

Store layout, same as gene_prediction_shared_abinitio (relative paths, relocatable):
  <out>/parameters.json, <out>/<out lower>/ (Augustus species), <out>/<out lower>.snap.hmm,
  <out>/glimmerhmm_stub/ (empty; keeps predict from running BUSCO for GlimmerHMM),
  <out>/provenance.json. GeneMark is left empty: GENEMARK_RUN supplies the GTF.

Dry run by default (prints per-genome status); --apply writes the stores and <store>/build_report.tsv.
"""
import argparse, csv, glob, json, os, re, shutil, datetime

ap = argparse.ArgumentParser()
ap.add_argument('--run-dir', required=True)
ap.add_argument('--samples', help='default <run-dir>/samples_wave1.csv')
ap.add_argument('--source', required=True, help='annotation tree with the earlier predictions')
ap.add_argument('--store', help='default <run-dir>/abinitio_own')
ap.add_argument('--min-busco', type=int, default=200)
ap.add_argument('--apply', action='store_true')
a = ap.parse_args()
run = os.path.abspath(a.run_dir)
store = os.path.abspath(a.store or os.path.join(run, 'abinitio_own'))
samples = a.samples or os.path.join(run, 'samples_wave1.csv')


def clean_strain(s):
    s = re.sub(r"['\"]", '', (s or '').strip()).split(';')[0].strip().replace(':', ' ')
    s = re.sub(r'^\s*\*+', '', s); s = re.sub(r'\*+\s*$', '', s); s = re.sub(r'\s*\*+\s*', '-', s)
    return s.strip()


def tag(sp, st):  # nextflow/modules/common/utils.nf makeSampleTag()
    sp = re.sub(r"['\"]", '', (sp or '').strip())
    return re.sub(r'[\s/#\[\]?{}]+', '_', '_'.join(x for x in (sp, clean_strain(st)) if x))


report, n = [], 0
for r in csv.DictReader(open(samples)):
    out = tag(r['SPECIES'], r['STRAIN']); low = out.lower()
    if glob.glob(f'{run}/genome_annotation/{out}/predict_results/{out}.gbk*'):
        continue
    n += 1
    spt = re.sub(r'\s+', '_', r['SPECIES'].strip())
    reads = f'{run}/rnaseq_reads/{spt}_norm_R1.fastq.gz'
    if os.path.exists(reads) and os.path.getsize(reads) > 0:
        continue                                    # RNA-seq: always retrain
    src = f'{a.source}/{out}'
    pj = f'{src}/predict_results/{low}.parameters.json'
    aug = f'{src}/predict_misc/ab_initio_parameters/augustus/species/{low}'
    snap = f'{src}/predict_misc/ab_initio_parameters/{low}.snap.hmm'
    status, nbusco, a0 = '', '', {}
    if not os.path.exists(pj):
        status = 'no_parameters_json'
    else:
        d = json.load(open(pj))
        a0, s0 = (d.get('augustus') or [{}])[0], (d.get('snap') or [{}])[0]
        if not (a0.get('source', '').startswith('BUCSCO') and s0.get('source', '').startswith('BUCSCO')):
            status = f"not_busco_trained(augustus={a0.get('source', 'none')[:20]})"
        elif not (os.path.isdir(aug) and os.listdir(aug) and os.path.exists(snap) and os.path.getsize(snap) > 0):
            status = 'source_files_missing'
        else:
            m = None
            log = f'{src}/logfiles/funannotate-predict.log'
            if os.path.exists(log):
                for line in open(log, errors='replace'):
                    x = re.search(r'([\d,]+) valid BUSCO predictions', line)
                    if x:
                        m = int(x.group(1).replace(',', ''))
            nbusco = '' if m is None else str(m)
            if m is None:
                status = 'busco_count_unknown'
            elif m < a.min_busco:
                status = f'busco_models_below_{a.min_busco}'
            else:
                status = 'ok'
    report.append([out, r['ASMID'], status, nbusco, a0.get('source', ''), a0.get('version', ''), a0.get('date', '')])
    if status != 'ok' or not a.apply:
        continue
    dest = os.path.join(store, out)
    tmp = os.path.join(store, f'.{out}.staging')
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    shutil.copytree(aug, os.path.join(tmp, low))
    shutil.copy2(snap, os.path.join(tmp, f'{low}.snap.hmm'))
    os.makedirs(os.path.join(tmp, 'glimmerhmm_stub'))
    src_tag = a0['source']
    json.dump({'augustus': [{'source': f'own-reuse {src_tag}', 'path': low}],
               'genemark': [{}],
               'snap': [{'source': f'own-reuse {src_tag}', 'path': f'{low}.snap.hmm'}],
               'codingquarry': [{}],
               'glimmerhmm': [{'source': 'suppressed-empty-stub', 'path': 'glimmerhmm_stub'}],
               'table': int(r.get('TRANSL_TABLE') or 1)},
              open(os.path.join(tmp, 'parameters.json'), 'w'), indent=2)
    json.dump({'out': out, 'asmid': r['ASMID'], 'source_predict_dir': src,
               'source_parameters_json': pj, 'source_version': a0.get('version', ''),
               'source_date': a0.get('date', ''), 'valid_busco_models': int(nbusco),
               'generated_at': datetime.datetime.now().isoformat(timespec='seconds')},
              open(os.path.join(tmp, 'provenance.json'), 'w'), indent=2)
    shutil.rmtree(dest, ignore_errors=True)
    os.rename(tmp, dest)

rep_path = os.path.join(store, "build_report.tsv")
counts = {}
for x in report:
    k = x[2].split('(')[0]
    counts[k] = counts.get(k, 0) + 1
print(f'genomes without a GBK: {n}; without RNA-seq: {len(report)}')
for k, v in sorted(counts.items()):
    print(f'  {k}: {v}')
if a.apply:
    with open(rep_path, 'w') as fh:
        fh.write('out\tasmid\tstatus\tvalid_busco_models\tsource\tversion\tdate\n')
        for x in report:
            fh.write('\t'.join(x) + '\n')
    print(f"wrote {counts.get('ok', 0)} stores under {store}; report {rep_path}")
else:
    for x in report:
        print('\t'.join(x))
    print('dry run; rerun with --apply to write the stores')
