import os, sys, glob, csv
sys.path.insert(0, os.environ['SRC'] + '/benchmark')
import benchmark_config as cfg

NEW = os.environ['NEW']
OUT = f'{NEW}/aldy_bench'
os.makedirs(OUT, exist_ok=True)

rows = []
for arm in ('cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen'):
    c = cfg.COHORTS[arm]
    dip = sorted(r['Sample_ID'].strip() for r in csv.DictReader(open(c['diploid_csv'])))
    suf = c['bam_suffix']
    bams = glob.glob(os.path.join(c['bam_dir'], '*', '*' + suf)) \
        or glob.glob(os.path.join(c['bam_dir'], '*' + suf))
    ids = {os.path.basename(b).replace(suf, ''): b for b in bams}
    valid = [s for s in dip if s in ids]
    for ref in valid:
        rows.append((arm, ref, ids[ref], str(len(ids) - 1)))
    print(f'{arm:16} refs={len(valid):3}  tests/ref={len(ids)-1:3}  runs={len(valid)*(len(ids)-1):5}')

with open(f'{OUT}/manifest_refs.txt', 'w') as f:
    for r in rows:
        f.write('\t'.join(r) + '\n')

print(f'\ntasks = {len(rows)}')
print(f'total runs = {sum(int(r[3]) for r in rows)}')
print(f'manifest: {OUT}/manifest_refs.txt')
