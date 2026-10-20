import os, sys, csv, glob
sys.path.insert(0, os.environ['NEW'] + '/aldy_bench')
sys.path.insert(0, os.environ['SRC'] + '/benchmark')
from aldy_norm import parse_aldy_file, clean_aldy_diplotype
from run_single_iteration import normalize_genotype, compare_genotype
import benchmark_config as cfg

NEW = os.environ['NEW']
WD = f'{NEW}/aldy_bench'

gold = {}
for arm in ('cohort1', 'cohort2'):
    p = f"{os.environ['SRC']}/data/CYP2D6_{arm}_gold_standard_v2.csv"
    gold[arm] = {r['Sample_ID'].strip(): r['Genotype'].strip()
                 for r in csv.DictReader(open(p))}

rows = []
for arm in ('cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen'):
    base = arm.replace('_dragen', '')
    aligner = 'DRAGEN' if arm.endswith('_dragen') else 'BWA-MEM'
    for d in sorted(glob.glob(f'{WD}/{arm}/ref_*')):
        ref = os.path.basename(d).replace('ref_', '')
        for f in sorted(glob.glob(f'{d}/*.cyp2d6.aldy')):
            s = os.path.basename(f).replace('.cyp2d6.aldy', '')
            major, nsol = parse_aldy_file(f)
            cleaned = clean_aldy_diplotype(major)
            g = gold[base].get(s)
            # 多解一律視為 discordant（與 CyriPanel 同標準）
            if nsol > 1:
                comp = 'discordant'
            elif cleaned is None or g is None:
                comp = 'no_call' if g is not None else 'no_gold'
            else:
                comp = compare_genotype(cleaned, g)
            rows.append(dict(arm=arm, cohort=base, aligner=aligner, reference=ref,
                             sample=s, n_solutions=nsol, raw_major=major or '',
                             cleaned=cleaned or '',
                             normalised=normalize_genotype(cleaned) if cleaned else '',
                             gold=g or '', gold_norm=normalize_genotype(g) if g else '',
                             comparison=comp))

out = f'{WD}/aldy_results_long.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f'wrote {len(rows)} rows -> {out}\n')

from collections import Counter, defaultdict
print(f'{"arm":18}{"runs":>7}{"concordant":>12}{"conc%":>9}{"multi-sol":>11}{"ms%":>7}')
print('-' * 66)
for arm in ('cohort1', 'cohort1_dragen', 'cohort2', 'cohort2_dragen'):
    r = [x for x in rows if x['arm'] == arm]
    c = sum(1 for x in r if x['comparison'] == 'concordant')
    m = sum(1 for x in r if x['n_solutions'] > 1)
    print(f'{arm:18}{len(r):>7}{c:>12}{c/len(r)*100:>8.1f}%{m:>11}{m/len(r)*100:>6.1f}%')
