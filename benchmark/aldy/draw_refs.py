"""Reproduce CyriPanel's reference-panel draw, and derive Aldy's single reference."""
import os, sys, glob, csv, random, json
sys.path.insert(0, os.environ['SRC'] + '/benchmark')
import benchmark_config as cfg


def build_pool(cohort):
    c = cfg.COHORTS[cohort]
    with open(c['diploid_csv']) as f:
        diploid_ids = sorted(r['Sample_ID'].strip() for r in csv.DictReader(f))
    suffix = c['bam_suffix']
    bams = glob.glob(os.path.join(c['bam_dir'], '*', '*' + suffix)) \
        or glob.glob(os.path.join(c['bam_dir'], '*' + suffix))
    all_ids = {os.path.basename(b).replace(suffix, ''): b for b in bams}
    valid = [s for s in diploid_ids if s in all_ids]
    return valid, all_ids


print('=== 驗證能否重現 CyriPanel 的抽樣 ===')
ROOT = os.environ['ROOT']
ok = bad = 0
for cohort in ('cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen'):
    valid, all_ids = build_pool(cohort)
    for p in (5, 10, 20):
        for it in (0, 7, 15, 29):
            fp = f'{ROOT}/{cohort}/panel{p}/iter_{it:03d}/reference_samples.json'
            if not os.path.exists(fp):
                continue
            actual = json.load(open(fp))
            actual = sorted(actual.get('references', actual) if isinstance(actual, dict) else actual)
            rng = random.Random(42 + it)
            mine = sorted(rng.sample(valid, p))
            if mine == actual:
                ok += 1
            else:
                bad += 1
                print(f'  MISMATCH {cohort} n={p} iter={it}')
                print(f'     actual: {actual}')
                print(f'     mine  : {mine}')
print(f'  matched {ok}, mismatched {bad}')

print()
print('=== Aldy 的單一 reference（seed=42+i, k=1）===')
for cohort in ('cohort1', 'cohort2'):
    valid, all_ids = build_pool(cohort)
    picks = [random.Random(42 + i).sample(valid, 1)[0] for i in range(30)]
    from collections import Counter
    c = Counter(picks)
    print(f'  {cohort}: pool={len(valid)}  30 iters 用到 {len(c)} 個不同 reference')
    print(f'     前 10 iter: {picks[:10]}')
    print(f'     最常抽到  : {c.most_common(3)}')
