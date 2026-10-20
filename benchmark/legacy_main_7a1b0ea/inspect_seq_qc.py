#!/usr/bin/env python3
"""
inspect_seq_qc.py
=================
Read CyriPanel_seq_qc.tsv and answer the three questions it was run for:

  1. Table 1 values: per-arm on-target %, duplicate rate, insert size, CYP2D6 depth.
  2. Does the PCR-free vs PCR-based library chemistry show up in the duplicate rate,
     and is cohort 1 really the higher one?
  3. Do the samples CyriPanel called wrongly have worse sequencing metrics than the
     rest of their cohort, and does cohort 1's sequencing batch matter?

No scipy. Percentile ranks rather than a test: with 5 discordant samples in cohort 1
a formal test is underpowered, and the rank tells you more.
"""
import os, sys
import pandas as pd

ROOT = '/work/u7715055/staging/biology/u7715055'
QC   = os.path.join(ROOT, 'CyriPanel_seq_qc.tsv')
CALLS= os.path.join(ROOT, 'cyripanel_benchmark_results/recompute_v2/table3/sample_level_calls.csv')

# cohort 1 sequencing runs (from the core facility notes + BAM listing)
RUNS = {
 'PG_20230217': set('HG00733 NA19240 NA24385 NA24631'.split()),
 'PG_20230428': set('''HG00735 HG00741 HG01071 HG01106 HG01175 HG01258 HG01358 HG01361
                       HG01891 HG01928 HG01952 HG01978 HG02148 HG02257 HG02572 HG02622
                       HG02630 HG02717 HG02886 HG03453 HG03516 HG03540 HG03579'''.split()),
 'PG_20230512': set('''HG00438 HG00621 HG00673 HG01109 HG01243 HG02055 HG02080 HG02109
                       HG02145 HG02723 HG02818 HG03098 HG03486 HG03492 NA18906 NA20129
                       NA21309'''.split()),
}
NUM = ['total_mapped','on_target','pct_on_target','dup_on_target','pct_dup_on_target',
       'insert_median','insert_mean','insert_sd','depth_locus','depth_CYP2D6','cov_pct_CYP2D6']

def med_iqr(s):
    s = s.dropna()
    if s.empty: return 'n/a'
    return '%8.1f  [%.1f – %.1f]  n=%d' % (s.median(), s.quantile(.25), s.quantile(.75), len(s))

qc = pd.read_csv(QC, sep='\t')
for c in NUM:
    if c in qc: qc[c] = pd.to_numeric(qc[c], errors='coerce')
qc['cohort']  = qc['arm'].str.replace('_dragen','', regex=False)
qc['aligner'] = qc['arm'].apply(lambda a: 'DRAGEN' if a.endswith('_dragen') else 'BWA-MEM')
qc['run'] = qc.apply(lambda r: next((k for k,v in RUNS.items() if r['sample'] in v), '—')
                     if r['cohort']=='cohort1' else 'PG_20250711', axis=1)

print('=' * 78)
print('  Rows: %d   (expected 232 = 44+72+44+72)' % len(qc))
bad = qc[qc[['pct_on_target','pct_dup_on_target','depth_CYP2D6']].isna().any(axis=1)]
if len(bad):
    print('  ROWS WITH MISSING VALUES (%d):' % len(bad))
    print(bad[['arm','sample','pct_on_target','pct_dup_on_target','insert_median','depth_CYP2D6']]
          .to_string(index=False))
for arm, g in qc.groupby('arm'):
    print('  %-16s %d samples' % (arm, len(g)))
if (qc['run']=='—').any():
    print('  UNMAPPED to a run:', sorted(qc.loc[qc['run']=='—','sample'].unique()))

print()
print('=' * 78)
print('  1. Table 1 values, median [IQR]')
print('=' * 78)
for arm in ['cohort1','cohort2','cohort1_dragen','cohort2_dragen']:
    g = qc[qc['arm']==arm]
    if g.empty: continue
    print('\n  %s  (n=%d)' % (arm, len(g)))
    print('     on-target %%        %s' % med_iqr(g['pct_on_target']))
    print('     duplicate %%        %s' % med_iqr(g['pct_dup_on_target']))
    print('     insert size (bp)   %s' % med_iqr(g['insert_median']))
    print('     depth, CYP2D locus %s' % med_iqr(g['depth_locus']))
    print('     depth, CYP2D6      %s' % med_iqr(g['depth_CYP2D6']))

print()
print('=' * 78)
print('  2. Library chemistry: PCR-free (cohort 1) vs KAPA HyperPlus (cohort 2)')
print('=' * 78)
b = qc[qc['aligner']=='BWA-MEM']
for c in ['cohort1','cohort2']:
    g = b[b['cohort']==c]
    print('  %-9s dup %5.1f %%   on-target %5.1f %%   insert %5.0f bp   CYP2D6 depth %7.0fx' % (
        c, g['pct_dup_on_target'].median(), g['pct_on_target'].median(),
        g['insert_median'].median(), g['depth_CYP2D6'].median()))
d1 = b.loc[b['cohort']=='cohort1','pct_dup_on_target'].median()
d2 = b.loc[b['cohort']=='cohort2','pct_dup_on_target'].median()
print('\n  duplicate-rate difference: %+.1f percentage points (cohort 1 - cohort 2)' % (d1-d2))
print('  Draft spot-check claimed cohort 1 88.9-91.5 %, cohort 2 65.6-72.3 %.')
print('  -> %s' % ('DIRECTION CONFIRMED' if d1 > d2 else
                   '*** DIRECTION REVERSED - the draft spot check was wrong ***'))

print()
print('  cohort 1 by sequencing run (BWA-MEM arm):')
for run, g in b[b['cohort']=='cohort1'].groupby('run'):
    print('    %-13s n=%-3d dup %5.1f %%   on-target %5.1f %%   insert %5.0f   CYP2D6 depth %7.0fx' % (
        run, len(g), g['pct_dup_on_target'].median(), g['pct_on_target'].median(),
        g['insert_median'].median(), g['depth_CYP2D6'].median()))

if not os.path.exists(CALLS):
    print('\n  (%s not found - skipping section 3)' % CALLS); sys.exit(0)

print()
print('=' * 78)
print('  3. Do the miscalled samples have worse sequencing metrics?')
print('=' * 78)
S = pd.read_csv(CALLS)
S = S[S['aligner']=='BWA-MEM'][['sample','cohort','stratum','gold','cyri_call','cyri_strict']]
m = S.merge(b[['sample','cohort','run','pct_dup_on_target','pct_on_target',
               'insert_median','depth_CYP2D6']], on=['sample','cohort'], how='left')

for c in ['cohort1','cohort2']:
    g = m[m['cohort']==c]
    ok, bad_ = g[g['cyri_strict']], g[~g['cyri_strict']]
    print('\n  --- %s: %d correct, %d wrong ---' % (c, len(ok), len(bad_)))
    for col, lab in [('depth_CYP2D6','CYP2D6 depth'), ('pct_dup_on_target','duplicate %'),
                     ('pct_on_target','on-target %'), ('insert_median','insert size')]:
        print('     %-14s correct %8.1f   wrong %8.1f' % (lab, ok[col].median(), bad_[col].median()))
    if len(bad_):
        print('\n     each miscalled sample, with its percentile rank in the cohort:')
        print('     %-9s %-22s %-13s %8s %8s %8s' % ('sample','gold','run','depth','pct','dup%'))
        for _, r in bad_.sort_values('depth_CYP2D6').iterrows():
            pr = (g['depth_CYP2D6'] < r['depth_CYP2D6']).mean() * 100
            print('     %-9s %-22s %-13s %8.0f %7.0f%% %8.1f' % (
                r['sample'], str(r['gold'])[:22], r['run'],
                r['depth_CYP2D6'] if pd.notna(r['depth_CYP2D6']) else -1, pr,
                r['pct_dup_on_target'] if pd.notna(r['pct_dup_on_target']) else -1))

print()
print('  cohort 1 concordance by sequencing run (batch hypothesis):')
g1 = m[m['cohort']=='cohort1']
for run, g in g1.groupby('run'):
    print('    %-13s %2d/%2d = %5.1f %%' % (run, g['cyri_strict'].sum(), len(g),
                                            100*g['cyri_strict'].mean()))
