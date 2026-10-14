#!/usr/bin/env python3
"""
recompute_concordance.py
========================
Score benchmark concordance against the reference-standard diplotypes in benchmark/data/,
using compare_genotype (with x/+ notation equivalence).

Does NOT re-run CyriPanel. Reads each iter's iteration_results.json, looks up the
called_genotype, compares it with the reference standard, and writes aggregated statistics.
comparison_single and concordance_single count a call as concordant only when it is a
single diplotype equal to the reference standard (single-candidate scoring). The
manuscript's primary, candidate-inclusive concordance is computed by recalculate.py.

Output:
  - Per-iteration summary CSV (one row per cohort/panel/iter)
  - Per-sample records CSV
  - Concordance, PASS-rate and per-sample tables (text)

Usage:
  python3 recompute_concordance.py
"""
import os
import sys
import json
import glob
import csv
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

# Set up paths and imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import benchmark_config as cfg
from run_single_iteration import compare_genotype, normalize_genotype

# Reference-standard diplotypes
_GOLD_DIR = os.path.join(SCRIPT_DIR, 'data')   # gold standards ship with the repo
COHORT_GOLD = {
    'cohort1_dragen': os.path.join(_GOLD_DIR, 'CYP2D6_cohort1_gold_standard_v2.csv'),
    'cohort2_dragen': os.path.join(_GOLD_DIR, 'CYP2D6_cohort2_gold_standard_v2.csv'),
    'cohort1': os.path.join(_GOLD_DIR, 'CYP2D6_cohort1_gold_standard_v2.csv'),
    'cohort2': os.path.join(_GOLD_DIR, 'CYP2D6_cohort2_gold_standard_v2.csv'),
}

OUT_DIR = '/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/recompute_v2'
os.makedirs(OUT_DIR, exist_ok=True)


def load_gold(csv_path):
    """Load gold standard CSV -> dict[sample_id] = genotype."""
    gold = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gold[row['Sample_ID']] = row['Genotype']
    return gold


def recompute_iteration(json_path, gold_map):
    """
    Read one iteration_results.json, compare each call with the reference standard,
    and return the summary stats + the per-sample records.
    """
    with open(json_path) as f:
        data = json.load(f)
    
    per_sample = data.get('per_sample', {})
    summary_old = data.get('summary', {})
    
    n_tested_with_gold = 0
    n_concordant = 0
    n_discordant = 0
    n_no_call = 0
    n_pass = 0
    
    sample_records = []
    for sid, ps in per_sample.items():
        called = ps.get('called_genotype')
        gold_gt = gold_map.get(sid)
        comp = compare_genotype(called, gold_gt)
        filt = ps.get('filter')
        
        # Track summary
        if gold_gt is not None:
            n_tested_with_gold += 1
        if comp == 'concordant':
            n_concordant += 1
        if comp == 'discordant':
            n_discordant += 1
        if called is None or comp == 'no_call':
            n_no_call += 1
        if filt == 'PASS':
            n_pass += 1
        
        sample_records.append({
            'sample_id': sid,
            'called_genotype': called,
            'gold': gold_gt,
            'comparison_single': comp,
            'filter': filt,
        })
    
    n_total = len([s for s in sample_records if s['called_genotype'] is not None or s['gold'] is not None])
    
    summary = {
        'cohort':              summary_old.get('cohort'),
        'panel_size':          summary_old.get('panel_size'),
        'iteration':           summary_old.get('iteration'),
        'seed':                summary_old.get('seed'),
        'n_test_samples':      summary_old.get('n_test_samples'),
        'n_tested_with_gold':  n_tested_with_gold,
        # Single-candidate scoring
        'n_concordant_single':     n_concordant,
        'n_discordant_single':     n_discordant,
        'concordance_single':      n_concordant / n_tested_with_gold if n_tested_with_gold else None,
        # Shared
        'n_no_call':           n_no_call,
        'n_pass':              n_pass,
        'pass_rate':           n_pass / n_tested_with_gold if n_tested_with_gold else None,
    }
    return summary, sample_records


def main():
    all_summaries = []
    all_sample_records = []
    
    for cohort in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
        gold_map = load_gold(COHORT_GOLD[cohort])
        print(f'Loaded the reference standard for {cohort}: {len(gold_map)} samples')
        
        for panel in [5, 10, 20]:
            base = os.path.join(cfg.BENCHMARK_ROOT, cohort, f'panel{panel}')
            files = sorted(glob.glob(os.path.join(base, 'iter_*/iteration_results.json')))
            print(f'  {cohort} panel={panel}: {len(files)} iterations')
            
            for f in files:
                summary, recs = recompute_iteration(f, gold_map)
                summary['cohort'] = cohort
                summary['panel_size'] = panel
                all_summaries.append(summary)
                for r in recs:
                    r['cohort'] = cohort
                    r['panel_size'] = panel
                    r['iteration'] = summary['iteration']
                    all_sample_records.append(r)
    
    df_sum = pd.DataFrame(all_summaries)
    df_samp = pd.DataFrame(all_sample_records)
    
    # Save raw recomputed data
    df_sum.to_csv(os.path.join(OUT_DIR, 'iteration_summary_v2.csv'), index=False)
    df_samp.to_csv(os.path.join(OUT_DIR, 'per_sample_v2.csv'), index=False)
    print()
    print(f'Saved iteration_summary_v2.csv ({len(df_sum)} rows) and per_sample_v2.csv ({len(df_samp)} rows)')
    print()
    
    # ─── Print single-candidate concordance ───────────────────────────
    print('=' * 110)
    print('  Single-candidate concordance (mean ± SD across iterations)')
    print('=' * 110)
    print()
    print('  %-16s %-8s | %-22s | %-12s' % ('Cohort', 'Panel', 'Concordance', '#iter'))
    print('  ' + '-' * 70)

    for cohort in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
        for panel in [5, 10, 20]:
            sub = df_sum[(df_sum['cohort'] == cohort) & (df_sum['panel_size'] == panel)]
            if sub.empty:
                continue
            print('  %-16s n=%-6d | %.3f ± %.3f        | %d' % (
                cohort, panel, sub['concordance_single'].mean(), sub['concordance_single'].std(), len(sub)))
    
    print()
    print('=' * 110)
    print('  PASS rate (independent of the reference standard)')
    print('=' * 110)
    print()
    for cohort in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
        for panel in [5, 10, 20]:
            sub = df_sum[(df_sum['cohort'] == cohort) & (df_sum['panel_size'] == panel)]
            if sub.empty:
                continue
            pm = sub['pass_rate'].mean()
            ps = sub['pass_rate'].std()
            ncm = sub['n_no_call'].mean()
            print('  %-10s n=%-6d  PASS rate: %.3f ± %.3f   No-call (mean count): %.2f' % (
                cohort, panel, pm, ps, ncm))
    
    # ─── Per-sample concordance under v2 ──────────────────────────────
    print()
    print('=' * 110)
    print('  Per-sample single-candidate concordance (samples with <100% concordance only)')
    print('=' * 110)
    
    for cohort in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
        print()
        print(f'--- {cohort.upper()} ---')
        df_c = df_samp[df_samp['cohort'] == cohort]
        # Group by sample
        per_sample = defaultdict(lambda: defaultdict(lambda: {'n_test': 0, 'n_concord': 0, 'gold': None}))
        for _, row in df_c.iterrows():
            sid = row['sample_id']
            panel = row['panel_size']
            if pd.notna(row['gold']):
                per_sample[sid][panel]['n_test'] += 1
                per_sample[sid][panel]['gold'] = row['gold']
                if row['comparison_single'] == 'concordant':
                    per_sample[sid][panel]['n_concord'] += 1
        
        # Print problem samples (any panel <100%)
        problem = []
        for sid, panels in per_sample.items():
            for p, st in panels.items():
                if st['n_test'] > 0 and st['n_concord'] < st['n_test']:
                    problem.append(sid)
                    break
        
        print('  %-12s %-25s %-15s %-15s %-15s' % ('Sample', 'Reference', 'panel=5', 'panel=10', 'panel=20'))
        print('  ' + '-' * 80)
        for sid in sorted(set(problem)):
            cells = []
            gold = ''
            for p in [5, 10, 20]:
                st = per_sample[sid].get(p, {'n_test': 0, 'n_concord': 0, 'gold': None})
                if st['n_test'] == 0:
                    cells.append('—')
                    continue
                gold = st['gold']
                rate = st['n_concord'] / st['n_test'] * 100
                cells.append('%d/%d (%d%%)' % (st['n_concord'], st['n_test'], rate))
            print('  %-12s %-25s %-15s %-15s %-15s' % (sid, gold or '?', cells[0], cells[1], cells[2]))


if __name__ == '__main__':
    main()
