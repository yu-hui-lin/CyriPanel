#!/usr/bin/env python3
"""
recompute_concordance.py
========================
Re-evaluate benchmark concordance using v2 gold standard (with *185 -> *4 corrections
from collaborator's decypher v2) and the patched compare_genotype function (with x/+
notation equivalence).

Does NOT re-run CyriPanel. Reads each iter's iteration_results.json, looks up the
called_genotype, re-compares against v2 gold, and writes new aggregated statistics.

Output:
  - Per-iteration recomputed summary CSV (one row per cohort/panel/iter)
  - Per-sample recomputed concordance rate CSV
  - Side-by-side v1 vs v2 statistics table (text)

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
SCRIPT_DIR = '/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark'
sys.path.insert(0, SCRIPT_DIR)
import benchmark_config as cfg
from run_single_iteration import compare_genotype, normalize_genotype

# v2 gold standards
COHORT_V2_GOLD = {
    'cohort1_dragen': '/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_gold_standard_v2.csv',
    'cohort2_dragen': '/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_gold_standard_v2.csv',
    'cohort1': '/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_gold_standard_v2.csv',
    'cohort2': '/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_gold_standard_v2.csv',
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


def recompute_iteration(json_path, v2_gold):
    """
    Read one iteration_results.json, recompute per-sample comparison against v2 gold,
    and return new summary stats + the per-sample dict with v1 / v2 comparison.
    """
    with open(json_path) as f:
        data = json.load(f)
    
    per_sample = data.get('per_sample', {})
    summary_old = data.get('summary', {})
    
    n_tested_with_gold = 0
    n_concordant_v1 = 0
    n_concordant_v2 = 0
    n_discordant_v1 = 0
    n_discordant_v2 = 0
    n_no_call = 0
    n_pass = 0
    
    sample_records = []
    for sid, ps in per_sample.items():
        called = ps.get('called_genotype')
        gold_v1 = ps.get('gold_genotype')
        gold_v2 = v2_gold.get(sid)
        comp_v1 = ps.get('comparison')  # original
        comp_v2 = compare_genotype(called, gold_v2)
        filt = ps.get('filter')
        
        # Track summary
        if gold_v2 is not None and called is not None:
            n_tested_with_gold += 1
        if comp_v2 == 'concordant':
            n_concordant_v2 += 1
        if comp_v2 == 'discordant':
            n_discordant_v2 += 1
        if comp_v1 == 'concordant':
            n_concordant_v1 += 1
        if comp_v1 == 'discordant':
            n_discordant_v1 += 1
        if called is None or comp_v2 == 'no_call':
            n_no_call += 1
        if filt == 'PASS':
            n_pass += 1
        
        sample_records.append({
            'sample_id': sid,
            'called_genotype': called,
            'gold_v1': gold_v1,
            'gold_v2': gold_v2,
            'comparison_v1': comp_v1,
            'comparison_v2': comp_v2,
            'filter': filt,
        })
    
    n_total = len([s for s in sample_records if s['called_genotype'] is not None or s['gold_v2'] is not None])
    
    summary = {
        'cohort':              summary_old.get('cohort'),
        'panel_size':          summary_old.get('panel_size'),
        'iteration':           summary_old.get('iteration'),
        'seed':                summary_old.get('seed'),
        'n_test_samples':      summary_old.get('n_test_samples'),
        'n_tested_with_gold':  n_tested_with_gold,
        # v1
        'n_concordant_v1':     n_concordant_v1,
        'n_discordant_v1':     n_discordant_v1,
        'concordance_v1':      n_concordant_v1 / n_tested_with_gold if n_tested_with_gold else None,
        # v2
        'n_concordant_v2':     n_concordant_v2,
        'n_discordant_v2':     n_discordant_v2,
        'concordance_v2':      n_concordant_v2 / n_tested_with_gold if n_tested_with_gold else None,
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
        v2_gold = load_gold(COHORT_V2_GOLD[cohort])
        print(f'Loaded v2 gold for {cohort}: {len(v2_gold)} samples')
        
        for panel in [5, 10, 20]:
            base = os.path.join(cfg.BENCHMARK_ROOT, cohort, f'panel{panel}')
            files = sorted(glob.glob(os.path.join(base, 'iter_*/iteration_results.json')))
            print(f'  {cohort} panel={panel}: {len(files)} iterations')
            
            for f in files:
                summary, recs = recompute_iteration(f, v2_gold)
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
    
    # ─── Print v1 vs v2 comparison ────────────────────────────────────
    print('=' * 110)
    print('  V1 vs V2 concordance comparison (mean ± SD across iterations)')
    print('=' * 110)
    print()
    print('  %-10s %-8s | %-22s | %-22s | %-22s | %-12s' % (
        'Cohort', 'Panel', 'v1 concordance', 'v2 concordance', 'Improvement', '#iter'))
    print('  ' + '-' * 100)
    
    for cohort in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
        for panel in [5, 10, 20]:
            sub = df_sum[(df_sum['cohort'] == cohort) & (df_sum['panel_size'] == panel)]
            if sub.empty:
                continue
            v1m = sub['concordance_v1'].mean()
            v1s = sub['concordance_v1'].std()
            v2m = sub['concordance_v2'].mean()
            v2s = sub['concordance_v2'].std()
            improvement = (v2m - v1m) * 100
            print('  %-10s n=%-6d | %.3f ± %.3f       | %.3f ± %.3f       | %+.2f percentage pts   | %d' % (
                cohort, panel, v1m, v1s, v2m, v2s, improvement, len(sub)))
    
    print()
    print('=' * 110)
    print('  PASS rate (unchanged from v1, as filter is independent of gold)')
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
    print('  Per-sample concordance under v2 gold (samples with <100% concordance only)')
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
            if pd.notna(row['gold_v2']) and row['called_genotype']:
                per_sample[sid][panel]['n_test'] += 1
                per_sample[sid][panel]['gold'] = row['gold_v2']
                if row['comparison_v2'] == 'concordant':
                    per_sample[sid][panel]['n_concord'] += 1
        
        # Print problem samples (any panel <100%)
        problem = []
        for sid, panels in per_sample.items():
            for p, st in panels.items():
                if st['n_test'] > 0 and st['n_concord'] < st['n_test']:
                    problem.append(sid)
                    break
        
        print('  %-12s %-25s %-15s %-15s %-15s' % ('Sample', 'v2 gold', 'panel=5', 'panel=10', 'panel=20'))
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
