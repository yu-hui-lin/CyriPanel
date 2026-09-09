#!/usr/bin/env python3
"""
regenerate_table3.py
====================
Regenerate Table 3 (3a/3b/3c) and Supplementary Tables S1/S3 for the CyriPanel
manuscript, from data, for BOTH tools, under ONE documented aggregation rule.

Why this exists
---------------
The hand-assembled draft mixed aggregation conventions between the strict and
lenient columns (strict 99/116 does not reproduce under any single rule; the
modal-call rule gives 106/115).  Every number below therefore comes from one
pass over the raw per-evaluation records, so 3a and 3b differ ONLY in the
scoring rule and in nothing else.

Aggregation rule (applies identically to CyriPanel and Aldy4)
-------------------------------------------------------------
For each (sample x tool x aligner) the MODAL call is taken across every
reference-panel configuration in that arm, with "no call" treated as its own
category.  Ties are broken deterministically (highest count, then
lexicographically smallest).  A sample whose modal outcome is "no call" is
scored as a failure for that tool and stays in the denominator.  This is what
keeps HG00232 (CyriPanel never produced a genotype) counted as a CyriPanel
failure rather than silently dropped -- dropping it would bias the paired
comparison toward CyriPanel.
A sensitivity mode (--modal-exclude-nocall) recomputes the modal call over
non-null evaluations only; the script reports whether any sample changes.

Scoring rules
-------------
strict  : the reported call must equal the gold genotype after normalisation.
          A multi-solution call ("A;B") is therefore always discordant.
lenient : additionally concordant when the gold genotype is one of the
          reported candidates.  Rationale (Methods): both tools derive the same
          candidate set at the major-allele stage -- Aldy4's ILP returns
          solutions with identical objective values and selects one at the
          later minor-allele stage -- so scoring a reported ambiguity as
          discordant while scoring a suppressed one as concordant measures
          reporting convention, not inferential accuracy.

Statistics
----------
exact McNemar : two-sided binomial test on the discordant pairs, p = 0.5.
95% CI on the paired risk difference : Newcombe (1998) Stat Med 17:2635-2650,
          method 10 ("square-and-add").  Verified by --selftest against the
          worked value in the draft (a=19,b=7,c=0,d=2 -> +8.7%, +42.1%).
No scipy / numpy dependency: exact binomial via math.comb, z hard-coded.

Usage
-----
  python3 regenerate_table3.py                 # writes CSVs + markdown
  python3 regenerate_table3.py --selftest      # check the CI implementation
  python3 regenerate_table3.py --modal-exclude-nocall
"""

import argparse
import collections
import csv
import math
import os
import re
import sys

import pandas as pd

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
ROOT       = '/work/u7715055/staging/biology/u7715055'
CYRI_DIR   = os.path.join(ROOT, 'CyriPanel')
PER_SAMPLE = os.path.join(ROOT, 'cyripanel_benchmark_results/recompute_v2/per_sample_v2.csv')
ALDY_LONG  = os.path.join(ROOT, 'aldy_bench/aldy_results_long.csv')
GOLD = {
    'cohort1': os.path.join(CYRI_DIR, 'benchmark/data/CYP2D6_cohort1_gold_standard_v2.csv'),
    'cohort2': os.path.join(CYRI_DIR, 'benchmark/data/CYP2D6_cohort2_gold_standard_v2.csv'),
}
OUT_DIR = os.path.join(ROOT, 'cyripanel_benchmark_results/recompute_v2/table3')

Z95 = 1.959963984540054

# ----------------------------------------------------------------------
# Genotype normalisation -- CyriPanel's own functions, with a local fallback
# ----------------------------------------------------------------------
sys.path.insert(0, os.path.join(CYRI_DIR, 'benchmark'))
try:
    from run_single_iteration import normalize_genotype as _norm_gt
    NORM_SOURCE = 'CyriPanel run_single_iteration.normalize_genotype'
except Exception as exc:                                    # pragma: no cover
    NORM_SOURCE = 'LOCAL FALLBACK (%s)' % exc

    def _sort_key(tok):
        nums = re.findall(r'(\d+)', tok)
        return tuple(int(n) for n in nums) if nums else (tok,)

    def _expand(tok):
        m = re.match(r'^(\*\d+(?:\.\d+)?)x(\d+)$', tok.strip())
        return [m.group(1)] * int(m.group(2)) if m else [tok.strip()]

    def _norm_hap(h):
        parts = []
        for s_ in (h.split('+') if '+' in h else [h]):
            parts.extend(_expand(s_))
        return parts[0] if len(parts) == 1 else '+'.join(sorted(parts, key=_sort_key))

    def _norm_gt(gt):
        if gt is None or str(gt) in ('None', '', 'nan', 'Error'):
            return None
        gt = str(gt).strip()
        if '/' in gt:
            haps = sorted((_norm_hap(p) for p in gt.split('/')), key=_sort_key)
            return '/'.join(haps)
        return _norm_hap(gt)


def norm(gt):
    """Normalise one diplotype; returns None for absent / unparsable."""
    if gt is None:
        return None
    s = str(gt).strip()
    if s in ('', 'nan', 'NaN', 'None', 'Error', 'NA'):
        return None
    try:
        return _norm_gt(s)
    except Exception:
        return s


def candidates(call):
    """Split a possibly multi-solution call into normalised candidates."""
    if call is None:
        return []
    return [c for c in (norm(x) for x in str(call).split(';')) if c]


# ----------------------------------------------------------------------
# Structural strata (definition as published in Table 1)
# ----------------------------------------------------------------------
HYBRID_CORES = {13, 36, 68, 83}
_TOK = re.compile(r'\*(\d+)(?:\.\d+)?(?:x(\d+))?')


def core_tokens(gt):
    """[(core allele number, copy multiplier), ...] over both haplotypes."""
    out = []
    for hap in str(gt).split('/'):
        for tok in hap.split('+'):
            m = _TOK.match(tok.strip())
            if m:
                out.append((int(m.group(1)), int(m.group(2) or 1)))
    return out


def stratum(gold):
    """no_SV | CNV_only | hybrid, assigned from the v2 gold diplotype."""
    toks = core_tokens(gold)
    if any(c in HYBRID_CORES for c, _ in toks):
        return 'hybrid'
    if any(c == 5 for c, _ in toks):                       # whole-gene deletion
        return 'CNV_only'
    total_cn = sum(m for c, m in toks if c != 5)
    if total_cn != 2:
        return 'CNV_only'
    return 'no_SV'


# ----------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------
def wilson(k, n, z=Z95):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    den = 1.0 + z * z / n
    ctr = p + z * z / (2 * n)
    hlf = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((ctr - hlf) / den, (ctr + hlf) / den)


def newcombe10(a, b, c, d, z=Z95):
    """95% CI for p1 - p2 with paired data. a=both+, b=1 only, c=2 only, d=both-."""
    n = a + b + c + d
    if n == 0:
        return (float('nan'), float('nan'))
    p1, p2 = (a + b) / n, (a + c) / n
    theta = p1 - p2
    l1, u1 = wilson(a + b, n, z)
    l2, u2 = wilson(a + c, n, z)
    A = (a + b) * (c + d) * (a + c) * (b + d)
    phi = 0.0 if A == 0 else (a * d - b * c) / math.sqrt(A)
    t1, t2 = p1 - l1, u2 - p2
    lo = theta - math.sqrt(max(0.0, t1 * t1 - 2 * phi * t1 * t2 + t2 * t2))
    t3, t4 = u1 - p1, p2 - l2
    hi = theta + math.sqrt(max(0.0, t3 * t3 - 2 * phi * t3 * t4 + t4 * t4))
    return (max(-1.0, lo), min(1.0, hi))


def mcnemar_exact(b, c):
    """Two-sided exact McNemar p-value on the discordant pairs."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


def selftest():
    lo, hi = newcombe10(19, 7, 0, 2)
    p = mcnemar_exact(7, 0)
    ok = (abs(lo - 0.0873) < 5e-4 and abs(hi - 0.4210) < 5e-4
          and abs(p - 0.015625) < 1e-9)
    print('Newcombe method 10, a=19 b=7 c=0 d=2 -> RD %+.1f%% (%+.1f%%, %+.1f%%)'
          % (100 * (7 - 0) / 28, 100 * lo, 100 * hi))
    print('exact McNemar b=7 c=0 -> p = %.6f' % p)
    print('expected from draft: +25.0%% (+8.7%%, +42.1%%), p = 0.0156')
    print('SELFTEST %s' % ('PASS' if ok else 'FAIL'))
    # stratum assignment spot checks
    cases = [('*1/*4', 'no_SV'), ('*27/*2', 'no_SV'), ('*1/*162', 'no_SV'),
             ('*5/*4+*4', 'CNV_only'), ('*10/*10+*10', 'CNV_only'),
             ('*1x2/*4', 'CNV_only'), ('*10+*36+*36/*71', 'hybrid'),
             ('*2/*4+*68+*68', 'hybrid'), ('*35/*4+*68', 'hybrid')]
    bad = [(g, stratum(g), e) for g, e in cases if stratum(g) != e]
    print('stratum spot checks: %s' % ('PASS' if not bad else 'FAIL %s' % bad))
    return 0 if (ok and not bad) else 1


# ----------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------
def load_gold():
    gold = {}
    cohort_of = {}
    for ch, path in GOLD.items():
        with open(path) as f:
            for row in csv.DictReader(f):
                sid = row['Sample_ID'].strip()
                gold[sid] = row['Genotype'].strip()
                cohort_of[sid] = ch
    return gold, cohort_of


def modal(values, include_nocall=True):
    """Deterministic modal value. values: iterable of str-or-None."""
    vals = list(values)
    if not include_nocall:
        vals = [v for v in vals if v is not None]
    if not vals:
        return None, 0, 0
    cnt = collections.Counter('\x00NOCALL' if v is None else v for v in vals)
    top = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0]
    call = None if top[0] == '\x00NOCALL' else top[0]
    return call, top[1], len(vals)


def load_cyripanel(include_nocall):
    df = pd.read_csv(PER_SAMPLE)
    df['aligner'] = df['cohort'].map(
        lambda c: 'DRAGEN' if str(c).endswith('_dragen') else 'BWA-MEM')
    rows = []
    for (aln, sid), g in df.groupby(['aligner', 'sample_id']):
        call, n_top, n_ev = modal(
            [None if pd.isna(v) else str(v) for v in g['called_genotype']],
            include_nocall)
        n_multi = int(g['called_genotype'].astype(str).str.contains(';', na=False).sum())
        rows.append(dict(tool='CyriPanel', aligner=aln, sample=sid, modal_call=call,
                         n_events=n_ev, modal_count=n_top, n_multi_events=n_multi))
    return pd.DataFrame(rows)


def load_aldy(include_nocall):
    df = pd.read_csv(ALDY_LONG)
    df['aligner'] = df['aligner'].astype(str).str.upper().map(
        lambda a: 'DRAGEN' if 'DRAGEN' in a else 'BWA-MEM')
    rows = []
    for (aln, sid), g in df.groupby(['aligner', 'sample']):
        call, n_top, n_ev = modal(
            [None if pd.isna(v) else str(v) for v in g['normalised']],
            include_nocall)
        n_multi = int(g['normalised'].astype(str).str.contains(';', na=False).sum())
        rows.append(dict(tool='Aldy4', aligner=aln, sample=sid, modal_call=call,
                         n_events=n_ev, modal_count=n_top, n_multi_events=n_multi,
                         n_sol_max=int(pd.to_numeric(g.get('n_solutions'),
                                                     errors='coerce').fillna(1).max())))
    return pd.DataFrame(rows)


def score(call, gold_norm):
    """-> (strict_bool, lenient_bool, is_multi, is_nocall)"""
    if call is None:
        return False, False, False, True
    cands = candidates(call)
    multi = len(cands) > 1
    strict = (len(cands) == 1 and cands[0] == gold_norm)
    lenient = gold_norm in cands
    return strict, lenient, multi, False


# ----------------------------------------------------------------------
# Table assembly
# ----------------------------------------------------------------------
def paired(df, aligner, rule, subset=None):
    """Return (a, b, c, d, n) with tool1 = CyriPanel, tool2 = Aldy4."""
    d = df[df['aligner'] == aligner]
    if subset is not None:
        d = d[subset(d)]
    a = int(((d['cyri_' + rule]) & (d['aldy_' + rule])).sum())
    b = int(((d['cyri_' + rule]) & (~d['aldy_' + rule])).sum())
    c = int(((~d['cyri_' + rule]) & (d['aldy_' + rule])).sum())
    e = int(((~d['cyri_' + rule]) & (~d['aldy_' + rule])).sum())
    return a, b, c, e, len(d)


def row_stats(a, b, c, e):
    n = a + b + c + e
    if n == 0:
        return dict(n=0, cyripanel=float('nan'), aldy4=float('nan'), rd=float('nan'),
                    ci_lo=float('nan'), ci_hi=float('nan'), p=float('nan'),
                    c_only=0, a_only=0, both=0, neither=0)
    p1, p2 = (a + b) / n, (a + c) / n
    lo, hi = newcombe10(a, b, c, e)
    return dict(n=n, cyripanel=100 * p1, aldy4=100 * p2, rd=100 * (p1 - p2),
                ci_lo=100 * lo, ci_hi=100 * hi, p=mcnemar_exact(b, c),
                c_only=b, a_only=c, both=a, neither=e)


def fmt_row(label_cells, s):
    ci = '(%+.1f, %+.1f)' % (s['ci_lo'], s['ci_hi']) if s['c_only'] + s['a_only'] else '—'
    return '| %s | %d | %.1f %% | %.1f %% | %+.1f | %s | %.4f |' % (
        ' | '.join(label_cells), s['n'], s['cyripanel'], s['aldy4'], s['rd'], ci, s['p'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--modal-exclude-nocall', action='store_true',
                    help='sensitivity: take the modal call over non-null evaluations only')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--out', default=OUT_DIR)
    args = ap.parse_args()
    if args.selftest:
        sys.exit(selftest())

    inc = not args.modal_exclude_nocall
    os.makedirs(args.out, exist_ok=True)
    print('normalize_genotype source : %s' % NORM_SOURCE)
    print('modal rule                : %s\n' %
          ('no-call counts as a category (primary)' if inc
           else 'no-call excluded (sensitivity)'))

    gold_raw, cohort_of = load_gold()
    gold_n = {s: norm(g) for s, g in gold_raw.items()}
    strat = {s: stratum(g) for s, g in gold_raw.items()}
    print('gold standard             : %d samples' % len(gold_raw))
    sc = collections.Counter(strat.values())
    print('strata                    : no_SV %d, CNV_only %d, hybrid %d  (expected 61/27/28)'
          % (sc['no_SV'], sc['CNV_only'], sc['hybrid']))

    cyri = load_cyripanel(inc)
    aldy = load_aldy(inc)
    print('CyriPanel arms            : %s' % sorted(cyri['aligner'].unique()))
    print('Aldy4 arms                : %s' % sorted(aldy['aligner'].unique()))
    if 'n_sol_max' in aldy:
        print('Aldy4 max n_solutions     : %d  (multi-solution rows: %d)'
              % (aldy['n_sol_max'].max(), int(aldy['n_multi_events'].sum())))

    # ---- merge into one sample-level frame -------------------------------
    recs = []
    for aln in ['BWA-MEM', 'DRAGEN']:
        cy = cyri[cyri['aligner'] == aln].set_index('sample')
        al = aldy[aldy['aligner'] == aln].set_index('sample')
        both = sorted(set(cy.index) | set(al.index))
        for sid in both:
            if sid not in gold_n:
                continue
            g = gold_n[sid]
            cc = cy['modal_call'].get(sid) if sid in cy.index else None
            ac = al['modal_call'].get(sid) if sid in al.index else None
            cc = None if (cc is None or pd.isna(cc)) else cc
            ac = None if (ac is None or pd.isna(ac)) else ac
            cs, cl, cm, cn = score(cc, g)
            as_, al_, am, an = score(ac, g)
            recs.append(dict(
                aligner=aln, sample=sid, cohort=cohort_of[sid],
                stratum=strat[sid], gold=gold_raw[sid], gold_norm=g,
                cyri_call=cc, aldy_call=ac,
                cyri_strict=cs, cyri_lenient=cl, cyri_multi=cm, cyri_nocall=cn,
                aldy_strict=as_, aldy_lenient=al_, aldy_multi=am, aldy_nocall=an,
                cyri_events=int(cy['n_events'].get(sid, 0)) if sid in cy.index else 0,
                cyri_multi_events=int(cy['n_multi_events'].get(sid, 0)) if sid in cy.index else 0,
                aldy_events=int(al['n_events'].get(sid, 0)) if sid in al.index else 0,
            ))
    S = pd.DataFrame(recs)
    S.to_csv(os.path.join(args.out, 'sample_level_calls.csv'), index=False)

    for aln in ['BWA-MEM', 'DRAGEN']:
        d = S[S['aligner'] == aln]
        miss_c = sorted(d.loc[d['cyri_call'].isna(), 'sample'])
        miss_a = sorted(d.loc[d['aldy_call'].isna(), 'sample'])
        print('\n%-8s n = %d   CyriPanel no-call: %s   Aldy4 no-call: %s'
              % (aln, len(d), miss_c or 'none', miss_a or 'none'))
        only_c = sorted(set(cyri[cyri.aligner == aln]['sample']) - set(aldy[aldy.aligner == aln]['sample']))
        only_a = sorted(set(aldy[aldy.aligner == aln]['sample']) - set(cyri[cyri.aligner == aln]['sample']))
        if only_c or only_a:
            print('  WARNING evaluated by one tool only  CyriPanel-only: %s  Aldy-only: %s'
                  % (only_c, only_a))

    # ---- Tables 3a / 3b ---------------------------------------------------
    md = ['# Table 3 (regenerated from data)', '',
          'Aggregation: modal call per sample per tool per aligner across all '
          'reference-panel configurations; no-call treated as a category and '
          'scored as a failure. Identical rule for both tools and both scoring '
          'rules, so 3a and 3b differ only in the scoring rule.', '']
    out_rows = []
    for rule, title in [('lenient', '3a. Lenient scoring'),
                        ('strict', '3b. Strict scoring (sensitivity analysis)')]:
        md += ['### %s' % title, '',
               '| Aligner | Stratum | n | CyriPanel | Aldy4 | RD | 95 % CI | p |',
               '|---|---|---|---|---|---|---|---|']
        for aln in ['BWA-MEM', 'DRAGEN']:
            for st in ['no_SV', 'CNV_only', 'hybrid']:
                a, b, c, e, n = paired(S, aln, rule, lambda d, st=st: d['stratum'] == st)
                s = row_stats(a, b, c, e)
                md.append(fmt_row([aln, '`%s`' % st], s))
                out_rows.append(dict(rule=rule, aligner=aln, stratum=st, **s))
            a, b, c, e, n = paired(S, aln, rule)
            s = row_stats(a, b, c, e)
            md.append(fmt_row([aln, '*all*'], s))
            out_rows.append(dict(rule=rule, aligner=aln, stratum='all', **s))
        md.append('')
    pd.DataFrame(out_rows).to_csv(os.path.join(args.out, 'table3ab.csv'), index=False)

    # ---- Table 3c per cohort ---------------------------------------------
    md += ['### 3c. Per cohort (strict scoring)', '',
           '| Cohort | Aligner | n | CyriPanel | Aldy4 | C-only | A-only | p |',
           '|---|---|---|---|---|---|---|---|']
    c_rows = []
    for ch in ['cohort1', 'cohort2']:
        for aln in ['BWA-MEM', 'DRAGEN']:
            a, b, c, e, n = paired(S, aln, 'strict', lambda d, ch=ch: d['cohort'] == ch)
            s = row_stats(a, b, c, e)
            md.append('| %s | %s | %d | %.1f %% | %.1f %% | %d | %d | %.4f |'
                      % (ch[-1], aln, s['n'], s['cyripanel'], s['aldy4'],
                         s['c_only'], s['a_only'], s['p']))
            c_rows.append(dict(cohort=ch, aligner=aln, **s))
    pd.DataFrame(c_rows).to_csv(os.path.join(args.out, 'table3c.csv'), index=False)
    md.append('')

    # ---- S1 multi-solution ------------------------------------------------
    flip = S[(S['cyri_lenient']) & (~S['cyri_strict'])]
    multi = S[S['cyri_multi']]
    md += ['### Supplementary Table S1 (regenerated). CyriPanel multi-solution calls', '',
           '| Aligner | Sample | Gold | Modal candidates | Stratum | multi events / total | gold among candidates |',
           '|---|---|---|---|---|---|---|']
    for _, r in multi.sort_values(['aligner', 'sample']).iterrows():
        md.append('| %s | %s | `%s` | `%s` | `%s` | %d / %d | %s |'
                  % (r['aligner'], r['sample'], r['gold'], r['cyri_call'],
                     r['stratum'], r['cyri_multi_events'], r['cyri_events'],
                     'yes' if r['cyri_lenient'] else 'NO'))
    md.append('')
    multi.to_csv(os.path.join(args.out, 'tableS1_multisolution.csv'), index=False)
    flip.to_csv(os.path.join(args.out, 'strict_lenient_flips.csv'), index=False)

    # ---- S3 hybrid discordant --------------------------------------------
    md += ['### Supplementary Table S3 (regenerated). Hybrid stratum, one tool correct only', '']
    s3 = []
    for aln in ['BWA-MEM', 'DRAGEN']:
        for rule in ['lenient', 'strict']:
            h = S[(S['aligner'] == aln) & (S['stratum'] == 'hybrid')]
            conly = h[h['cyri_' + rule] & ~h['aldy_' + rule]]
            aonly = h[~h['cyri_' + rule] & h['aldy_' + rule]]
            md += ['**%s, %s** — CyriPanel only n = %d; Aldy4 only n = %d' %
                   (aln, rule, len(conly), len(aonly)), '',
                   '| Winner | Sample | Gold | CyriPanel call | Aldy4 call |',
                   '|---|---|---|---|---|']
            for tag, sub in [('CyriPanel', conly), ('Aldy4', aonly)]:
                for _, r in sub.sort_values('sample').iterrows():
                    md.append('| %s | %s | `%s` | `%s` | `%s` |'
                              % (tag, r['sample'], r['gold'], r['cyri_call'], r['aldy_call']))
                    s3.append(dict(aligner=aln, rule=rule, winner=tag, sample=r['sample'],
                                   gold=r['gold'], cyri=r['cyri_call'], aldy=r['aldy_call']))
            md.append('')
    pd.DataFrame(s3).to_csv(os.path.join(args.out, 'tableS3_hybrid.csv'), index=False)

    # ---- Table 2 PASS rate, one documented denominator (Open item 1) -------
    it_sum = os.path.join(os.path.dirname(PER_SAMPLE), 'iteration_summary_v2.csv')
    if os.path.exists(it_sum):
        isum = pd.read_csv(it_sum)
        md += ['### Table 2 PASS rate (denominator = n_tested_with_gold)', '',
               'PASS rate is the proportion of *evaluable* samples (those with both a '
               'call and a gold-standard genotype) whose Genotype filter was PASS; the '
               'no-call rate is reported separately. This is the same denominator used '
               'for concordance, so the two columns are directly comparable.', '',
               '| Cohort | Aligner | Panel size | Concordance | PASS rate | mean no-calls |',
               '|---|---|---|---|---|---|']
        p_rows = []
        for coh in ['cohort1', 'cohort2', 'cohort1_dragen', 'cohort2_dragen']:
            for panel in [5, 10, 20]:
                sub = isum[(isum['cohort'] == coh) & (isum['panel_size'] == panel)]
                if sub.empty:
                    continue
                aln = 'DRAGEN' if coh.endswith('_dragen') else 'BWA-MEM'
                md.append('| %s | %s | %d | %.1f ± %.1f %% | %.1f ± %.1f %% | %.2f |'
                          % (coh.replace('_dragen', '')[-1], aln, panel,
                             100 * sub['concordance_v2'].mean(), 100 * sub['concordance_v2'].std(),
                             100 * sub['pass_rate'].mean(), 100 * sub['pass_rate'].std(),
                             sub['n_no_call'].mean()))
                p_rows.append(dict(cohort=coh, aligner=aln, panel_size=panel,
                                   concordance_mean=sub['concordance_v2'].mean(),
                                   concordance_sd=sub['concordance_v2'].std(),
                                   pass_rate_mean=sub['pass_rate'].mean(),
                                   pass_rate_sd=sub['pass_rate'].std(),
                                   n_no_call_mean=sub['n_no_call'].mean(),
                                   n_iterations=len(sub)))
        md.append('')
        pd.DataFrame(p_rows).to_csv(os.path.join(args.out, 'table2_pass_rate.csv'), index=False)
    else:
        print('NOTE: %s not found -- Table 2 PASS-rate block skipped' % it_sum)

    md_path = os.path.join(args.out, 'table3_regenerated.md')
    with open(md_path, 'w') as f:
        f.write('\n'.join(md) + '\n')

    # ---- console summary ---------------------------------------------------
    print('\n' + '=' * 78)
    print('  Strict vs lenient: samples that change status (CyriPanel)')
    print('=' * 78)
    if flip.empty:
        print('  none')
    else:
        print(flip[['aligner', 'sample', 'gold', 'cyri_call', 'stratum']]
              .to_string(index=False))
    print('\n  hybrid stratum, CyriPanel:')
    for aln in ['BWA-MEM', 'DRAGEN']:
        h = S[(S['aligner'] == aln) & (S['stratum'] == 'hybrid')]
        print('    %-8s strict %d/%d   lenient %d/%d'
              % (aln, h['cyri_strict'].sum(), len(h), h['cyri_lenient'].sum(), len(h)))
    print('\n' + '\n'.join(l for l in md if l.startswith('|') or l.startswith('###')))
    print('\nWrote: %s' % args.out)
    for f_ in sorted(os.listdir(args.out)):
        print('   %s' % f_)


if __name__ == '__main__':
    main()

