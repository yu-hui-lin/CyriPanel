#!/usr/bin/env python3
"""Memory of the CyriPanel benchmark array tasks that produced the published results.
A. every array job, B. the task behind each final iteration_results.json,
C. final tasks above HIGH GiB, D. sample states and E. reference samples that
separate the high-memory iterations from the others."""
import os, json, time, subprocess, statistics as st
from collections import defaultdict

ROOT = os.environ.get('CYRI_RESULTS', '/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results')
HIGH = 6.0  # GiB
names = ['cyri_%s%s_n%d' % (c, a, n) for c in ('cohort1', 'cohort2') for a in ('', '_dragen') for n in (5, 10, 20)]
out = subprocess.run(['sacct', '-u', os.environ.get('USER', ''), '-S', '2026-01-01', '-P', '-n',
                      '--name=' + ','.join(names),
                      '--format=JobID,JobName,State,Start,End,ElapsedRaw,TotalCPU,MaxRSS'],
                     stdout=subprocess.PIPE, universal_newlines=True).stdout

U = {'K': 1 / 1024 ** 2, 'M': 1 / 1024, 'G': 1.0, 'T': 1024.0}
def gib(m):
    if not m: return None
    return float(m[:-1]) * U[m[-1]] if m[-1] in U else float(m) / 1024 ** 3
def secs(t):
    d = 0
    if '-' in t:
        d, t = t.split('-'); d = int(d)
    s = 0.0
    for x in t.split(':'): s = s * 60 + float(x)
    return d * 86400 + s
def epoch(ts):
    try: return time.mktime(time.strptime(ts, '%Y-%m-%dT%H:%M:%S'))
    except ValueError: return None
def q(x, p):
    x = sorted(x); i = (len(x) - 1) * p; lo = int(i); hi = min(lo + 1, len(x) - 1)
    return x[lo] + (x[hi] - x[lo]) * (i - lo)

jobname, steps = {}, []
for line in out.splitlines():
    f = line.split('|')
    if len(f) != 8: continue
    if f[0].endswith('.batch'): steps.append([f[0][:-6]] + f[2:])
    elif '.' not in f[0]: jobname[f[0]] = f[1]

rows = []
for parent, state, start, end, el, cpu, rss in steps:
    jn = jobname.get(parent, '')
    if not jn.startswith('cyri_') or '_' not in parent: continue
    arr, idx = parent.split('_', 1)
    if not idx.isdigit(): continue
    arm, _, n = jn[5:].rpartition('_n')
    path = '%s/%s/panel%s/iter_%03d/iteration_results.json' % (ROOT, arm, n, int(idx))
    mt = os.path.getmtime(path) if os.path.exists(path) else None
    e = epoch(end)
    rows.append(dict(arm=arm, n=int(n), it=int(idx), arr=arr, state=state, start=start[:10], end=end,
                     el=int(el or 0), cpu=secs(cpu) if cpu else 0.0, rss=gib(rss), path=path,
                     match=(mt is not None and e is not None and -120 <= e - mt <= 3600)))
ok = [r for r in rows if r['state'] == 'COMPLETED' and r['el'] >= 600 and r['rss'] is not None]

print('== A. Array jobs (completed tasks of at least 10 min)')
by_arr = defaultdict(list)
for r in ok: by_arr[(r['arm'], r['n'], r['arr'])].append(r)
for k in sorted(by_arr):
    v = by_arr[k]; g = [r['rss'] for r in v]
    print('%-15s n=%-2d job %-9s tasks=%-2d started %s to %s  final=%-2d  MaxRSS median %.2f max %.2f GiB  >%.0f GiB: %d'
          % (k[0], k[1], k[2], len(v), min(r['start'] for r in v), max(r['start'] for r in v),
             sum(r['match'] for r in v), st.median(g), max(g), HIGH, sum(x > HIGH for x in g)))

final = {}
for r in ok:
    k = (r['arm'], r['n'], r['it'])
    if k not in final or (r['match'], r['end']) > (final[k]['match'], final[k]['end']): final[k] = r
print('\n== B. Tasks behind the final iteration_results.json (matched by file time, else the latest run)')
grp = defaultdict(list)
for k, r in final.items(): grp[(k[0], k[1])].append(r)
for k in sorted(grp):
    v = grp[k]; g = [r['rss'] for r in v]; c = [r['cpu'] / r['el'] for r in v]
    print('%-15s n=%-2d iterations=%-2d (by file time %-2d)  MaxRSS median %.2f IQR %.2f-%.2f max %.2f GiB  >%.0f GiB: %-2d  CPU/wall %.2f-%.2f'
          % (k[0], k[1], len(v), sum(r['match'] for r in v), q(g, .5), q(g, .25), q(g, .75), max(g), HIGH,
             sum(x > HIGH for x in g), min(c), max(c)))

print('\n== C. Final tasks above %.0f GiB' % HIGH)
high = {k: r for k, r in final.items() if r['rss'] > HIGH}
for k in sorted(high):
    r = high[k]
    print('%-15s n=%-2d iter %02d  job %s_%d  MaxRSS %.1f GiB  elapsed %.0f min' % (k[0], k[1], k[2], r['arr'], k[2], r['rss'], r['el'] / 60))

print('\n== D/E. What the high-memory iterations share, per arm (h = high iterations with it, o = other iterations with it)')
for arm in sorted({k[0] for k in high}):
    H = [k for k in final if k[0] == arm and k in high]; O = [k for k in final if k[0] == arm and k not in high]
    state, ref = defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
    for keys, col in ((H, 0), (O, 1)):
        for k in keys:
            try: d = json.load(open(final[k]['path']))
            except Exception: continue
            for sid, ps in d.get('per_sample', {}).items():
                state[(sid, ps.get('total_cn'), ps.get('spacer_cn'), ps.get('cnv_group'), ps.get('filter'), ps.get('called_genotype'))][col] += 1
            for sid in d.get('summary', {}).get('reference_ids', []): ref[sid][col] += 1
    enrich = lambda kv: -(kv[1][0] / max(len(H), 1) - kv[1][1] / max(len(O), 1))
    print('-- %s: %d high, %d other iterations' % (arm, len(H), len(O)))
    print('   test-sample states (sample, total_cn, spacer_cn, cnv_group, filter, call):')
    for s_, (h, o) in sorted(state.items(), key=enrich)[:8]: print('     h=%-2d o=%-3d %s' % (h, o, s_))
    print('   samples in the reference panel:')
    for s_, (h, o) in sorted(ref.items(), key=enrich)[:5]: print('     h=%-2d o=%-3d %s' % (h, o, s_))
