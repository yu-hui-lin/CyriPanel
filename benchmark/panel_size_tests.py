#!/usr/bin/env python3
"""Friedman and Wilcoxon tests for reference-panel size, from the archived records.

Reads recalculate_output/common_sample_rates.csv, which recalculate.py writes, so run
recalculate.py first. For each cohort and aligner it uses the samples evaluated at every
panel size and tests their candidate-inclusive per-sample concordance at 5, 10 and 20
references, as Supplemental Methods S6 describes.

  * Friedman test across the three panel sizes (chi-square approximation). It is reported
    as undefined when no sample's concordance differs between panel sizes.
  * Pairwise two-sided Wilcoxon signed-rank tests (zero_method="wilcox"). A comparison
    with no non-zero paired difference gets p = 1.
  * Holm correction across the three pairwise comparisons.

The test logic is the same as paired_tests() in common_sample_analysis.py, which read the
per-iteration outputs on the cluster. This script needs scipy and writes
panel_size_tests.csv next to it.
"""
import csv
import itertools
import os
import statistics
import sys

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "recalculate_output", "common_sample_rates.csv")
OUT = os.path.join(HERE, "panel_size_tests.csv")
PANELS = (5, 10, 20)

if not os.path.exists(SRC):
    sys.exit(f"{SRC} not found. Run recalculate.py first.")

rates = {}
for r in csv.DictReader(open(SRC)):
    rates.setdefault(r["cohort"], {}).setdefault(r["sample_id"], {})[int(r["panel_size"])] = float(r["concordance_pct"])

rows = []
for cohort in sorted(rates):
    ids = sorted(s for s, v in rates[cohort].items() if all(p in v for p in PANELS))
    col = {p: [rates[cohort][s][p] for s in ids] for p in PANELS}
    varies = any(len({rates[cohort][s][p] for p in PANELS}) > 1 for s in ids)
    if varies:
        stat, p = stats.friedmanchisquare(*[col[p] for p in PANELS])
        rows.append([cohort, "Friedman", "5 vs 10 vs 20", len(ids), "", f"{stat:.4f}", f"{p:.4f}", f"{p:.4f}"])
    else:
        rows.append([cohort, "Friedman", "5 vs 10 vs 20", len(ids), "", "undefined", "undefined", "undefined"])
    raw = []
    for a, b in itertools.combinations(PANELS, 2):
        d = [y - x for x, y in zip(col[a], col[b])]
        if all(abs(v) < 1e-12 for v in d):
            s_, p_ = None, 1.0
        else:
            s_, p_ = stats.wilcoxon(col[a], col[b], zero_method="wilcox", alternative="two-sided")
        raw.append((a, b, s_, float(p_), statistics.median(d), sum(1 for v in d if abs(v) >= 1e-12)))
    order = sorted(range(len(raw)), key=lambda i: raw[i][3])
    adj, running = [None] * len(raw), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(raw) - rank) * raw[i][3]))
        adj[i] = running
    for i, (a, b, s_, p_, med, nz) in enumerate(raw):
        rows.append([cohort, "Wilcoxon signed-rank", f"{a} vs {b}", len(ids),
                     f"median {med:+.2f} pp, {nz} non-zero",
                     "" if s_ is None else f"{s_:.1f}", f"{p_:.4f}", f"{adj[i]:.4f}"])

head = ["cohort", "test", "panel_sizes", "n_common_samples", "paired_difference", "statistic", "p_value", "p_holm"]
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(head)
    w.writerows(rows)
for r in [head] + rows:
    print("  ".join(str(x) for x in r))
defined = [float(r[6]) for r in rows if r[1] == "Friedman" and r[6] != "undefined"]
holm = [float(r[7]) for r in rows if r[1] != "Friedman"]
print(f"\nSmallest Friedman p: {min(defined):.3f}. Smallest Holm-adjusted Wilcoxon p: {min(holm):.3f}.")
