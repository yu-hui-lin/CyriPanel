#!/usr/bin/env python3
"""
common_sample_analysis.py
=========================
Removes the test-set-composition confound from the reference-panel-size comparison.

THE PROBLEM
-----------
Each iteration draws n reference samples and EXCLUDES them from that iteration's test
set. A larger panel therefore removes more samples from testing, so the panel5, panel10
and panel20 arms are not evaluated on the same samples. Any difference between arms mixes
two things: the effect of reference-panel size, and the effect of which samples happened
to remain testable. That confound is the most likely explanation for concordance being
flat-to-slightly-worse at n=20 while PASS rate rises monotonically.

WHAT THIS SCRIPT DOES
---------------------
1. Restricts every arm of a cohort to the samples evaluable in ALL THREE panel sizes.
2. Recomputes arm concordance on that common set two ways:
     (a) per-iteration concordance, then mean +/- SD  -- comparable to the headline table
     (b) per-sample concordance rate, then mean over samples -- the unbiased estimate
3. Tests panel sizes against each other with PAIRED tests on per-sample rates
   (Friedman across all three; Wilcoxon signed-rank pairwise, Holm-corrected), which is
   the statistically correct comparison now that the samples are the same in every arm.

Outputs (written to OUT_DIR):
    common_samples_per_cohort.csv     which samples survived, and why others dropped
    arm_summary_common.csv            per-arm concordance restricted to common samples
    per_sample_rates_common.csv       per-sample concordance rate in each arm
    paired_tests_common.csv           Friedman + pairwise Wilcoxon results
  plus a formatted report on stdout.

Usage:
    python3 common_sample_analysis.py
    python3 common_sample_analysis.py --min-obs 5      # require >=5 evaluable iters/arm
    python3 common_sample_analysis.py --cohorts cohort1 cohort2
    python3 common_sample_analysis.py --gold-version v2

Run it on the HPC, where BENCHMARK_ROOT lives. It re-reads iteration_results.json only;
it never re-runs CyriPanel.
"""

import os
import sys
import json
import glob
import csv
import argparse
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import benchmark_config as cfg              # noqa: E402
from run_single_iteration import compare_genotype   # noqa: E402

DATA_DIR = os.path.join(cfg.CYRIPANEL_SRC_DIR, "data")

GOLD_FILES = {
    "v1": {
        "cohort1": "CYP2D6_cohort1_gold_standard.csv",
        "cohort2": "CYP2D6_cohort2_gold_standard.csv",
    },
    "v2": {
        "cohort1": "CYP2D6_cohort1_gold_standard_v2.csv",
        "cohort2": "CYP2D6_cohort2_gold_standard_v2.csv",
    },
}

# A dragen cohort shares its base cohort's gold standard.
def base_cohort(cohort):
    return cohort.replace("_dragen", "")

OUT_DIR = os.path.join(cfg.BENCHMARK_ROOT, "common_sample_analysis")


# ──────────────────────────────────────────────────────────────────────
# Loading
# ──────────────────────────────────────────────────────────────────────

def load_gold(cohort, version):
    fname = GOLD_FILES[version][base_cohort(cohort)]
    path = os.path.join(DATA_DIR, fname)
    with open(path) as f:
        return {r["Sample_ID"]: r["Genotype"] for r in csv.DictReader(f)}


def load_long_table(cohorts, panel_sizes, gold_version):
    """Return a long DataFrame: one row per (cohort, panel, iter, sample)."""
    rows = []
    for cohort in cohorts:
        gold = load_gold(cohort, gold_version)
        for panel in panel_sizes:
            base = os.path.join(cfg.BENCHMARK_ROOT, cohort, "panel%d" % panel)
            files = sorted(glob.glob(os.path.join(base, "iter_*/iteration_results.json")))
            if not files:
                print("  WARNING: no iterations found for %s panel%d (%s)"
                      % (cohort, panel, base))
                continue
            print("  %-16s panel=%-3d %3d iterations" % (cohort, panel, len(files)))
            for path in files:
                with open(path) as f:
                    data = json.load(f)
                it = data.get("summary", {}).get("iteration")
                seed = data.get("summary", {}).get("seed")
                for sid, ps in data.get("per_sample", {}).items():
                    called = ps.get("called_genotype")
                    g = gold.get(sid)
                    comp = compare_genotype(called, g)
                    rows.append({
                        "cohort": cohort,
                        "panel_size": panel,
                        "iteration": it,
                        "seed": seed,
                        "sample_id": sid,
                        "called_genotype": called,
                        "gold": g,
                        "comparison": comp,
                        "filter": ps.get("filter"),
                        # "evaluable" mirrors the denominator used in the headline table:
                        # a sample counts only when it has both a call and a gold genotype.
                        "evaluable": bool(called is not None and g is not None),
                        "concordant": comp == "concordant",
                    })
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────
# Common-sample restriction
# ──────────────────────────────────────────────────────────────────────

def find_common_samples(df, cohort, panel_sizes, min_obs):
    """
    Samples evaluable at least `min_obs` times in EVERY panel size of this cohort.
    Returns (common_set, diagnostics_rows).
    """
    sub = df[(df["cohort"] == cohort) & (df["evaluable"])]
    counts = (sub.groupby(["sample_id", "panel_size"]).size()
                 .unstack(fill_value=0)
                 .reindex(columns=panel_sizes, fill_value=0))
    ok = (counts >= min_obs).all(axis=1)
    common = set(counts.index[ok])

    diags = []
    all_ids = set(df[df["cohort"] == cohort]["sample_id"])
    for sid in sorted(all_ids):
        row = {"cohort": cohort, "sample_id": sid}
        for p in panel_sizes:
            row["n_evaluable_panel%d" % p] = int(counts.loc[sid, p]) if sid in counts.index else 0
        row["in_common_set"] = sid in common
        if sid not in common:
            zero = [p for p in panel_sizes
                    if (counts.loc[sid, p] if sid in counts.index else 0) == 0]
            if zero:
                row["reason"] = "never evaluable at panel " + ",".join(str(p) for p in zero)
            else:
                row["reason"] = "fewer than %d evaluable iterations in some arm" % min_obs
        else:
            row["reason"] = ""
        diags.append(row)
    return common, diags


# ──────────────────────────────────────────────────────────────────────
# Statistics
# ──────────────────────────────────────────────────────────────────────

def paired_tests(rates, cohort, panel_sizes):
    """
    rates: DataFrame indexed by sample_id, one column per panel size, values = concordance rate.
    Friedman across all arms; pairwise Wilcoxon signed-rank with Holm correction.
    """
    out = []
    try:
        from scipy import stats
    except ImportError:
        print("  scipy not available - skipping paired tests. "
              "Activate the cyrius_env (see benchmark_config.VENV_ACTIVATE).")
        return out

    cols = [p for p in panel_sizes if p in rates.columns]
    mat = rates[cols].dropna()
    n = len(mat)
    if n < 3 or len(cols) < 2:
        print("  Too few paired samples (n=%d) for testing." % n)
        return out

    if len(cols) >= 3:
        stat, p = stats.friedmanchisquare(*[mat[c].values for c in cols])
        out.append({"cohort": cohort, "test": "Friedman",
                    "comparison": "panel " + " vs ".join(str(c) for c in cols),
                    "n_samples": n, "statistic": stat, "p_value": p,
                    "p_adjusted": p, "method": "chi-square approximation"})

    raw = []
    for a, b in itertools.combinations(cols, 2):
        d = mat[b] - mat[a]
        if np.allclose(d, 0):
            stat, p = np.nan, 1.0
        else:
            stat, p = stats.wilcoxon(mat[a], mat[b], zero_method="wilcox",
                                     alternative="two-sided")
        raw.append((a, b, stat, p, float(np.median(d))))

    # Holm-Bonferroni over the pairwise family
    order = sorted(range(len(raw)), key=lambda i: raw[i][3])
    m = len(raw)
    adj = [None] * m
    running = 0.0
    for rank, idx in enumerate(order):
        val = min(1.0, (m - rank) * raw[idx][3])
        running = max(running, val)
        adj[idx] = running

    for i, (a, b, stat, p, med) in enumerate(raw):
        out.append({"cohort": cohort, "test": "Wilcoxon signed-rank",
                    "comparison": "panel%d vs panel%d" % (a, b),
                    "n_samples": n, "statistic": stat, "p_value": p,
                    "p_adjusted": adj[i],
                    "method": "median difference (panel%d - panel%d) = %+.4f" % (b, a, med)})
    return out


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohorts", nargs="+",
                    default=["cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen"])
    ap.add_argument("--panel-sizes", nargs="+", type=int,
                    default=cfg.REFERENCE_PANEL_SIZES)
    ap.add_argument("--gold-version", choices=["v1", "v2"], default="v2")
    ap.add_argument("--min-obs", type=int, default=1,
                    help="minimum evaluable iterations per arm for a sample to be common "
                         "(default 1; use 5 for a stricter set)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading iteration results (gold standard %s)..." % args.gold_version)
    df = load_long_table(args.cohorts, args.panel_sizes, args.gold_version)
    if df.empty:
        print("No iteration results found. Check BENCHMARK_ROOT in benchmark_config.py.")
        return 1
    print("Loaded %d sample-observations.\n" % len(df))

    all_diags, all_arms, all_rates, all_tests = [], [], [], []

    for cohort in args.cohorts:
        if df[df["cohort"] == cohort].empty:
            continue
        common, diags = find_common_samples(df, cohort, args.panel_sizes, args.min_obs)
        all_diags.extend(diags)

        n_total = df[df["cohort"] == cohort]["sample_id"].nunique()
        print("=" * 100)
        print("  %s - %d of %d samples evaluable in all %d arms (min %d iterations each)"
              % (cohort.upper(), len(common), n_total, len(args.panel_sizes), args.min_obs))
        print("=" * 100)
        if not common:
            print("  No common samples; skipping.\n")
            continue

        sub = df[(df["cohort"] == cohort) & (df["sample_id"].isin(common)) & df["evaluable"]]

        # (a) per-iteration concordance on the common set, then mean +/- SD
        per_iter = (sub.groupby(["panel_size", "iteration"])["concordant"]
                       .mean().reset_index(name="concordance"))
        # (b) per-sample concordance rate, then mean over samples
        per_samp = (sub.groupby(["panel_size", "sample_id"])["concordant"]
                       .mean().reset_index(name="rate"))
        rates = per_samp.pivot(index="sample_id", columns="panel_size", values="rate")

        print()
        print("  %-8s | %-24s | %-24s | %-6s" % (
            "Panel", "per-iteration mean+/-SD", "per-sample mean+/-SD", "#iter"))
        print("  " + "-" * 78)
        for p in args.panel_sizes:
            pi = per_iter[per_iter["panel_size"] == p]["concordance"]
            psx = per_samp[per_samp["panel_size"] == p]["rate"]
            if pi.empty:
                continue
            print("  n=%-6d | %.4f +/- %.4f        | %.4f +/- %.4f        | %d"
                  % (p, pi.mean(), pi.std(), psx.mean(), psx.std(), len(pi)))
            all_arms.append({
                "cohort": cohort, "panel_size": p,
                "n_common_samples": len(common),
                "n_iterations": len(pi),
                "concordance_per_iteration_mean": pi.mean(),
                "concordance_per_iteration_sd": pi.std(),
                "concordance_per_sample_mean": psx.mean(),
                "concordance_per_sample_sd": psx.std(),
                "gold_version": args.gold_version,
                "min_obs": args.min_obs,
            })

        r = rates.reset_index()
        r.insert(0, "cohort", cohort)
        all_rates.append(r)

        print()
        tests = paired_tests(rates, cohort, args.panel_sizes)
        if tests:
            print("  Paired tests on per-sample concordance rates (same samples in every arm):")
            for t in tests:
                print("    %-22s %-24s n=%-4d p=%.4g  p_adj=%.4g  %s"
                      % (t["test"], t["comparison"], t["n_samples"],
                         t["p_value"], t["p_adjusted"], t["method"]))
        all_tests.extend(tests)

        # Samples whose rate moves most across panel sizes - candidates for the case study
        if rates.shape[1] >= 2:
            spread = (rates.max(axis=1) - rates.min(axis=1)).sort_values(ascending=False)
            movers = spread[spread > 0].head(8)
            if len(movers):
                print()
                print("  Largest movers across panel sizes (candidate case studies):")
                for sid, sp in movers.items():
                    cells = "  ".join("n=%d:%.2f" % (p, rates.loc[sid, p])
                                      for p in args.panel_sizes if p in rates.columns
                                      and pd.notna(rates.loc[sid, p]))
                    print("    %-12s spread=%.2f   %s" % (sid, sp, cells))
        print()

    # ── Write outputs ────────────────────────────────────────────────
    pd.DataFrame(all_diags).to_csv(
        os.path.join(args.out_dir, "common_samples_per_cohort.csv"), index=False)
    pd.DataFrame(all_arms).to_csv(
        os.path.join(args.out_dir, "arm_summary_common.csv"), index=False)
    if all_rates:
        pd.concat(all_rates, ignore_index=True).to_csv(
            os.path.join(args.out_dir, "per_sample_rates_common.csv"), index=False)
    pd.DataFrame(all_tests).to_csv(
        os.path.join(args.out_dir, "paired_tests_common.csv"), index=False)

    print("=" * 100)
    print("Wrote 4 CSVs to %s" % args.out_dir)
    print("""
Reading the output
------------------
* "per-iteration mean+/-SD" is directly comparable to the headline table, but now every
  arm is scored on the same samples, so a difference between arms is attributable to the
  reference panel rather than to which samples remained testable.
* "per-sample mean+/-SD" weights every sample equally rather than every iteration, and is
  the estimate to quote if a reviewer asks for accuracy per sample.
* The SD in either column is dispersion across iterations or samples - it is NOT a
  confidence interval on accuracy. Say so in the Methods.
* The paired tests are the defensible way to claim "panel size does not matter": the same
  samples appear in every arm, so the pairing is real.
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
