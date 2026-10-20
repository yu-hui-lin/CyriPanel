#!/usr/bin/env python3
"""Reconstruct exploratory panel-size tests from archived evaluation records.

Requires reproduce_tables.py, NumPy 1.26.4 and SciPy 1.15.2.
No caller execution. Outputs go to a new directory.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import importlib.util
import io
import itertools
import json
from pathlib import Path
import statistics as st
import sys

import numpy as np
import scipy
from scipy import stats

HELPER_SHA = "8f51f22ad5b397bb62c92424784dd4ea8e7a86e5a31cae4f5d61261fc62a0d04"
EXPECTED = "8213ba6939d719776ca7279d7548bf0c82e500dd51cc3553bcad687f9d4fbaba"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fmt(value):
    return "NA" if value is None else f"{value:.6f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    require(not args.outdir.exists() and not args.outdir.is_symlink(), "Output directory already exists.")
    require(np.__version__ == "1.26.4" and scipy.__version__ == "1.15.2",
            "Activate cyripanel-rebuild: expected NumPy 1.26.4 and SciPy 1.15.2.")
    bench = Path(__file__).resolve().parent
    helper = bench / "reproduce_tables.py"
    require(sha(helper.read_bytes()) == HELPER_SHA, "Unexpected reproduce_tables.py version.")
    spec = importlib.util.spec_from_file_location("table_reproduction", helper)
    module = importlib.util.module_from_spec(spec)
    sys.dont_write_bytecode = True
    spec.loader.exec_module(module)
    names = ["published_results/per_sample_v2.csv"] + [
        f"data/CYP2D6_{c}_gold_standard_v2.csv" for c in module.ARMS[:2]]
    data = {}
    for name in names:
        raw = (bench / name).read_bytes()
        count, digest = module.INPUTS[name]
        require(sha(raw) == digest, "Input checksum mismatch: " + name)
        data[name] = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        require(len(data[name]) == count, "Unexpected input row count.")
    gold = {c: {r["Sample_ID"]: r["Genotype"] for r in
                data[f"data/CYP2D6_{c}_gold_standard_v2.csv"]} for c in module.ARMS[:2]}
    outcomes, keys, iterations = defaultdict(list), set(), Counter()
    for row in data[names[0]]:
        arm, sample, n, i = row["cohort"], row["sample_id"], int(row["panel_size"]), int(row["iteration"])
        key = (arm, sample, n, i)
        require(key not in keys, "Duplicate evaluation.")
        keys.add(key)
        iterations[arm, n, i] += 1
        target = gold[arm.replace("_dragen", "")][sample]
        require(row["gold_v2"] == target, "Reference-standard mismatch.")
        candidates = [module.normalize(x) for x in row["called_genotype"].split(";")]
        outcomes[arm, sample, n].append(module.normalize(target) in candidates)
    wanted = {(a, n, i): len(gold[a.replace("_dragen", "")]) - n
              for a in module.ARMS for n in module.SIZES for i in range(30)}
    require(dict(iterations) == wanted, "Incomplete iteration records.")
    rates = {key: Fraction(sum(values), len(values)) for key, values in outcomes.items()}
    rate_rows, member_rows, test_rows, signatures, summary = [], [], [], [], []
    for arm in module.ARMS:
        for sample in sorted(gold[arm.replace("_dragen", "")]):
            counts = [len(outcomes[arm, sample, n]) for n in module.SIZES]
            member_rows.append([arm, sample, *counts, min(counts) >= 1, min(counts) >= 5])
            for n in module.SIZES:
                values = outcomes[arm, sample, n]
                rate_rows.append([arm, sample, n, len(values), sum(values), float(rates[arm, sample, n])])
    for minimum in (1, 5):
        for arm in module.ARMS:
            all_samples = sorted(gold[arm.replace("_dragen", "")])
            common = [s for s in all_samples if all(len(outcomes[arm, s, n]) >= minimum for n in module.SIZES)]
            require(bool(common), "No common samples.")
            excluded = sorted(set(all_samples) - set(common))
            matrix = [[float(rates[arm, s, n]) for s in common] for n in module.SIZES]
            varying = sum(len({rates[arm, s, n] for n in module.SIZES}) > 1 for s in common)
            if varying:
                result = stats.friedmanchisquare(*matrix)
                fstat, fp = float(result.statistic), float(result.pvalue)
            else:
                fstat, fp = None, None
            test_rows.append([minimum, arm, "Friedman", "5/10/20", len(common), varying,
                              fstat, fp, None, None, "chi-square approximation" if varying else "undefined: all tied"])
            pairs = []
            for a, b in itertools.combinations(module.SIZES, 2):
                # Subtract exact ratios before converting to float, preserving true ties.
                difference = [float(rates[arm, s, b] - rates[arm, s, a]) for s in common]
                nonzero = sum(d != 0 for d in difference)
                if nonzero:
                    result = stats.wilcoxon(difference, zero_method="wilcox", correction=False,
                                           alternative="two-sided", method="asymptotic")
                    w, p = float(result.statistic), float(result.pvalue)
                else:
                    w, p = None, 1.0
                pairs.append([a, b, nonzero, w, p, 100 * st.median(difference)])
            adjusted, previous = [0.0] * 3, 0.0
            for rank, index in enumerate(sorted(range(3), key=lambda j: pairs[j][4])):
                previous = max(previous, min(1.0, (3 - rank) * pairs[index][4]))
                adjusted[index] = previous
            signature = [str(minimum), arm, str(len(common)), ",".join(excluded), fmt(fstat), fmt(fp)]
            for (a, b, nonzero, w, p, median), adj in zip(pairs, adjusted):
                status = "asymptotic; wilcox; no continuity correction" if nonzero else "all zero: p=1 convention"
                test_rows.append([minimum, arm, "Wilcoxon", f"{a}/{b}", len(common),
                                  nonzero, w, p, adj, median, status])
                signature.extend([str(nonzero), fmt(w), fmt(p), fmt(adj), fmt(median)])
            signatures.append(signature)
            summary.append([str(minimum), arm, str(len(common)), fmt(fp),
                            *[fmt(p) for p in adjusted], ",".join(excluded) or "none"])
    require(sha(json.dumps(signatures, separators=(",", ":")).encode()) == EXPECTED,
            "Results differ from the independently verified reconstruction.")
    files = {}
    def add_csv(name, header, rows):
        stream = io.StringIO(newline="")
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
        files[name] = stream.getvalue()
    add_csv("per_sample_rates.csv", ["arm", "sample", "panel_size", "evaluations", "concordant", "rate"], rate_rows)
    add_csv("sample_membership.csv", ["arm", "sample", "n5_evaluations", "n10_evaluations", "n20_evaluations",
                                    "included_min1", "included_min5"], member_rows)
    add_csv("paired_tests.csv", ["min_evaluations", "arm", "test", "panels", "samples", "informative_samples",
                               "statistic", "p_raw", "p_holm", "median_difference_pp", "method"], test_rows)
    header = ["Minimum evaluations", "Arm", "Samples", "Friedman p",
              "Holm p: 5/10", "Holm p: 5/20", "Holm p: 10/20", "Excluded"]
    files["summary.md"] = "\n".join(["# Exploratory panel-size reconstruction", "",
        "| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |",
        *["| " + " | ".join(row) + " |" for row in summary], "",
        "All 24 pairwise median differences are zero (12 per minimum-evaluation setting).",
        "NA: Friedman is undefined when every sample has identical rates across panel sizes.",
        "These are exploratory tests; shared reference samples create dependence.",
        "Nonsignificance does not establish equivalence or non-inferiority.", ""])
    files["INPUT_SHA256SUMS"] = "".join(module.INPUTS[name][1] + "  " + name + "\n" for name in names)
    files["validation.txt"] = (
        f"NumPy: {np.__version__}\nSciPy: {scipy.__version__}\nPython: {sys.version.split()[0]}\n"
        f"Analysis script SHA256: {sha(Path(__file__).read_bytes())}\nScoring helper SHA256: {HELPER_SHA}\n"
        "Scope: new reconstruction from the archived CSV, not proof of the original HPC execution.\n"
        "Candidate-inclusive scoring; no-calls retained; common samples across all three panel sizes.\n"
        "Exact Fraction differences converted to float preserve true ties for Wilcoxon.\n"
        "Wilcoxon uses asymptotic p-values, wilcox zeros, two-sided, no continuity correction.\n"
        "Holm correction: three panel pairs per arm and minimum-evaluation setting.\n"
        "All-zero Wilcoxon: statistic undefined, p=1 by convention; all-tied Friedman: undefined.\n"
        "Floating-rate subtraction previously broke a true tie for cohort2 BWA 5/10.\n"
        "Correcting that tie changes its raw p from 0.345231 to 0.278517,\n"
        "and Holm p from 1.000000 to 0.835551, for both minimum-evaluation settings.\n"
        "The manuscript's rounded minimum p-value summaries remain consistent.\n"
        "With minimum 5 evaluations, HG01530 (24/17/4 evaluations) is excluded in both cohort2 arms.\n"
        "All 8 arm/settings and 32 test records match the independent reconstruction.\n"
    )
    files["SHA256SUMS"] = "".join(sha(text.encode()) + "  " + name + "\n" for name, text in sorted(files.items()))
    args.outdir.mkdir(parents=True, exist_ok=False)
    for name, text in files.items():
        with (args.outdir / name).open("x", encoding="utf-8", newline="") as handle:
            handle.write(text)
    for row in summary:
        print(" ".join(row[:4]) + " Holm=" + ",".join(row[4:7]) + " excluded=" + row[7])
    print("PANEL_SIZE_TESTS_REPRODUCED:", args.outdir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, ImportError) as error:
        sys.exit("STOP: " + str(error))
