#!/usr/bin/env python3
"""Table 5, CyriPanel per-iteration results by reference-panel size, from the archived records.

Reads published_results/per_sample_v2.csv, which has one row per CyriPanel evaluation,
and prints for each cohort and reference-panel size:
  * candidate-inclusive concordance per iteration, as mean and SD over the 30 iterations,
    for each aligner. A call is concordant when any reported candidate equals the
    reference standard after normalization. No-calls stay in the denominator and count
    as discordant.
  * the mean number of discordant samples per iteration, for each aligner.
  * genotype-filter outcomes as percentages of all scheduled evaluations. They are
    identical between aligners, which the script checks before printing them once.

The normalizer is the one in recalculate.py and caller_comparison_n5.py. SD is the
sample standard deviation.

recalculate.py, the script used for the preprint, computes the concordance columns. This
script was written for release v1.0.2 to print every column of Table 5 in one place,
including the discordant counts and the genotype-filter outcomes. It reproduces every
value in Table 5.

Standard library only. Run from the repository's benchmark/ directory, or pass the
benchmark directory as the first argument:
    python3 table5.py  [path/to/CyriPanel/benchmark]
"""
import collections
import csv
import os
import re
import statistics
import sys

BENCH = sys.argv[1] if len(sys.argv) > 1 else "."
PER_SAMPLE = os.path.join(BENCH, "published_results", "per_sample_v2.csv")
COHORTS = ("cohort1", "cohort2")
PANELS = (5, 10, 20)
ITERATIONS = 30
ARMS = (("BWA-MEM", ""), ("DRAGEN", "_dragen"))
FILTERS = (("PASS", "PASS"),
           ("Ambiguous", "More_than_one_possible_genotype"),
           ("Not assigned", "Not_assigned_to_haplotypes"),
           ("No-call", ""))
FAILS = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(label)


def numeric_key(value):
    return tuple(int(x) for x in re.findall(r"\d+", value))


def normalize(value):
    """Same normaliser as recalculate.py and caller_comparison_n5.py."""
    if value is None or value.strip() in ("", "None", "Error", "NA", "nan", "NaN"):
        return None
    haplotypes = []
    for hap in value.strip().split("/"):
        copies = []
        for token in hap.split("+"):
            token = token.strip()
            mult = re.fullmatch(r"(\*\d+(?:\.\d+)?)x(\d+)", token)
            copies.extend([mult[1]] * int(mult[2]) if mult else [token])
        haplotypes.append("+".join(sorted(copies, key=numeric_key)))
    return "/".join(sorted(haplotypes, key=numeric_key))


def main():
    with open(PER_SAMPLE, newline="") as fh:
        rows = list(csv.DictReader(fh))

    tally = collections.defaultdict(lambda: [0, 0])   # (arm, panel, iteration) -> [concordant, scheduled]
    outcome = {}                                       # (aligner, cohort, panel, iteration, sample) -> filter
    for r in rows:
        arm, panel, it = r["cohort"], int(r["panel_size"]), int(r["iteration"])
        target = normalize(r["gold"])
        candidates = [c for c in (normalize(x) for x in r["called_genotype"].split(";")) if c is not None]
        tally[(arm, panel, it)][0] += target is not None and target in candidates
        tally[(arm, panel, it)][1] += 1
        aligner = "DRAGEN" if arm.endswith("_dragen") else "BWA-MEM"
        outcome[(aligner, arm.replace("_dragen", ""), panel, it, r["sample_id"])] = r["filter"]

    print("=" * 88)
    print("INPUT CHECKS")
    print("=" * 88)
    check("per_sample_v2.csv has 16,680 evaluations", len(rows) == 16680, str(len(rows)))
    complete = all((c + suffix, n, i) in tally
                   for c in COHORTS for _, suffix in ARMS for n in PANELS for i in range(ITERATIONS))
    check("30 iterations for every cohort, aligner and panel size", complete)
    sizes = {}
    for c in COHORTS:
        for n in PANELS:
            sizes[(c, n)] = {tally[(c + suffix, n, i)][1] for _, suffix in ARMS for i in range(ITERATIONS)}
    check("test-set size is constant within each cohort and panel size",
          all(len(v) == 1 for v in sizes.values()))
    pairs = collections.Counter()
    differ = 0
    for (aligner, c, n, i, s), f in outcome.items():
        if aligner != "BWA-MEM":
            continue
        other = outcome.get(("DRAGEN", c, n, i, s))
        pairs[c] += 1
        differ += other != f
    check("filter outcomes identical between aligners in every paired evaluation", differ == 0,
          f"{sum(pairs.values()):,} pairs, {pairs['cohort1']:,} in cohort 1 and "
          f"{pairs['cohort2']:,} in cohort 2, {differ} differ")
    known = {code for _, code in FILTERS}
    other_filters = sorted({f for f in outcome.values() if f not in known})
    check("no filter values other than PASS, ambiguous, not assigned and no-call",
          not other_filters, ", ".join(other_filters))

    print()
    print("=" * 88)
    print("TABLE 5  CyriPanel per-iteration results by reference-panel size")
    print("=" * 88)
    head = ["Cohort", "Refs", "Tests", "Conc BWA-MEM", "Conc DRAGEN", "Disc BWA", "Disc DRG"]
    head += [label for label, _ in FILTERS]
    print("  " + "".join(f"{h:>14}" if k > 2 else f"{h:>7}" for k, h in enumerate(head)))
    for c in COHORTS:
        for n in PANELS:
            line = [c[-1], str(n), str(next(iter(sizes[(c, n)])))]
            disc = []
            for _, suffix in ARMS:
                values = [tally[(c + suffix, n, i)] for i in range(ITERATIONS)]
                rates = [100.0 * k / t for k, t in values]
                line.append(f"{statistics.mean(rates):.1f} ± {statistics.stdev(rates):.1f}")
                disc.append(f"{statistics.mean(t - k for k, t in values):.2f}")
            line += disc
            counts = collections.Counter(f for (a, cc, nn, _, _), f in outcome.items()
                                         if a == "BWA-MEM" and cc == c and nn == n)
            total = sum(counts.values())
            line += [f"{100.0 * counts[code] / total:.1f}" for _, code in FILTERS]
            print("  " + "".join(f"{v:>14}" if k > 2 else f"{v:>7}" for k, v in enumerate(line)))
    print("\n  Concordance is candidate-inclusive, mean ± SD over 30 iterations (%).")
    print("  Disc is the mean number of discordant samples per iteration.")
    print("  Filter outcomes are percentages of all scheduled evaluations.")

    print()
    print("=" * 88)
    print("FAILED CHECKS: " + (", ".join(FAILS) if FAILS else "none"))
    print("=" * 88)


if __name__ == "__main__":
    main()
