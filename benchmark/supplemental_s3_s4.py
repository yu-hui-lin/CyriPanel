#!/usr/bin/env python3
"""Supplemental Tables S3 and S4, from the archived records.

S3  Aldy 4 concordance with each single reference. For every aligner, cohort and
    reference, the number of concordant test samples over all test samples, overall and
    within the three structural strata.
S4  Samples with two candidate diplotypes. For each sample whose CyriPanel output listed
    more than one candidate, the candidates, the number of CyriPanel evaluations with each
    aligner over all panel sizes, the Aldy 4 major diplotype as reported, and the
    concordant Aldy 4 runs with each aligner.

Reads aldy/aldy_results_long.csv, published_results/per_sample_v2.csv and the reference
standard in data/. The strata and the normalizer are those of caller_comparison_n5.py.
Written for release v1.0.2.

Standard library only. Run from the repository's benchmark/ directory, or pass the
benchmark directory as the first argument:
    python3 supplemental_s3_s4.py  [path/to/CyriPanel/benchmark]
"""
import collections
import csv
import os
import re
import sys

BENCH = sys.argv[1] if len(sys.argv) > 1 else "."
PER_SAMPLE = os.path.join(BENCH, "published_results", "per_sample_v2.csv")
ALDY = os.path.join(BENCH, "aldy", "aldy_results_long.csv")
GOLD = {c: os.path.join(BENCH, "data", f"CYP2D6_{c}_gold_standard_v2.csv") for c in ("cohort1", "cohort2")}
HYB = {13, 36, 68, 83}
STRATA = (("all", "All samples"), ("no_SV", "No structural variation"),
          ("CNV_only", "Copy-number variation only"), ("hybrid", "Hybrid"))


def numeric_key(value):
    return tuple(int(x) for x in re.findall(r"\d+", value))


def normalize(value):
    """Same normaliser as caller_comparison_n5.py."""
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


def stratum(gold):
    """Same strata as caller_comparison_n5.py."""
    cores = {int(x) for x in re.findall(r"\*(\d+)", gold)}
    if cores & HYB:
        return "hybrid"
    if 5 in cores or "x" in gold.lower() or "+" in gold:
        return "CNV_only"
    return "no_SV"


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    gold = {}
    for c, p in GOLD.items():
        for r in read(p):
            gold[(c, r["Sample_ID"].strip())] = r["Genotype"].strip()
    aldy = read(ALDY)

    print("=" * 88)
    print("SUPPLEMENTAL TABLE S3  Aldy 4 concordance with each single reference")
    print("=" * 88)
    tally = collections.defaultdict(lambda: [0, 0])   # (aligner, cohort, reference, stratum) -> [concordant, tested]
    for r in aldy:
        cohort = r["cohort"]
        for key in ("all", stratum(gold[(cohort, r["sample"])])):
            t = tally[(r["aligner"], cohort, r["reference"], key)]
            t[0] += r["comparison"] == "concordant"
            t[1] += 1
    print("  " + "".join(f"{h:<10}" for h in ("Cohort", "Aligner", "Reference"))
          + "".join(f"{label:>28}" for _, label in STRATA))
    n_rows = 0
    for aligner in ("BWA-MEM", "DRAGEN"):
        for cohort in ("cohort1", "cohort2"):
            refs = sorted({k[2] for k in tally if k[0] == aligner and k[1] == cohort})
            for ref in refs:
                cells = []
                for key, _ in STRATA:
                    k, n = tally.get((aligner, cohort, ref, key), (0, 0))
                    cells.append(f"{k}/{n}")
                print("  " + "".join(f"{v:<10}" for v in (cohort[-1], aligner, ref))
                      + "".join(f"{v:>28}" for v in cells))
                n_rows += 1
    print(f"\n  {n_rows} references, {len(aldy):,} genotyping runs")

    print()
    print("=" * 88)
    print("SUPPLEMENTAL TABLE S4  Samples with two candidate diplotypes")
    print("=" * 88)
    rows = read(PER_SAMPLE)
    multi = set()
    for r in rows:
        candidates = [x for x in r["called_genotype"].split(";") if normalize(x) is not None]
        if len(candidates) > 1:
            multi.add((r["cohort"].replace("_dragen", ""), r["sample_id"]))
    for cohort, sample in sorted(multi, key=lambda k: (k[0], k[1])):
        calls = collections.Counter()
        evaluations = collections.Counter()
        for r in rows:
            if r["sample_id"] == sample and r["cohort"].replace("_dragen", "") == cohort:
                aligner = "DRAGEN" if r["cohort"].endswith("_dragen") else "BWA-MEM"
                evaluations[aligner] += 1
                calls[" and ".join(r["called_genotype"].split(";"))] += 1
        majors = collections.Counter(r["raw_major"] for r in aldy if r["sample"] == sample)
        runs = {}
        for aligner in ("BWA-MEM", "DRAGEN"):
            sub = [r for r in aldy if r["sample"] == sample and r["aligner"] == aligner]
            runs[aligner] = f"{sum(r['comparison'] == 'concordant' for r in sub)}/{len(sub)}"
        print(f"  {sample}  cohort {cohort[-1]}  reference standard {gold[(cohort, sample)]}")
        print(f"    CyriPanel candidates: " + "; ".join(f"{c} ({n} evaluations)" for c, n in calls.most_common()))
        print(f"    CyriPanel evaluations, BWA-MEM and DRAGEN: {evaluations['BWA-MEM']} and {evaluations['DRAGEN']}")
        print(f"    Aldy 4 major diplotype: " + ", ".join(f"{m} ({n} runs)" for m, n in majors.most_common()))
        print(f"    Aldy 4 concordant runs, BWA-MEM and DRAGEN: {runs['BWA-MEM']} and {runs['DRAGEN']}")
    print(f"\n  {len(multi)} samples")


if __name__ == "__main__":
    main()
