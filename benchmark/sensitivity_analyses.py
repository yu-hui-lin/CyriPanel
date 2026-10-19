#!/usr/bin/env python3
"""Two sensitivity analyses from the archived records. Nothing is re-run.

Run recalculate.py first, because part 1 reads recalculate_output/common_sample_rates.csv.
Needs numpy, because it imports the scoring functions from run_single_iteration.py.

1. The g.42129180A>T exclusion.
   This variant was removed from the D7-contamination filter after review of benchmark
   calls (caller/call_variants.py, D7_CONTAMINATION_PRONE_VARS). The filter acts only in
   configurations with *68 or an exon 9 hybrid, so the exclusion can change calls only in
   samples that carry *49 together with a hybrid allele. The script finds those samples in
   the reference standard and recomputes mean per-sample concordance with five references
   as if every evaluation of the samples it rescued had been discordant. That is the
   largest effect the exclusion could have had on these estimates. It also counts *49
   calls in samples without *49, which the exclusion could have produced.

2. Discordance types.
   Every discordant evaluation is classified by comparing copy structure (gene copies per
   haplotype, hybrid alleles and deletions) with the reference standard. Differences in
   allele definitions between allele-resource versions can explain a discordant call only
   when the structure matches and the allele names differ, or when no call is made.
"""
import csv
import collections
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_single_iteration import normalize_genotype  # noqa: E402

HYBRID = {"*36", "*68", "*13", "*83"}
RATES = os.path.join(HERE, "recalculate_output", "common_sample_rates.csv")
if not os.path.exists(RATES):
    sys.exit(f"{RATES} not found. Run recalculate.py first.")

gold = {}
for c in ("cohort1", "cohort2"):
    for r in csv.DictReader(open(os.path.join(HERE, "data", f"CYP2D6_{c}_gold_standard_v2.csv"))):
        gold[r["Sample_ID"].strip()] = r["Genotype"].strip()


def alleles(gt):
    return [a for h in gt.split("/") for a in h.split("+")]


def core(a):
    m = re.match(r"(\*\d+)", a.strip())
    return m.group(1) if m else a.strip()


hybrid_samples = {s for s, g in gold.items() if HYBRID & {core(a) for a in alleles(g)}}

# ---------------------------------------------------------------- 1. g.42129180A>T exclusion
print("=" * 80)
print("1. g.42129180A>T exclusion")
print("=" * 80)
carriers = sorted(s for s, g in gold.items() if "*49" in {core(a) for a in alleles(g)})
affected = [s for s in carriers if s in hybrid_samples]
print("*49 carriers:", ", ".join(f"{s} {gold[s]}" for s in carriers))
print("*49 with a hybrid allele, where the filter applies:", ", ".join(affected))

rates = collections.defaultdict(dict)
for r in csv.DictReader(open(RATES)):
    if r["panel_size"] == "5":
        rates[r["cohort"]][r["sample_id"]] = float(r["concordance_pct"])

for aligner, arms in (("BWA-MEM", ("cohort1", "cohort2")), ("DRAGEN", ("cohort1_dragen", "cohort2_dragen"))):
    per = {s: v for a in arms for s, v in rates[a].items()}
    rescued = [s for s in affected if per.get(s, 0) > 0]
    bound = {s: (0.0 if s in rescued else v) for s, v in per.items()}
    for label, keep in (("all samples", lambda s: True), ("hybrid stratum", lambda s: s in hybrid_samples)):
        ids = [s for s in per if keep(s)]
        print(f"  {aligner:8} {label:15} n={len(ids):3}  archived {statistics.mean(per[s] for s in ids):6.2f}%"
              f"  with {', '.join(rescued)} discordant {statistics.mean(bound[s] for s in ids):6.2f}%")

false49 = collections.Counter()
for r in csv.DictReader(open(os.path.join(HERE, "published_results", "per_sample_v2.csv"))):
    if "*49" in (r["called_genotype"] or "") and "*49" not in r["gold"]:
        false49[r["sample_id"]] += 1
print("CyriPanel evaluations calling *49 in a sample without *49:", sum(false49.values()))

# ---------------------------------------------------------------- 2. discordance types
print()
print("=" * 80)
print("2. Discordance types")
print("=" * 80)


def structure(gt):
    n = normalize_genotype(gt)
    if not n:
        return None
    haps = []
    for h in n.split("/"):
        cls = []
        for a in h.split("+"):
            a = core(a)
            cls.append(a if a in HYBRID else ("deletion" if a == "*5" else "CYP2D6"))
        haps.append(tuple(sorted(cls)))
    return tuple(sorted(haps))


def kind(called, ref):
    if not called or called in ("None", ""):
        return "no call"
    if any(structure(c) == structure(ref) for c in called.split(";")):
        return "same structure, other alleles"
    return "different copy structure"


counts, samples = collections.Counter(), collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(HERE, "published_results", "per_sample_v2.csv"))):
    cands = r["called_genotype"].split(";") if r["called_genotype"] else []
    if any(normalize_genotype(c) == normalize_genotype(r["gold"]) for c in cands):
        continue
    aligner = "DRAGEN" if r["cohort"].endswith("_dragen") else "BWA-MEM"
    k = kind(r["called_genotype"], r["gold"])
    counts[("CyriPanel", aligner, k)] += 1
    samples[("CyriPanel", k)].add(r["sample_id"])
for r in csv.DictReader(open(os.path.join(HERE, "aldy", "aldy_results_long.csv"))):
    if r["comparison"] == "concordant":
        continue
    k = kind(r["normalised"], r["gold"])
    counts[("Aldy 4", r["aligner"], k)] += 1
    samples[("Aldy 4", k)].add(r["sample"])
for key in sorted(counts):
    print(f"  {key[0]:9} {key[1]:8} {key[2]:30} {counts[key]:5} evaluations")
for key in sorted(samples):
    if key[1] != "different copy structure":
        print(f"  {key[0]:9} {key[1]:30} samples: {', '.join(sorted(samples[key]))}")
