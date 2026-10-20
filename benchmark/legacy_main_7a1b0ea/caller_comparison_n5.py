#!/usr/bin/env python3
"""Recalculate every number in the Results section "Diplotype concordance and
comparison with Aldy 4" under the Methods in "Concordance and statistical analysis".

Rules implemented (Methods, verbatim intent):
  * CyriPanel: evaluations with a five-sample reference panel only (panel_size == 5).
  * Aldy 4: every eligible single-reference run (aldy_results_long.csv).
  * Per-sample concordance = concordant evaluations / scheduled evaluations of that
    sample. No-calls stay in the denominator and count as discordant.
  * Candidate-inclusive (primary): any reported candidate equals the reference standard.
  * Single-candidate (sensitivity): exactly one candidate, equal to the reference standard.
  * Mean concordance = unweighted mean of per-sample rates over samples evaluated by
    both callers (each sample equal weight, cohorts pooled).
  * Paired difference = mean within-sample difference, CyriPanel minus Aldy 4, in pp.
  * Strata from the reference-standard diplotype: hybrid = any *13/*36/*68/*83;
    CNV_only = deletion (*5) or duplication, no hybrid; no_SV = everything else.

Standard library only. Run from the repository's benchmark/ directory, or pass the
benchmark directory as the first argument:
    python3 caller_comparison_n5.py  [path/to/CyriPanel/benchmark]

Release v1.0.2 renamed two columns of per_sample_v2.csv, gold_v2 to gold and
comparison_v2 to comparison_single. This copy reads the new names. Nothing else changed.
"""
import collections
import csv
import os
import re
import statistics
import sys

BENCH = sys.argv[1] if len(sys.argv) > 1 else "."
PER_SAMPLE = os.path.join(BENCH, "published_results", "per_sample_v2.csv")
ALDY = os.path.join(BENCH, "aldy", "aldy_results_long.csv")
GOLD = {c: os.path.join(BENCH, "data", f"CYP2D6_{c}_gold_standard_v2.csv") for c in ("cohort1", "cohort2")}
DIPLOIDS = {"cohort1": os.path.join(BENCH, "data", "CYP2D6_cohort1_33diploids.csv"),
            "cohort2": os.path.join(BENCH, "data", "CYP2D6_cohort2_28diploids.csv")}
PANEL = 5
ALIGNERS = {"BWA-MEM": ("cohort1", "cohort2"), "DRAGEN": ("cohort1_dragen", "cohort2_dragen")}
HYB = {13, 36, 68, 83}
FAILS = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(label)


def numeric_key(value):
    return tuple(int(x) for x in re.findall(r"\d+", value))


def normalize(value):
    """Same normaliser as recalculate.py."""
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
    cores = {int(x) for x in re.findall(r"\*(\d+)", gold)}
    if cores & HYB:
        return "hybrid"
    if 5 in cores or "x" in gold.lower() or "+" in gold:
        return "CNV_only"
    return "no_SV"


def multicopy_hybrid_tandem(gold):
    """True if one haplotype carries a hybrid allele together with at least one more copy."""
    for hap in normalize(gold).split("/"):
        parts = hap.split("+")
        if len(parts) >= 2 and any(int(re.findall(r"\d+", p)[0]) in HYB for p in parts):
            return True
    return False


def pct(x):
    return 100.0 * x


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    # ------------------------------------------------------------------ inputs
    gold = {c: {r["Sample_ID"].strip(): r["Genotype"].strip() for r in read(p)} for c, p in GOLD.items()}
    diploids = {}
    for c, p in DIPLOIDS.items():
        rows = read(p)
        col = next(k for k in rows[0] if "sample" in k.lower() or "id" in k.lower())
        diploids[c] = {r[col].strip() for r in rows}
    strata = {(c, s): stratum(g) for c in gold for s, g in gold[c].items()}

    print("=" * 88)
    print("INPUT CHECKS")
    print("=" * 88)
    n_str = collections.Counter(strata.values())
    check("116 samples, strata 61 / 27 / 28",
          (len(strata), n_str["no_SV"], n_str["CNV_only"], n_str["hybrid"]) == (116, 61, 27, 28),
          f"n={len(strata)} no_SV={n_str['no_SV']} CNV_only={n_str['CNV_only']} hybrid={n_str['hybrid']}")
    check("33 and 28 diploid references", (len(diploids["cohort1"]), len(diploids["cohort2"])) == (33, 28))

    ps = read(PER_SAMPLE)
    check("per_sample_v2.csv has 16,680 evaluations", len(ps) == 16680, str(len(ps)))
    cyri = collections.defaultdict(lambda: dict(n=0, prim=0, strict=0, nocall=0, amb=0, amb_gold_in=0))
    mism = 0
    for r in ps:
        if int(r["panel_size"]) != PANEL:
            continue
        arm, s = r["cohort"], r["sample_id"]
        base = arm.replace("_dragen", "")
        target = normalize(gold[base][s])
        assert target is not None and r["gold"] == gold[base][s]
        cands = [normalize(x) for x in r["called_genotype"].split(";") if normalize(x) is not None]
        strict = len(cands) == 1 and cands[0] == target
        prim = target in cands
        archived = r["comparison_single"]
        expect = "no_call" if not cands else ("concordant" if normalize(r["called_genotype"]) == target else "discordant")
        mism += expect != archived
        d = cyri[(arm, s)]
        d["n"] += 1
        d["prim"] += prim
        d["strict"] += strict
        d["nocall"] += not cands
        if len(cands) > 1:
            d["amb"] += 1
            d["amb_gold_in"] += prim
            d.setdefault("amb_ncands", set()).add(len(cands))
    check("strict labels reproduce archived comparison_single at panel size 5", mism == 0, f"{mism} mismatches")

    al = read(ALDY)
    check("aldy_results_long.csv has 6,814 runs", len(al) == 6814, str(len(al)))
    check("Aldy reported one diplotype in every run", all(r["n_solutions"] == "1" for r in al))
    aldy = collections.defaultdict(lambda: dict(n=0, conc=0))
    per_ref = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    bad_self = 0
    for r in al:
        arm, s, ref = r["arm"], r["sample"], r["reference"]
        base = arm.replace("_dragen", "")
        assert ref in diploids[base]
        bad_self += ref == s
        ok = r["comparison"] == "concordant"
        aldy[(arm, s)]["n"] += 1
        aldy[(arm, s)]["conc"] += ok
        for key in ("all", strata[(base, s)]):
            per_ref[(arm, key)][ref][0] += ok
            per_ref[(arm, key)][ref][1] += 1
    check("no Aldy run used the test sample as its own reference", bad_self == 0)
    ok_n = all(aldy[(a, s)]["n"] == len(diploids[a.replace("_dragen", "")]) - (s in diploids[a.replace("_dragen", "")])
               for a in ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen")
               for s in gold[a.replace("_dragen", "")])
    check("every sample has one Aldy run per eligible reference (33/32 and 28/27)", ok_n)
    both = all((a, s) in cyri and cyri[(a, s)]["n"] > 0 and aldy[(a, s)]["n"] > 0
               for a in ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen")
               for s in gold[a.replace("_dragen", "")])
    check("all 116 samples evaluated by both callers in every arm", both)

    def rate(d, k):
        return d[k] / d["n"]

    def samples(aligner, strat="all", cohort=None):
        out = []
        for arm in ALIGNERS[aligner]:
            base = arm.replace("_dragen", "")
            if cohort and base != cohort:
                continue
            for s in sorted(gold[base]):
                if strat == "all" or strata[(base, s)] == strat:
                    out.append((arm, s))
        return out

    def compare(aligner, strat="all", key="prim", cohort=None):
        ss = samples(aligner, strat, cohort)
        c = [rate(cyri[x], key) for x in ss]
        a = [rate(aldy[x], "conc") for x in ss]
        diff = [ci - ai for ci, ai in zip(c, a)]
        return len(ss), pct(statistics.mean(c)), pct(statistics.mean(a)), pct(statistics.mean(diff))

    # -------------------------------------------------------- main comparison
    print()
    print("=" * 88)
    print("MEAN PER-SAMPLE CONCORDANCE  (CyriPanel n=5 vs Aldy 4 all single references)")
    print("=" * 88)
    res = {}
    for key, lbl in (("prim", "candidate-inclusive (primary)"), ("strict", "single-candidate (sensitivity)")):
        print(f"\n  {lbl}")
        print(f"  {'aligner':9}{'stratum':10}{'cohort':8}{'n':>5}{'CyriPanel':>11}{'Aldy 4':>9}{'diff (pp)':>11}")
        for aligner in ALIGNERS:
            for strat in ("all", "no_SV", "CNV_only", "hybrid"):
                for cohort in (None, "cohort1", "cohort2"):
                    if cohort and strat != "all":
                        continue
                    n, c, a, d = compare(aligner, strat, key, cohort)
                    res[(key, aligner, strat, cohort)] = (n, c, a, d)
                    print(f"  {aligner:9}{strat:10}{(cohort or 'pooled'):8}{n:>5}{c:>10.2f}%{a:>8.2f}%{d:>+11.2f}")

    # ------------------------------------------------------- hybrid samples
    print()
    print("=" * 88)
    print("HYBRID SAMPLES, per-sample concordance  (Cyri = CyriPanel n=5, candidate-inclusive)")
    print("=" * 88)
    print(f"  {'cohort':8}{'sample':9}{'gold':26}{'tandem':>7}{'Aldy BWA':>10}{'Cyri BWA':>10}{'Aldy DRG':>10}{'Cyri DRG':>10}")
    hyb = [(c, s) for (c, s), t in sorted(strata.items()) if t == "hybrid"]
    for c, s in hyb:
        ab, cb = aldy[(c, s)], cyri[(c, s)]
        ad, cd = aldy[(c + "_dragen", s)], cyri[(c + "_dragen", s)]
        print(f"  {c:8}{s:9}{gold[c][s]:26}{'yes' if multicopy_hybrid_tandem(gold[c][s]) else 'no':>7}"
              f"{ab['conc']:>4}/{ab['n']:<2} {cb['prim']:>3}/{cb['n']:<2}  {ad['conc']:>3}/{ad['n']:<2}  {cd['prim']:>3}/{cd['n']:<2}")
    aldy_zero_bwa = [(c, s) for c, s in hyb if aldy[(c, s)]["conc"] == 0]
    print(f"\n  Hybrid samples with Aldy 4 discordant for EVERY reference under BWA-MEM: {len(aldy_zero_bwa)}")
    for c, s in aldy_zero_bwa:
        print(f"    {c} {s:9} {gold[c][s]:26} CyriPanel BWA {pct(rate(cyri[(c, s)], 'prim')):.1f}%")
    seven = [(c, s) for c, s in aldy_zero_bwa if rate(cyri[(c, s)], "prim") > 0]
    print(f"  ...of which CyriPanel was concordant in at least one evaluation: {len(seven)}")
    for c, s in seven:
        cores = {int(x) for x in re.findall(r"\*(\d+)", gold[c][s])}
        print(f"    {s:9} {gold[c][s]:26} has*36={36 in cores} has*68={68 in cores} "
              f"tandem={multicopy_hybrid_tandem(gold[c][s])}  Cyri BWA {pct(rate(cyri[(c, s)], 'prim')):.1f}%  "
              f"Aldy DRAGEN {pct(rate(aldy[(c + '_dragen', s)], 'conc')):.1f}%  Cyri DRAGEN {pct(rate(cyri[(c + '_dragen', s)], 'prim')):.1f}%")
    if seven:
        cr = [pct(rate(cyri[x], "prim")) for x in seven]
        print(f"  CyriPanel BWA range on these: {min(cr):.1f}% to {max(cr):.1f}%")
    for s in ("HG02071", "HG02056"):
        c = next(cc for cc, ss in strata if ss == s)
        print(f"  {s}: Aldy BWA {pct(rate(aldy[(c, s)], 'conc')):.1f}%  Cyri BWA {pct(rate(cyri[(c, s)], 'prim')):.1f}%  |  "
              f"Aldy DRAGEN {pct(rate(aldy[(c + '_dragen', s)], 'conc')):.1f}%  Cyri DRAGEN {pct(rate(cyri[(c + '_dragen', s)], 'prim')):.1f}%")

    # ------------------------------------------------ Aldy reference choice
    print()
    print("=" * 88)
    print("ALDY 4 CONCORDANCE PER REFERENCE  (value: number of references giving it)")
    print("=" * 88)
    for arm in ("cohort1", "cohort1_dragen", "cohort2", "cohort2_dragen"):
        for strat in ("all", "no_SV", "CNV_only", "hybrid"):
            d = per_ref.get((arm, strat))
            if not d:
                continue
            vals = collections.Counter((v[0], v[1]) for v in d.values())
            desc = "  ".join(f"{k}/{t}={pct(k / t):.1f}% x{m}" for (k, t), m in sorted(vals.items(), key=lambda kv: -kv[1]))
            minority = sorted(r for r, v in d.items() if (v[0], v[1]) != vals.most_common(1)[0][0])
            print(f"  {arm:15}{strat:9} refs={len(d):<3} {desc}")
            if minority:
                print(f"  {'':24}minority value from: {', '.join(minority)}")

    # --------------------------------------------- no-calls and candidates
    print()
    print("=" * 88)
    print("NO-CALLS AND TWO-CANDIDATE OUTPUTS  (CyriPanel, panel size 5)")
    print("=" * 88)
    for arm in ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen"):
        for s in sorted(gold[arm.replace("_dragen", "")]):
            d = cyri[(arm, s)]
            if d["nocall"]:
                print(f"  no-call  {arm:15}{s:9} {d['nocall']}/{d['n']} = {pct(d['nocall'] / d['n']):.1f}% of its evaluations")
    amb = collections.defaultdict(list)
    for arm in ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen"):
        for s in sorted(gold[arm.replace("_dragen", "")]):
            d = cyri[(arm, s)]
            if d["amb"]:
                amb[s].append(arm)
                print(f"  multi    {arm:15}{s:9} {d['amb']}/{d['n']} evaluations, candidates per output {sorted(d['amb_ncands'])}, "
                      f"reference standard among candidates in {d['amb_gold_in']}/{d['amb']}, stratum {strata[(arm.replace('_dragen', ''), s)]}")
    c1 = sorted(s for s in amb if (("cohort1", s) in strata))
    print(f"  samples with multiple candidates: {len(amb)}; in cohort 1: {len(c1)} ({', '.join(c1)})")

    # ------------------------------------------- values in order of the draft
    print()
    print("=" * 88)
    print("[N] VALUES IN ORDER OF APPEARANCE  (differences are CyriPanel minus Aldy 4, pp)")
    print("=" * 88)
    P = lambda k, a, s: res[(k, a, s, None)]
    seven_c = [pct(rate(cyri[x], "prim")) for x in seven]
    six = [x for x in seven if x[1] != "HG02071"]
    six_aldy_dragen = [pct(rate(aldy[(c + "_dragen", s)], "conc")) for c, s in six]
    hg02071 = next(c for c, s in strata if s == "HG02071")
    hg00232 = cyri[("cohort2", "HG00232")]
    vals = [
        ("P1 CyriPanel, BWA-MEM, all", P("prim", "BWA-MEM", "all")[1]),
        ("P1 Aldy 4, BWA-MEM, all", P("prim", "BWA-MEM", "all")[2]),
        ("P1 difference, BWA-MEM, all", P("prim", "BWA-MEM", "all")[3]),
        ("P1 CyriPanel, DRAGEN, all", P("prim", "DRAGEN", "all")[1]),
        ("P1 Aldy 4, DRAGEN, all", P("prim", "DRAGEN", "all")[2]),
        ("P1 difference, DRAGEN, all", P("prim", "DRAGEN", "all")[3]),
        ("P2 CyriPanel, BWA-MEM, hybrid", P("prim", "BWA-MEM", "hybrid")[1]),
        ("P2 Aldy 4, BWA-MEM, hybrid", P("prim", "BWA-MEM", "hybrid")[2]),
        ("P2 difference, BWA-MEM, hybrid", P("prim", "BWA-MEM", "hybrid")[3]),
        ("P2 CyriPanel, DRAGEN, hybrid", P("prim", "DRAGEN", "hybrid")[1]),
        ("P2 Aldy 4, DRAGEN, hybrid", P("prim", "DRAGEN", "hybrid")[2]),
        ("P2 difference, DRAGEN, hybrid", P("prim", "DRAGEN", "hybrid")[3]),
        ("P3 CyriPanel BWA-MEM on the seven, minimum", min(seven_c)),
        ("P3 CyriPanel BWA-MEM on the seven, maximum", max(seven_c)),
        ("P3 Aldy 4 DRAGEN on the other six, maximum", max(six_aldy_dragen)),
        ("P3 HG02071 Aldy 4, DRAGEN", pct(rate(aldy[(hg02071 + "_dragen", "HG02071")], "conc"))),
        ("P3 HG02071 CyriPanel, DRAGEN", pct(rate(cyri[(hg02071 + "_dragen", "HG02071")], "prim"))),
        ("P4 difference, BWA-MEM, no_SV", P("prim", "BWA-MEM", "no_SV")[3]),
        ("P4 difference, DRAGEN, no_SV", P("prim", "DRAGEN", "no_SV")[3]),
        ("P4 difference, BWA-MEM, CNV_only", P("prim", "BWA-MEM", "CNV_only")[3]),
        ("P4 difference, DRAGEN, CNV_only", P("prim", "DRAGEN", "CNV_only")[3]),
        ("P6 HG00232 no-call share of its evaluations", pct(hg00232["nocall"] / hg00232["n"])),
        ("P6 samples with two candidate diplotypes", len(amb)),
        ("P6 single-candidate difference, BWA-MEM, all", P("strict", "BWA-MEM", "all")[3]),
        ("P6 single-candidate difference, DRAGEN, all", P("strict", "DRAGEN", "all")[3]),
        ("P6 single-candidate difference, BWA-MEM, no_SV", P("strict", "BWA-MEM", "no_SV")[3]),
        ("P6 single-candidate difference, DRAGEN, no_SV", P("strict", "DRAGEN", "no_SV")[3]),
        ("P6 single-candidate difference, BWA-MEM, hybrid", P("strict", "BWA-MEM", "hybrid")[3]),
        ("P6 single-candidate difference, DRAGEN, hybrid", P("strict", "DRAGEN", "hybrid")[3]),
        ("P6 two-candidate samples in cohort 1", len(c1)),
        ("P6 two-candidate samples in total", len(amb)),
    ]
    for i, (lbl, v) in enumerate(vals, 1):
        shown = f"{v}" if isinstance(v, int) else f"{v:.1f}"
        print(f"  [{i:>2}] {lbl:52} {shown:>7}")

    print()
    print("=" * 88)
    print("FAILED CHECKS: " + (", ".join(FAILS) if FAILS else "none"))
    print("=" * 88)


if __name__ == "__main__":
    main()
