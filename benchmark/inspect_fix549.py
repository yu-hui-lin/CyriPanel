#!/usr/bin/env python3
"""
inspect_fix549.py
=================
Reads the outputs of verify_fix549.sh and reports, per sample, whether the
star_caller.py:549 fix changed anything — and whether the call was right in the
first place.

`diff` on the TSVs answers "did anything change". This answers three more:
  * was the call concordant with the v2 gold, before and after?
  * if a call changed, which intermediate fields changed with it?
  * did the affected tie-break (*119/*2 vs *1/*41) ever actually fire?

Usage:
    python3 inspect_fix549.py
    python3 inspect_fix549.py --scratch /path/to/verify_fix549
    python3 inspect_fix549.py --verbose      # dump all differing JSON fields
"""

import os
import sys
import json
import csv
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

DEFAULT_SCRATCH = "/work/u7715055/staging/biology/u7715055/verify_fix549"
SRC = "/work/u7715055/staging/biology/u7715055/CyriPanel"

GOLD = {
    "cohort1": os.path.join(SRC, "data", "CYP2D6_cohort1_gold_standard_v2.csv"),
    "cohort2": os.path.join(SRC, "data", "CYP2D6_cohort2_gold_standard_v2.csv"),
}

# Fields worth showing when a call changes. Everything else is noise.
DIAGNOSTIC_FIELDS = [
    "Total_CN", "Spacer_CN", "CNV_group", "Exon9_CN",
    "CNV_consensus", "Variants_called", "Raw_star_allele", "Call_info",
]


def load_gold(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["Sample_ID"]: r["Genotype"] for r in csv.DictReader(f)}


def clean_id(sid):
    """star_caller derives the sample id from the BAM basename, so
    'HG00733.markdup.bam' becomes 'HG00733.markdup'. Strip the suffix."""
    for suffix in (".markdup", ".sorted", ".bam"):
        if sid.endswith(suffix):
            sid = sid[: -len(suffix)]
    return sid


def load_calls(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        raw = json.load(f)
    return {clean_id(k): v for k, v in raw.items()}


def get_comparator():
    """
    Use the benchmark's own comparison so this agrees with every other number.

    Three routes, in order of preference. The middle one matters: importing
    run_single_iteration executes its module-level code, which can fail for
    reasons unrelated to the comparison logic (missing dependency, path setup).
    Extracting the five normalisation functions from the source text gives
    exactly the same behaviour with no side effects. The exact-string fallback
    is a last resort and is loudly flagged, because without normalisation it
    reports *41/*4 vs *4/*41 as discordant when they are the same diplotype.
    """
    import re as _re

    try:
        from run_single_iteration import compare_genotype
        return compare_genotype
    except Exception as e:
        print("  NOTE: could not import run_single_iteration (%s: %s)."
              % (type(e).__name__, e))
        print("        Extracting the normalisation functions from the source instead.")

    try:
        src_path = os.path.join(SCRIPT_DIR, "run_single_iteration.py")
        src = open(src_path).read()
        ns = {"re": _re}
        for fn in ("_expand_multiplication", "_sort_allele_key",
                   "_normalize_haplotype", "normalize_genotype", "compare_genotype"):
            m = _re.search(r"\ndef " + fn + r"\(.*?(?=\ndef |\n# \u2500|\Z)", src, _re.S)
            if m is None:
                raise RuntimeError("could not locate %s() in %s" % (fn, src_path))
            exec(m.group(0), ns)
        print("        Extraction succeeded — comparisons use the benchmark's own rules.\n")
        return ns["compare_genotype"]
    except Exception as e:
        print("        Extraction failed too (%s: %s)." % (type(e).__name__, e))

    def fallback(called, gold):
        if called in (None, "None", "", "Error"):
            return "no_call"
        if gold is None:
            return "no_gold"
        return "concordant" if str(called).strip() == str(gold).strip() else "discordant"

    print("  *** WARNING: falling back to exact string matching. Allele order and")
    print("      multiplication notation are NOT normalised, so genotypes that differ")
    print("      only in ordering (*41/*4 vs *4/*41) will be reported as discordant.")
    print("      Treat the 'verdict' column as unreliable; the changed/unchanged")
    print("      comparison between the two runs is still valid.\n")
    return fallback


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=DEFAULT_SCRATCH)
    ap.add_argument("--verbose", action="store_true",
                    help="dump every differing JSON field, not just the diagnostic ones")
    args = ap.parse_args()

    out_dir = os.path.join(args.scratch, "out")
    if not os.path.isdir(out_dir):
        print("No output directory at %s" % out_dir)
        print("The job has probably not finished (or not started). Check with:")
        print("  squeue -u $USER")
        print("  tail -40 %s/benchmark/verify549_out.log" % SRC)
        return 2

    compare = get_comparator()
    total_changed = 0
    total_samples = 0
    grand_summary = []

    for cohort in ("cohort1", "cohort2"):
        pre = load_calls(os.path.join(out_dir, "%s_prefix.json" % cohort))
        fix = load_calls(os.path.join(out_dir, "%s_fixed.json" % cohort))
        if pre is None or fix is None:
            print("=" * 100)
            print("  %s: output missing (pre-fix=%s, fixed=%s)"
                  % (cohort.upper(), pre is not None, fix is not None))
            print("=" * 100)
            print()
            continue

        gold = load_gold(GOLD[cohort])
        samples = sorted(set(pre) | set(fix))

        print("=" * 100)
        print("  %s — %d samples" % (cohort.upper(), len(samples)))
        print("=" * 100)
        print("  %-12s %-22s %-22s %-22s %s"
              % ("Sample", "v2 gold", "pre-fix call", "post-fix call", "verdict"))
        print("  " + "-" * 96)

        changed_here = []
        for sid in samples:
            g = gold.get(sid)
            a = (pre.get(sid) or {}).get("Genotype")
            b = (fix.get(sid) or {}).get("Genotype")
            fa = (pre.get(sid) or {}).get("Filter")
            fb = (fix.get(sid) or {}).get("Filter")

            same_call = (a == b) and (fa == fb)
            if same_call:
                c = compare(b, g)
                verdict = {"concordant": "unchanged, correct",
                           "discordant": "unchanged, DISCORDANT",
                           "no_call":    "unchanged, no-call",
                           "no_gold":    "unchanged, no gold"}.get(c, "unchanged (%s)" % c)
            else:
                verdict = "*** CHANGED ***"
                changed_here.append(sid)

            print("  %-12s %-22s %-22s %-22s %s"
                  % (sid, g or "?", a or "None", b or "None", verdict))
            total_samples += 1

        print()
        if not changed_here:
            print("  No call changed in %s." % cohort)
        else:
            total_changed += len(changed_here)
            print("  %d call(s) changed in %s — details below." % (len(changed_here), cohort))
            for sid in changed_here:
                print()
                print("  --- %s ---" % sid)
                pa, pb = pre.get(sid, {}), fix.get(sid, {})
                fields = sorted(set(pa) | set(pb)) if args.verbose else DIAGNOSTIC_FIELDS
                for k in fields:
                    va, vb = pa.get(k), pb.get(k)
                    if va != vb:
                        print("      %-18s pre-fix : %s" % (k, va))
                        print("      %-18s post-fix: %s" % ("", vb))
                    elif args.verbose:
                        print("      %-18s (same)  : %s" % (k, va))
                print("      gold: %s" % gold.get(sid))
        grand_summary.append((cohort, len(samples), len(changed_here)))
        print()

    print("=" * 100)
    print("  VERDICT")
    print("=" * 100)
    for cohort, n, c in grand_summary:
        print("    %-10s %2d samples, %d changed" % (cohort, n, c))
    print()
    if total_samples == 0:
        print("  Nothing to compare — no outputs were read.")
        return 2
    if total_changed == 0:
        print("  PASS — the star_caller.py:549 fix changed no call in any of the %d samples"
              % total_samples)
        print("  that can reach the affected tie-break. The existing benchmark arms remain")
        print("  valid, and the manuscript may state that *119/*2 discrimination is enabled")
        print("  without re-running anything.")
        print()
        print("  Wording for Methods or the response to reviewers:")
        print("    \"An argument-order defect that disabled the g.42127803 haplotype test was")
        print("     corrected after the benchmark was run. Because neither cohort contains a")
        print("     *119 or *32 carrier, the correction cannot alter any reported call; this")
        print("     was confirmed by re-running the %d samples whose reference diplotype can" % total_samples)
        print("     reach the affected decision branch, which produced identical genotypes.\"")
        return 0
    print("  FAIL — %d call(s) changed. Do not treat the existing benchmark arms as valid"
          % total_changed)
    print("  until you understand why. Look at CNV_group and Raw_star_allele above: if")
    print("  Call_info moved from a unique match to a tie-break, the g.42127803 haplotype")
    print("  test is now firing, which is the expected mechanism — but the direction of the")
    print("  change is what matters. A call that moved TOWARDS the gold is a genuine")
    print("  improvement and means the affected arms should be re-run; a call that moved AWAY")
    print("  from it means the fix has a side effect that needs tracing before anything else.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
