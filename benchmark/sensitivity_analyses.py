#!/usr/bin/env python3
"""Audit sensitivity analyses from pinned archived calls; does not rerun callers.

Run: python -I benchmark/sensitivity_analyses.py --outdir NEW_DIRECTORY
Requires reproduce_tables.py beside this script. Standard library only.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import sys

HELPER_SHA = "8f51f22ad5b397bb62c92424784dd4ea8e7a86e5a31cae4f5d61261fc62a0d04"
DISPLAY_SHA = "d9a0bca70f9958e9a2fdec83a01ce15da07426fb4b6d34c1c9705a92b7094f88"
CATEGORIES = ("no_call", "phased_copy_distribution", "unphased_total_copy",
              "phased_hybrid_status", "unphased_unresolved", "other_unresolved")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def alleles(value):
    return set(map(int, re.findall(r"\*(\d+)", value)))


def structure(value):
    # Ordinary star identities are collapsed. *5 has zero copies.
    # Each haplotype is represented as (normal copies, hybrid copies).
    if not value:
        return None
    phased = "_" not in value
    groups = value.split("/") if phased else [value]
    require(len(groups) == (2 if phased else 1), "Unexpected phase notation: " + value)
    counts = []
    for group in groups:
        normal = hybrid = 0
        for token in re.split(r"[+_]", group):
            match = re.fullmatch(r"\*(\d+)(?:x(\d+))?", token.strip())
            require(match is not None, "Unsupported allele token: " + token)
            star, copies = int(match[1]), int(match[2] or 1)
            require(copies > 0, "Nonpositive gene-copy count.")
            if star in {13, 36, 68, 83}:
                hybrid += copies
            elif star != 5:
                normal += copies
        counts.append((normal, hybrid))
    return phased, tuple(sorted(counts))


def classify(call, gold):
    if call is None:
        return "no_call"
    phased, counts = call
    _, target = gold
    if not phased:
        return ("unphased_total_copy" if sum(map(sum, counts)) != sum(map(sum, target))
                else "unphased_unresolved")
    if sorted(map(sum, counts)) != sorted(map(sum, target)):
        return "phased_copy_distribution"
    return "phased_hybrid_status" if counts != target else "other_unresolved"


def csv_text(header, rows):
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows([header, *rows])
    return stream.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    require(not args.outdir.exists() and not args.outdir.is_symlink(), "Output path already exists.")
    helper = Path(__file__).resolve().with_name("reproduce_tables.py")
    require(sha(helper.read_bytes()) == HELPER_SHA, "reproduce_tables.py checksum mismatch.")
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("table_helpers", helper)
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    data = {}
    for name, (count, digest) in h.INPUTS.items():
        raw = (args.benchmark / name).read_bytes()
        require(sha(raw) == digest, "Input checksum mismatch: " + name)
        data[name] = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        require(len(data[name]) == count, "Input row count mismatch: " + name)
    gold = {c: {r["Sample_ID"]: r["Genotype"] for r in
                data[f"data/CYP2D6_{c}_gold_standard_v2.csv"]} for c in h.ARMS[:2]}
    strata = {s: h.stratum(g) for cohort in gold.values() for s, g in cohort.items()}
    carrier49 = {s for cohort in gold.values() for s, g in cohort.items() if 49 in alleles(g)}
    hybrid49 = {s for s in carrier49 if strata[s] == "Hybrid"}
    require(len(strata) == 116 and len(carrier49) == 4 and
            hybrid49 == {"HG00597", "HG00609", "NA24631"}, "Unexpected sample membership.")
    cyri, aldy = defaultdict(list), defaultdict(list)
    details, classes = [], Counter()
    by_sample, carriers = defaultdict(Counter), defaultdict(Counter)
    false49 = eligible49 = 0
    for row in data["published_results/per_sample_v2.csv"]:
        arm, sid = row["cohort"], row["sample_id"]
        truth = gold[arm.replace("_dragen", "")][sid]
        require(row["gold_v2"] == truth, "CyriPanel truth mismatch.")
        candidates = [v for part in row["called_genotype"].split(";")
                      if (v := h.normalize(part)) is not None]
        correct = h.normalize(truth) in candidates
        strict = correct and len(candidates) == 1
        scenario = correct and sid not in hybrid49
        if int(row["panel_size"]) == 5:
            cyri[arm, sid].append((correct, strict, scenario))
        if sid not in carrier49:
            eligible49 += 1
            false49 += int(49 in alleles(row["called_genotype"]))
        if sid in hybrid49:
            carriers[arm, sid]["evaluations"] += 1
            carriers[arm, sid]["correct"] += int(correct)
        if not correct:
            require(len(candidates) <= 1, "Unexpected multiple-candidate discordance.")
            call = structure(candidates[0]) if candidates else None
            target = structure(truth)
            category = classify(call, target)
            classes[category] += 1
            by_sample[sid][category] += 1
            details.append([arm, sid, row["panel_size"], row["iteration"], truth,
                            row["called_genotype"], row["filter"], category,
                            json.dumps(call, separators=(",", ":")),
                            json.dumps(target, separators=(",", ":"))])
    for row in data["aldy/aldy_results_long.csv"]:
        arm, sid = row["arm"], row["sample"]
        target = h.normalize(gold[arm.replace("_dragen", "")][sid])
        call = h.normalize(h.clean_aldy(row["raw_major"]))
        require(row["n_solutions"] == "1" and call == row["normalised"] and
                target == row["gold_norm"], "Aldy retained-call mismatch.")
        aldy[arm, sid].append(call == target)
    expected_keys = {(a, s) for a in h.ARMS for s in gold[a.replace("_dragen", "")]}
    require(set(cyri) == set(aldy) == expected_keys, "Missing sample/aligner rates.")
    rates, sample_rows = {}, []
    for arm, sid in sorted(expected_keys):
        values, acalls = cyri[arm, sid], aldy[arm, sid]
        counts = [sum(v[i] for v in values) for i in range(3)]
        fractions = [Fraction(n, len(values)) for n in counts] + [Fraction(sum(acalls), len(acalls))]
        rates[arm, sid] = fractions
        sample_rows.append([arm, sid, strata[sid], len(values), *counts, len(acalls),
                            sum(acalls), *map(str, fractions)])
    display = []
    for label, suffix in (("BWA-MEM", ""), ("DRAGEN", "_dragen")):
        for group in ("All samples", "No structural variation", "Copy-number variation only", "Hybrid"):
            keys = [(c + suffix, s) for c in h.ARMS[:2] for s in gold[c]
                    if group == "All samples" or strata[s] == group]
            means = [sum((rates[k][i] for k in keys), Fraction()) * 100 / len(keys) for i in range(4)]
            display.append([label, group, str(len(keys)), *[f"{float(v):.4f}" for v in means]])
    require(sha(json.dumps(display, ensure_ascii=False, separators=(",", ":")).encode()) == DISPLAY_SHA,
            "Sensitivity rows differ from independent reconstruction.")
    require([classes[c] for c in CATEGORIES] == [106, 430, 172, 270, 8, 0] and len(details) == 986,
            "Discordance categories differ from independent reconstruction.")
    require(eligible49 == 15960 and false49 == 0, "Unexpected false-*49 call audit.")
    unknown = {(r[0], r[1], int(r[2]), int(r[3])) for r in details if r[7] == "unphased_unresolved"}
    require(unknown == {("cohort2", "HG00597", 5, i) for i in (2, 10, 11, 14, 16, 23, 24)} |
            {("cohort2", "HG00597", 10, 26)}, "Unexpected unresolved evaluations.")
    files = {}
    files["sensitivity.csv"] = csv_text(
        ["aligner", "stratum", "samples", "cyri_candidate_inclusive_pct", "cyri_single_candidate_pct",
         "cyri_star49_recode_pct", "aldy_retained_major_pct"], display)
    files["per_sample_rates.csv"] = csv_text(
        ["arm", "sample", "stratum", "cyri_evaluations", "cyri_primary_correct", "cyri_single_correct",
         "cyri_star49_recode_correct", "aldy_evaluations", "aldy_correct", "primary_fraction",
         "single_fraction", "star49_recode_fraction", "aldy_fraction"], sample_rows)
    files["discordance_evaluations.csv"] = csv_text(
        ["arm", "sample", "panel_size", "iteration", "gold", "call", "filter", "category",
         "call_structure", "gold_structure"], details)
    files["discordance_by_sample.csv"] = csv_text(
        ["sample", *CATEGORIES, "total"],
        [[s, *[by_sample[s][c] for c in CATEGORIES], sum(by_sample[s].values())] for s in sorted(by_sample)])
    files["star49_hybrid_samples.csv"] = csv_text(
        ["arm", "sample", "gold", "all_panel_evaluations", "primary_correct", "hypothetical_correct"],
        [[a, s, gold[a.replace("_dragen", "")][s], counts["evaluations"], counts["correct"], 0]
         for (a, s), counts in sorted(carriers.items())])
    files["summary.md"] = """# Sensitivity analyses from archived calls

`sensitivity.csv` uses panel size 5 and equal weighting of sample-level rates.
All evaluated observations, including no-calls, remain in denominators.
Primary scoring accepts any matching candidate. Single-candidate scoring requires
exactly one candidate and a match. Exact rational rates are in `per_sample_rates.csv`.
Aldy scores refer only to the retained major calls in the archived CSV; this audit
does not verify candidate selection against original Aldy solution files.

The *49 scenario forces every evaluation of HG00597, HG00609 and NA24631 to be
incorrect while preserving denominators. This is hypothetical rescoring, not a
filter-removal experiment or a demonstrated worst-case bound on filter effects.
Across all panel sizes and both aligners, no *49 call was observed in 15,960
archived evaluations of the 112 gold-non-*49 samples (no-calls included).

Across all 16,680 CyriPanel evaluations, candidate-inclusive discordance comprises:

- 106 no-calls.
- 430 phased calls with different per-haplotype gene-copy counts.
- 172 unphased calls with a different total gene-copy count.
- 270 phased calls with equal per-haplotype copy counts but different normal/hybrid counts.
- 8 unphased HG00597 calls whose total copy count does not resolve haplotype structure.

The structure signature is [phased, sorted [normal_count, hybrid_count] groups].
*5 contributes zero copies; *13, *36, *68 and *83 are collapsed to one hybrid class.
Other star identities are ignored. These coarse categories do not establish
biological equivalence, a causal mechanism, or validated structural truth.
The eight unresolved calls cannot support a claim that all discordance is structural.
This audit re-scores archived results; it does not rerun genotyping or CNV estimation.
In INPUT_SHA256SUMS, CSV paths are relative to --benchmark; script paths are
relative to the directory containing sensitivity_analyses.py.
"""
    files["validation.txt"] = (
        "SENSITIVITY_ROWS_MATCH: 8\nPER_SAMPLE_RATES: 232\nDISCORDANT_EVALUATIONS: 986\n" +
        "".join(f"{c}: {classes[c]}\n" for c in CATEGORIES) +
        "FALSE_STAR49_CALLS: 0/15960 gold-non-*49 evaluations\n" +
        "INPUT_SOURCE_COMMIT: " + h.SOURCE_COMMIT + "\n" +
        "SENSITIVITY_DISPLAY_SHA256: " + DISPLAY_SHA + "\n")
    files["INPUT_SHA256SUMS"] = "".join(d + "  " + name + "\n" for name, (_, d) in h.INPUTS.items())
    files["INPUT_SHA256SUMS"] += HELPER_SHA + "  reproduce_tables.py\n"
    files["INPUT_SHA256SUMS"] += sha(Path(__file__).read_bytes()) + "  sensitivity_analyses.py\n"
    files["SHA256SUMS"] = "".join(sha(text.encode()) + "  " + name + "\n" for name, text in sorted(files.items()))
    args.outdir.mkdir(parents=True)
    for name, text in files.items():
        (args.outdir / name).write_bytes(text.encode())
    print("SENSITIVITY_ROWS_MATCH: 8")
    print("DISCORDANCE_CLASSIFICATION:", ", ".join(f"{c}={classes[c]}" for c in CATEGORIES))
    print("FALSE_STAR49_CALLS: 0/15960 gold-non-*49 evaluations")
    print("SENSITIVITY_ANALYSES_REPRODUCED:", args.outdir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        raise SystemExit("STOP: " + str(exc))
