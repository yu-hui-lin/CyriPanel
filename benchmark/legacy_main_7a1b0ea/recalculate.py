"""Recalculate descriptive concordance from archived evaluations; no calling rerun.

Uses only the Python standard library. Reads published_results/per_sample_v2.csv and the
reference standard in data/. SOURCE_SHA256 identifies the per_sample_v2.csv of release
v1.0.2, and provenance.json records the hashes of all inputs. Run this file without
arguments to write the outputs to recalculate_output/.

This is the script used for the preprint, adapted for release v1.0.2. Its inputs used to
be copies in inputs/ next to the script, taken from release v1.0.0, and its outputs were
written beside it. v1.0.2 renamed two columns of per_sample_v2.csv, gold_v2 to gold and
comparison_v2 to comparison_single, which also changed the file hash. The output labels
now follow the preprint's table numbering. The calculations are unchanged.
"""

import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import statistics
import sys

ROOT = Path(__file__).resolve().parent
PER_SAMPLE = ROOT / "published_results" / "per_sample_v2.csv"
DATA = ROOT / "data"
OUT = ROOT / "recalculate_output"
COHORTS = ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen")
SIZES = (5, 10, 20)
SOURCE_SHA256 = "897d8286610c564379c4282a792f570ddc1b15a29f93d6742db7c32367b40e41"


def numeric_key(value):
    return tuple(int(x) for x in re.findall(r"\d+", value))


def normalize(value):
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


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    source = PER_SAMPLE
    gold_files = [DATA / f"CYP2D6_{c}_gold_standard_v2.csv" for c in ("cohort1", "cohort2")]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == SOURCE_SHA256
    truth = {
        c: {r["Sample_ID"]: r["Genotype"] for r in read_csv(path)}
        for c, path in zip(("cohort1", "cohort2"), gold_files)
    }
    rows = read_csv(source)
    assert len(rows) == 16680
    keys = set()
    groups = collections.defaultdict(list)
    label_mismatches = []
    ambiguous_total = 0
    for row in rows:
        cohort, sample = row["cohort"], row["sample_id"]
        panel, iteration = int(row["panel_size"]), int(row["iteration"])
        key = (cohort, panel, iteration, sample)
        assert key not in keys
        keys.add(key)
        gold = truth[cohort.replace("_dragen", "")][sample]
        assert gold == row["gold"]
        called = normalize(row["called_genotype"])
        target = normalize(gold)
        row["evaluable"] = target is not None
        candidates = [normalize(x) for x in row["called_genotype"].split(";") if normalize(x) is not None]
        row["strict"] = row["evaluable"] and called == target
        row["primary"] = row["evaluable"] and target in candidates
        row["ambiguous"] = len(candidates) > 1
        ambiguous_total += row["ambiguous"]
        expected_label = "no_call" if called is None else ("concordant" if row["strict"] else "discordant")
        if expected_label != row["comparison_single"]:
            label_mismatches.append(key)
        groups[(cohort, panel, iteration)].append(row)
    assert not label_mismatches
    assert len(groups) == 360

    per_iteration, summary, common_samples, sample_rates, pairwise = [], [], [], [], []
    findings = {}
    for cohort in COHORTS:
        cohort_rows = [r for r in rows if r["cohort"] == cohort]
        ids = {r["sample_id"] for r in cohort_rows}
        evaluable_sets = [{r["sample_id"] for r in cohort_rows if int(r["panel_size"]) == n and r["evaluable"]} for n in SIZES]
        common = set.intersection(*evaluable_sets)
        findings[cohort] = {"n_common": len(common), "excluded": sorted(ids - common)}
        for sample in sorted(ids):
            entry = {"cohort": cohort, "sample_id": sample, "in_common_set": sample in common}
            for n in SIZES:
                entry[f"evaluable_iterations_n{n}"] = sum(r["evaluable"] for r in cohort_rows if r["sample_id"] == sample and int(r["panel_size"]) == n)
            common_samples.append(entry)
        rates = {}
        for n in SIZES:
            for iteration in range(30):
                group = groups[cohort, n, iteration]
                expected_n = (44 if cohort.startswith("cohort1") else 72) - n
                assert len(group) == expected_n
                ev = [r for r in group if r["evaluable"]]
                restricted = [r for r in ev if r["sample_id"] in common]
                per_iteration.append({
                    "cohort": cohort, "panel_size": n, "iteration": iteration,
                    "n_test": len(group), "n_evaluable": len(ev),
                    "n_primary_concordant": sum(r["primary"] for r in ev),
                    "n_strict_concordant": sum(r["strict"] for r in ev),
                    "primary_concordance_pct": statistics.mean(r["primary"] for r in ev) * 100,
                    "strict_concordance_pct": statistics.mean(r["strict"] for r in ev) * 100,
                    "common_iteration_concordance_pct": statistics.mean(r["primary"] for r in restricted) * 100,
                    "pass_rate_pct": statistics.mean(r["filter"] == "PASS" for r in ev) * 100,
                    "n_no_call": sum(normalize(r["called_genotype"]) is None for r in group),
                    "n_ambiguous": sum(r["ambiguous"] for r in ev),
                })
            rates[n] = {}
            for sample in sorted(common):
                ev = [r for r in cohort_rows if int(r["panel_size"]) == n and r["sample_id"] == sample and r["evaluable"]]
                rate = statistics.mean(r["primary"] for r in ev)
                rates[n][sample] = rate
                sample_rates.append({"cohort": cohort, "panel_size": n, "sample_id": sample,
                                     "n_evaluable": len(ev), "n_concordant": sum(r["primary"] for r in ev),
                                     "concordance_pct": rate * 100})
            series = [x for x in per_iteration if x["cohort"] == cohort and x["panel_size"] == n]
            record = {"cohort": cohort, "panel_size": n, "n_iterations": len(series), "n_common": len(common)}
            for col in ("primary_concordance_pct", "strict_concordance_pct", "common_iteration_concordance_pct", "pass_rate_pct"):
                record[col + "_mean"] = statistics.mean(x[col] for x in series)
                record[col + "_sd"] = statistics.stdev(x[col] for x in series)
            record["no_call_count_mean"] = statistics.mean(x["n_no_call"] for x in series)
            record["ambiguous_count_mean"] = statistics.mean(x["n_ambiguous"] for x in series)
            record["equal_sample_mean_concordance_pct"] = statistics.mean(rates[n].values()) * 100
            summary.append(record)
        findings[cohort]["samples_with_varying_rates"] = [s for s in sorted(common) if len({rates[n][s] for n in SIZES}) > 1]
        findings[cohort]["always_discordant_common_samples"] = [s for s in sorted(common) if all(rates[n][s] == 0 for n in SIZES)]
        for a, b in itertools.combinations(SIZES, 2):
            diffs = [100 * (rates[b][s] - rates[a][s]) for s in sorted(common)]
            pairwise.append({"cohort": cohort, "smaller_panel": a, "larger_panel": b,
                             "n_common": len(common), "median_paired_difference_pp": statistics.median(diffs),
                             "n_nonzero_differences": sum(d != 0 for d in diffs)})

    for name, data in (("iteration_summary.csv", per_iteration), ("condition_summary.csv", summary),
                       ("common_samples.csv", common_samples), ("common_sample_rates.csv", sample_rates),
                       ("paired_descriptive_summary.csv", pairwise)):
        write_csv(name, data)
    provenance = {"source": "published_results/ and data/ of release v1.0.2", "python": sys.version,
                  "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [source] + gold_files},
                  "n_records": len(rows), "n_iterations": len(groups), "n_ambiguous_evaluations": ambiguous_total,
                  "strict_label_mismatches": len(label_mismatches), "findings": findings,
                  "inferential_tests": "Not recomputed. This output contains descriptive summaries only."}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    header = ["# Table 5. CyriPanel concordance under candidate-inclusive scoring", "",
              "Concordance counts a call as correct when any reported candidate matches the reference-standard diplotype. Values are mean ± SD across 30 iterations. Concordance and PASS rates are computed over every evaluation with a reference standard, so a no-call is retained and scored as a failure rather than excluded from the denominator. No-call counts per iteration are also reported.", "",
              "| Cohort | Reference-panel size | BWA-MEM concordance (%) | DRAGEN concordance (%) | PASS rate (%) | No-calls per iteration |",
              "|---|---|---|---|---|---|"]
    lookup = {(x["cohort"], x["panel_size"]): x for x in summary}
    for c in (1, 2):
        for n in SIZES:
            bwa, dragen = lookup[f"cohort{c}", n], lookup[f"cohort{c}_dragen", n]
            assert abs(bwa["pass_rate_pct_mean"] - dragen["pass_rate_pct_mean"]) < 1e-10
            assert bwa["no_call_count_mean"] == dragen["no_call_count_mean"]
            fm = lambda x: f'{x["primary_concordance_pct_mean"]:.1f} ± {x["primary_concordance_pct_sd"]:.1f}'
            header.append(f'| {c} | {n} | {fm(bwa)} | {fm(dragen)} | {bwa["pass_rate_pct_mean"]:.1f} ± {bwa["pass_rate_pct_sd"]:.1f} | {bwa["no_call_count_mean"]:.2f} |')
    header += ["", "## Common-sample summaries", "",
               "The common set requires at least one evaluable record at each panel size. The equally weighted mean below averages the per-sample concordance rates, assigning each common sample equal weight.", "",
               "| Cohort | Aligner | Common samples | n = 5 (%) | n = 10 (%) | n = 20 (%) |",
               "|---|---|---|---|---|---|"]
    for cohort in COHORTS:
        values = [lookup[cohort, n]["equal_sample_mean_concordance_pct"] for n in SIZES]
        header.append(f'| {cohort[6]} | {"DRAGEN" if "dragen" in cohort else "BWA-MEM"} | {findings[cohort]["n_common"]} | ' + " | ".join(f"{v:.2f}" for v in values) + " |")
    header += ["", "These are recalculated descriptive results. The Friedman and Wilcoxon tests are in common_sample_analysis.py. Every saved output was retained, with unresolved and empty outputs alike scored as discordant evaluations, so no sample is excluded on the basis of the caller's own output.", ""]
    (OUT / "Table5_candidate_inclusive.md").write_text("\n".join(header))
    print(json.dumps({"integrity": {"rows": len(rows), "strict_label_mismatches": len(label_mismatches)}, "findings": findings}, indent=2))


if __name__ == "__main__":
    main()
