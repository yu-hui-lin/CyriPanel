#!/usr/bin/env python3
"""Reproduce descriptive manuscript Tables 3-5 from the archived evaluation CSVs.

Standard library only. Run from any directory:
  python reproduce_tables.py --outdir NEW_DIRECTORY
Inputs default to this script's benchmark directory; --benchmark overrides it.
This re-scores retained calls; it does not rerun callers or inferential tests,
verify original Aldy output files, or reconstruct CNVPanelizer depth matrices.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import statistics as st
import sys

INPUTS = {
    "published_results/per_sample_v2.csv": (16680, "fb83f267f6d142a17a64765e1f8e9165e89be430d72add9f8c234d5a3993b92f"),
    "aldy/aldy_results_long.csv": (6814, "5ab195d6b612148cb5a578a3237cc40b6264205067cfe715e45aed41a47a5c24"),
    "data/CYP2D6_cohort1_gold_standard_v2.csv": (44, "a9bfb3f6b2af094a1a12ff023bc2d2be4b97f1fe17975f688a2a41bf06c84e2a"),
    "data/CYP2D6_cohort2_gold_standard_v2.csv": (72, "918efa9334788228d103b28fcfd460708e09968dfbf564b8412ae5ad5525084f"),
    "data/CYP2D6_cohort1_33diploids.csv": (33, "7f841778411ac2b6c95bbec566b5a36ca188fb8713392da4ed8d9102605e3e46"),
    "data/CYP2D6_cohort2_28diploids.csv": (28, "bb4623b51a6f02a756defa28a502aaefb84d47d4b6e9ef31887a57f16d4d7040"),
}
# Expected display-row hashes, independently checked against the manuscript.
# Canonical form: UTF-8 JSON, ensure_ascii=False, separators=(",", ":").
EXPECTED = {
    "Table3": "faa9462d750c8c9e860216c787f6994d899f678ec3360af329941ed3d5456752",
    "Table4": "9e11e85e8e32a0e063107b597ba35b816e88eceacbc254d8b120d6e98ba53e81",
    "Table5": "4e39770b744c56a9a39281b9c89203faed71c1dc489224e740f20a07076f8817",
}
MANUSCRIPT_SHA256 = "fe094417392cab03c3021c557e41ae5e9a8c41889f72d5135bfc833a2eafe42b"
SOURCE_COMMIT = "5d33b820dd2ba8b874f1bbea51f411b5f081856f"
ARMS = ("cohort1", "cohort2", "cohort1_dragen", "cohort2_dragen")
SIZES = (5, 10, 20)
HYBRIDS = {13, 36, 68, 83}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def numeric(value):
    return tuple(map(int, re.findall(r"\d+", value)))


def normalize(value):
    if value is None or value.strip() in ("", "None", "Error", "NA", "nan", "NaN"):
        return None
    haps = []
    for hap in value.strip().split("/"):
        copies = []
        for part in hap.split("+"):
            part = part.strip()
            match = re.fullmatch(r"(\*\d+(?:\.\d+)?)x(\d+)", part)
            copies.extend([match[1]] * int(match[2]) if match else [part])
        haps.append("+".join(sorted(copies, key=numeric)))
    return "/".join(sorted(haps, key=numeric))


def clean_aldy(value):
    # Reproduce the archived core-star parser, retaining each listed gene copy.
    haps = []
    for hap in value.split("/"):
        copies = [m[1] for part in hap.split("+")
                  if (m := re.match(r"^(\*\d+)", part.strip()))]
        require(bool(copies), "Aldy major haplotype has no core star allele.")
        haps.append("+".join(copies))
    return "/".join(haps)


def stratum(value):
    alleles = set(map(int, re.findall(r"\*(\d+)", value)))
    if alleles & HYBRIDS:
        return "Hybrid"
    if 5 in alleles or "x" in value.lower() or "+" in value:
        return "Copy-number variation only"
    return "No structural variation"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--outdir", type=Path, required=True, help="New output directory; must not exist.")
    args = parser.parse_args()
    require(not args.outdir.exists() and not args.outdir.is_symlink(), "Output path already exists.")
    data = {}
    for name, (count, digest) in INPUTS.items():
        raw = (args.benchmark / name).read_bytes()
        require(sha(raw) == digest, "Input checksum mismatch: " + name)
        data[name] = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        require(len(data[name]) == count, "Input row count mismatch: " + name)

    gold = {c: {r["Sample_ID"]: r["Genotype"] for r in
                data[f"data/CYP2D6_{c}_gold_standard_v2.csv"]} for c in ARMS[:2]}
    refs = {c: {r["Sample_ID"] for r in data[f"data/CYP2D6_{c}_{n}diploids.csv"]}
            for c, n in (("cohort1", 33), ("cohort2", 28))}
    samples = [(c, s) for c in ARMS[:2] for s in sorted(gold[c])]
    require(Counter(stratum(gold[c][s]) for c, s in samples) ==
            {"No structural variation": 61, "Copy-number variation only": 27, "Hybrid": 28},
            "Unexpected sample strata.")
    cyri, aldy = defaultdict(list), defaultdict(list)
    groups, keys, filters = defaultdict(list), set(), {}
    for row in data["published_results/per_sample_v2.csv"]:
        arm, sample = row["cohort"], row["sample_id"]
        base = arm.replace("_dragen", "")
        n, iteration = int(row["panel_size"]), int(row["iteration"])
        key = (arm, n, iteration, sample)
        require(key not in keys, "Duplicate CyriPanel evaluation.")
        keys.add(key)
        target = normalize(gold[base][sample])
        require(row["gold_v2"] == gold[base][sample], "CyriPanel reference-standard mismatch.")
        candidates = [v for part in row["called_genotype"].split(";")
                      if (v := normalize(part)) is not None]
        row["_correct"] = target in candidates
        row["_nocall"] = not candidates
        archived = ("no_call" if not candidates else
                    "concordant" if normalize(row["called_genotype"]) == target else "discordant")
        require(row["comparison_v2"] == archived, "Archived CyriPanel strict label mismatch.")
        groups[arm, n, iteration].append(row)
        filters[key] = row["filter"]
        if n == 5:
            cyri[arm, sample].append(row["_correct"])
    require(set(groups) == {(a, n, i) for a in ARMS for n in SIZES for i in range(30)},
            "Missing or unexpected CyriPanel iteration.")
    for (arm, n, iteration), rows in groups.items():
        base = arm.replace("_dragen", "")
        observed = {r["sample_id"] for r in rows}
        excluded = set(gold[base]) - observed
        require(len(rows) == len(gold[base]) - n and len(excluded) == n and excluded <= refs[base],
                "Unexpected held-out reference set.")
    for (arm, n, iteration, sample), value in filters.items():
        if "_dragen" not in arm:
            require(filters.get((arm + "_dragen", n, iteration, sample)) == value,
                    "Aligners differ in test membership or filter status.")
    akeys = set()
    for row in data["aldy/aldy_results_long.csv"]:
        arm, sample, ref = row["arm"], row["sample"], row["reference"]
        base = arm.replace("_dragen", "")
        key = (arm, sample, ref)
        require(key not in akeys, "Duplicate Aldy evaluation.")
        akeys.add(key)
        cleaned = clean_aldy(row["raw_major"])
        call, target = normalize(cleaned), normalize(gold[base][sample])
        correct = call == target
        require(row["n_solutions"] == "1" and row["gold"] == gold[base][sample] and
                row["cleaned"] == cleaned and row["normalised"] == call and
                row["gold_norm"] == target and
                row["comparison"] == ("concordant" if correct else "discordant"),
                "Aldy retained-call rescoring differs from archived fields.")
        aldy[arm, sample].append(correct)
    expected_aldy = {(arm, s, ref) for arm in ARMS
                     for s in gold[arm.replace("_dragen", "")]
                     for ref in refs[arm.replace("_dragen", "")] if s != ref}
    require(akeys == expected_aldy, "Aldy eligible-reference schedule is incomplete.")
    require(set(cyri) == set(aldy) == {(a, s) for a in ARMS for s in gold[a.replace("_dragen", "")]},
            "A sample is missing from the caller comparison.")
    tables = {"Table3": [], "Table4": [], "Table5": []}
    for label, suffix in (("BWA-MEM", ""), ("DRAGEN", "_dragen")):
        for group in ("All samples", "No structural variation", "Copy-number variation only", "Hybrid"):
            selected = [(c + suffix, s) for c, s in samples
                        if group == "All samples" or stratum(gold[c][s]) == group]
            rates = [f"{100 * st.mean(st.mean(calls[k]) for k in selected):.1f}"
                     for calls in (cyri, aldy)]
            tables["Table3"].append([label, group, str(len(selected)), *rates])
    hybrids = [(c, s) for c, s in samples if stratum(gold[c][s]) == "Hybrid"]
    copies = lambda c, s: max(len(h.split("+")) for h in normalize(gold[c][s]).split("/"))
    for c, s in sorted(hybrids, key=lambda k: (-copies(*k), *k)):
        fractions = [f"{sum(calls[arm, s])}/{len(calls[arm, s])}"
                     for calls in (cyri, aldy) for arm in (c, c + "_dragen")]
        tables["Table4"].append([s, c[-1], normalize(gold[c][s]), str(copies(c, s)), *fractions])
    for c in ARMS[:2]:
        for n in SIZES:
            concordance, failures = [], []
            for arm in (c, c + "_dragen"):
                runs = [groups[arm, n, i] for i in range(30)]
                series = [100 * st.mean(r["_correct"] for r in rows) for rows in runs]
                concordance.append(f"{st.mean(series):.1f} ± {st.stdev(series):.1f}")
                failures.append(f"{st.mean(sum(not r['_correct'] for r in rows) for rows in runs):.2f}")
            rows = [r for i in range(30) for r in groups[c, n, i]]
            statuses = [f"{100 * st.mean(r['filter'] == status for r in rows):.1f}"
                        for status in ("PASS", "More_than_one_possible_genotype", "Not_assigned_to_haplotypes")]
            statuses.append(f"{100 * st.mean(r['_nocall'] for r in rows):.1f}")
            tables["Table5"].append([c[-1], str(n), str(len(gold[c]) - n),
                                     *concordance, *failures, *statuses])
    headers = {
        "Table3": ["Aligner", "Samples", "n", "CyriPanel (%)", "Aldy 4 (%)"],
        "Table4": ["Sample", "Cohort", "Reference-standard diplotype", "Copies on the larger haplotype",
                   "CyriPanel, BWA-MEM", "CyriPanel, DRAGEN", "Aldy 4, BWA-MEM", "Aldy 4, DRAGEN"],
        "Table5": ["Cohort", "References", "Test samples", "Concordance, BWA-MEM (%)", "Concordance, DRAGEN (%)",
                   "Discordant samples, BWA-MEM", "Discordant samples, DRAGEN", "PASS (%)",
                   "Ambiguous (%)", "Not assigned (%)", "No-call (%)"],
    }
    files, markdown = {}, []
    for name, rows in tables.items():
        digest = sha(json.dumps(rows, ensure_ascii=False, separators=(",", ":")).encode())
        require(digest == EXPECTED[name], name + " differs from the manuscript display rows.")
        stream = io.StringIO(newline="")
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows([headers[name], *rows])
        files[name + ".csv"] = stream.getvalue()
        markdown.extend(["## " + name, "", "| " + " | ".join(headers[name]) + " |",
                         "| " + " | ".join(["---"] * len(headers[name])) + " |"])
        markdown.extend("| " + " | ".join(row) + " |" for row in rows)
        markdown.append("")
        print(f"{name.upper()}_MATCH: {len(rows)} rows")
    files["tables.md"] = "\n".join(markdown)
    files["INPUT_SHA256SUMS"] = "".join(digest + "  " + name + "\n"
                                       for name, (_, digest) in INPUTS.items())
    files["validation.txt"] = (
        f"Input source commit: {SOURCE_COMMIT}\n"
        f"Expected manuscript file SHA256: {MANUSCRIPT_SHA256}\n"
        f"Analysis script SHA256: {sha(Path(__file__).read_bytes())}\n"
        f"Python: {sys.version.split()[0]}\n"
        "Six input file checksums matched; 360 CyriPanel iterations checked.\n"
        "16680 CyriPanel records re-scored; 6814 Aldy retained major calls re-scored.\n"
        "Tables 3, 4 and 5 match all 42 manuscript display rows.\n"
        "No-calls remain in denominators; samples have equal weight in Table 3.\n"
        "Table 5 SD uses sample SD (n-1); discordance includes no-calls.\n"
        "Scope: descriptive evaluation-table reconstruction, not caller execution,\n"
        "raw Aldy solution-file verification, depth-matrix recovery or inferential testing.\n"
    )
    files["SHA256SUMS"] = "".join(sha(text.encode()) + "  " + name + "\n"
                                 for name, text in sorted(files.items()))
    args.outdir.mkdir(parents=True, exist_ok=False)
    for name, text in files.items():
        with (args.outdir / name).open("x", encoding="utf-8", newline="") as handle:
            handle.write(text)
    print("TABLES_3_4_5_REPRODUCED:", args.outdir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        sys.exit("STOP: " + str(error))
