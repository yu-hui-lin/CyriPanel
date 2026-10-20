#!/usr/bin/env python3
"""
run_single_iteration.py
========================
Execute ONE benchmark iteration for a given cohort and reference panel size.

Architecture: Each iteration gets its own ISOLATED work directory containing
symlinks to the CyriPanel source tree. This allows parallel SLURM array tasks
to run concurrently without clobbering each other's ref_dir/ or
cnv_panelizer_results/ directories.

Usage:
    python run_single_iteration.py \\
        --cohort cohort1 \\
        --panel-size 5 \\
        --iteration 0 \\
        --seed 42 \\
        --work-dir /path/to/iter_000/work \\
        --out-dir  /path/to/iter_000

Work directory layout created:
    work_dir/
    ├── star_caller.py          → SYMLINK to CyriPanel/star_caller.py
    ├── caller                  → SYMLINK to CyriPanel/caller/
    ├── data                    → SYMLINK to CyriPanel/data/
    ├── depth_calling/          ← real dir (so cnv_panelizer_results is isolated)
    │   ├── *.py                → SYMLINKs to CyriPanel/depth_calling/*.py
    │   ├── run_CNVPanelizer.R  → SYMLINK
    │   └── cnv_panelizer_results/  ← real, empty (written into at runtime)
    └── ref_dir/                ← real dir with symlinks to chosen reference BAMs

When `python work_dir/star_caller.py` runs, its __file__ resolves to
work_dir/star_caller.py, so CyriPanel's hard-coded
`script_dir = os.path.dirname(os.path.abspath(__file__))`
points INTO work_dir, where it finds the iteration-specific ref_dir/.
"""

import os
import sys
import argparse
import csv
import json
import glob
import shutil
import random
import subprocess
import logging
import re
from datetime import datetime

import numpy as np

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# -- Resolve benchmark_config --
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import benchmark_config as cfg


# ──────────────────────────────────────────────────────────────────────
# Argument parsing
# ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Run one CyriPanel benchmark iteration")
    p.add_argument("--cohort", required=True, choices=list(cfg.COHORTS.keys()))
    p.add_argument("--panel-size", type=int, required=True,
                   help="Number of diploid references")
    p.add_argument("--iteration", type=int, required=True,
                   help="Iteration index (0-based)")
    p.add_argument("--seed", type=int, required=True,
                   help="Random seed for this iteration")
    p.add_argument("--work-dir", required=True,
                   help="Isolated work directory (CyriPanel symlinks + ref_dir go here)")
    p.add_argument("--out-dir", required=True,
                   help="Output directory for results (JSON, CSV, manifests)")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────
# CSV readers
# ──────────────────────────────────────────────────────────────────────

def read_csv_map(csv_path, key_col, val_col):
    """Read a two-column CSV into {key: value} dict (for gold standard)."""
    d = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            d[row[key_col].strip()] = row[val_col].strip()
    return d


def read_sample_list(csv_path, id_col="Sample_ID"):
    """
    Read a single-column CSV (or the first column of a multi-column CSV)
    into a list of sample IDs.
    """
    ids = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.append(row[id_col].strip())
    return ids


# ──────────────────────────────────────────────────────────────────────
# BAM discovery
# ──────────────────────────────────────────────────────────────────────

def find_bam(bam_dir, sample_id, suffix=".markdup.bam"):
    """
    Locate BAM and its index for a sample.

    Supports two directory layouts:
      Flat:   bam_dir/{sample_id}{suffix}
      Subdir: bam_dir/{sample_id}/{sample_id}{suffix}   <- Taiwania default
    """
    bam_path = os.path.join(bam_dir, sample_id, sample_id + suffix)
    if not os.path.exists(bam_path):
        bam_path = os.path.join(bam_dir, sample_id + suffix)
    if not os.path.exists(bam_path):
        candidates = glob.glob(
            os.path.join(bam_dir, "**", sample_id + "*" + suffix),
            recursive=True,
        )
        if candidates:
            bam_path = candidates[0]
        else:
            return None, None

    # Index: try .markdup.bai first, then .markdup.bam.bai
    bai_path = bam_path.replace(suffix, suffix[:-1] + "i")
    if not os.path.exists(bai_path):
        bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path):
        bai_path = None

    return bam_path, bai_path


# ──────────────────────────────────────────────────────────────────────
# Isolated work-dir construction (symlink CyriPanel tree)
# ──────────────────────────────────────────────────────────────────────

def _symlink(src, dst):
    """Create a symlink, removing any pre-existing dst."""
    if os.path.islink(dst) or os.path.exists(dst):
        if os.path.islink(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
        else:
            os.remove(dst)
    os.symlink(os.path.abspath(src), dst)


def setup_work_dir(work_dir, cyripanel_src):
    """
    Build an isolated CyriPanel execution environment under work_dir.

    Top-level contents of CyriPanel to symlink:
      - star_caller.py    (file symlink)
      - caller/           (directory symlink – contains star-allele tables)
      - data/             (directory symlink – contains BED, SNP, haplotype tables)

    depth_calling/ is reconstructed as a REAL directory with individual
    symlinks so that cnv_panelizer_results/ inside it can be an iteration-
    specific real directory (prevents concurrent tasks overwriting each
    other's CNV reports).
    """
    os.makedirs(work_dir, exist_ok=True)

    # ---- Top-level star_caller.py, caller/, data/ ----
    for name in ("star_caller.py", "caller", "data"):
        src = os.path.join(cyripanel_src, name)
        if not os.path.exists(src):
            logging.error("CyriPanel source missing: %s", src)
            sys.exit(1)
        _symlink(src, os.path.join(work_dir, name))

    # ---- depth_calling/ : real dir with individual symlinks ----
    dc_src = os.path.join(cyripanel_src, "depth_calling")
    dc_dst = os.path.join(work_dir, "depth_calling")
    os.makedirs(dc_dst, exist_ok=True)

    # Every .py and .R file in depth_calling becomes a symlink
    for fname in os.listdir(dc_src):
        fpath = os.path.join(dc_src, fname)
        if os.path.isfile(fpath) and (fname.endswith(".py") or fname.endswith(".R")):
            _symlink(fpath, os.path.join(dc_dst, fname))

    # Real (empty) cnv_panelizer_results/ — this is where panel_cn.py writes
    cnv_out = os.path.join(dc_dst, "cnv_panelizer_results")
    if os.path.exists(cnv_out):
        shutil.rmtree(cnv_out)
    os.makedirs(cnv_out, exist_ok=True)

    logging.info("  Work dir set up: %s", work_dir)


def setup_ref_dir(work_dir, ref_bam_paths):
    """
    Create work_dir/ref_dir/ and symlink reference BAMs + indices into it.
    """
    ref_dir = os.path.join(work_dir, "ref_dir")
    if os.path.exists(ref_dir):
        shutil.rmtree(ref_dir)
    os.makedirs(ref_dir, exist_ok=True)

    for bam_path, bai_path in ref_bam_paths:
        bam_name = os.path.basename(bam_path)
        _symlink(bam_path, os.path.join(ref_dir, bam_name))

        if bai_path and os.path.exists(bai_path):
            bai_name = os.path.basename(bai_path)
            _symlink(bai_path, os.path.join(ref_dir, bai_name))
            # Also create .bam.bai companion (some tools look for this form)
            companion = os.path.join(ref_dir, bam_name + ".bai")
            if not os.path.exists(companion):
                _symlink(bai_path, companion)

    logging.info("  Symlinked %d reference BAMs into %s",
                 len(ref_bam_paths), ref_dir)
    return ref_dir


# ──────────────────────────────────────────────────────────────────────
# Manifest + CyriPanel execution
# ──────────────────────────────────────────────────────────────────────

def write_manifest(manifest_path, bam_paths):
    """Write one BAM path per line."""
    with open(manifest_path, "w") as f:
        for p in sorted(bam_paths):
            f.write(p + "\n")
    logging.info("  Manifest written: %s (%d samples)",
                 manifest_path, len(bam_paths))


def run_cyripanel(work_dir, manifest_path, out_dir, prefix, genome=38, threads=1):
    """
    Invoke work_dir/star_caller.py (which symlinks to the real one but resolves
    __file__ to work_dir, so it finds the iteration-specific ref_dir).

    Python 3.6 compatible (uses subprocess.PIPE instead of text=True).
    """
    star_caller = os.path.join(work_dir, "star_caller.py")
    cmd = [
        sys.executable, star_caller,
        "--manifest", manifest_path,
        "--genome", str(genome),
        "--outDir", out_dir,
        "--prefix", prefix,
        "--threads", str(threads),
    ]
    logging.info("  CMD: %s", " ".join(cmd))
    t0 = datetime.now()
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=work_dir,
    )
    elapsed = (datetime.now() - t0).total_seconds()
    stdout_str = result.stdout.decode("utf-8", errors="replace") if result.stdout else ""
    stderr_str = result.stderr.decode("utf-8", errors="replace") if result.stderr else ""
    logging.info("  CyriPanel finished in %.1fs  (exit code %d)",
                 elapsed, result.returncode)
    if result.returncode != 0:
        logging.error("  STDERR:\n%s", stderr_str[-2000:])
    return result.returncode, stdout_str, stderr_str, elapsed


def collect_mean_ratios(cnvpanelizer_results_dir, test_sample_ids):
    """
    Parse CNVPanelizer CSV reports for each test sample.
    Returns dict: {sample_id: {region: MeanRatio, ...}}
    """
    mean_ratios = {}
    target_regions = [
        "CYP2D6plusREP6_hapcn2",
        "CYP2D7",
        "D7spacer_hapcn1",
    ]

    for sid in test_sample_ids:
        # Try various filename patterns (e.g. HG00438.markdup_CNV_exon_level_report.csv)
        candidates = glob.glob(
            os.path.join(cnvpanelizer_results_dir, "%s*_report.csv" % sid)
        )
        report = candidates[0] if candidates else None

        if report and os.path.exists(report) and HAS_PANDAS:
            try:
                df = pd.read_csv(report, index_col=0)
                ratios = {}
                for region in target_regions:
                    if region in df.index and "MeanRatio" in df.columns:
                        ratios[region] = float(df.loc[region, "MeanRatio"])
                    else:
                        ratios[region] = None
                mean_ratios[sid] = ratios
            except Exception as e:
                logging.warning("  Failed to parse MeanRatio for %s: %s", sid, e)
                mean_ratios[sid] = {r: None for r in target_regions}
        else:
            mean_ratios[sid] = {r: None for r in target_regions}

    return mean_ratios


# ──────────────────────────────────────────────────────────────────────
# Genotype normalization
# ──────────────────────────────────────────────────────────────────────

def _sort_allele_key(allele):
    """
    Generate a sort key for a single allele or haplotype string.
    Extracts ALL numeric parts as a tuple for natural ordering with tiebreaking.
    """
    nums = re.findall(r'(\d+)', allele)
    return tuple(int(n) for n in nums) if nums else (allele,)


def _expand_multiplication(allele):
    """
    Expand multiplication notation: '*1x2' -> ['*1', '*1'].
    Returns a list of alleles (for joining with '+').
    Non-multiplication alleles return a single-element list: '*4' -> ['*4'].
    """
    m = re.match(r'^(\*\d+(?:\.\d+)?)x(\d+)$', allele.strip())
    if m:
        base, count = m.group(1), int(m.group(2))
        return [base] * count
    return [allele]


def _normalize_haplotype(hap_str):
    """
    Sort sub-alleles within a single haplotype (the '+'-delimited components).
    Expands multiplication notation (*AxN) into repeated '+'-form for comparison
    so that '*1+*1' and '*1x2' compare as equivalent.
    Examples:
      '*36+*10'  -> '*10+*36'
      '*1x2'     -> '*1+*1'
      '*4+*1x2'  -> '*1+*1+*4'
    """
    # Split by '+' if present, then expand x-notation in each component
    if "+" in hap_str:
        sub_alleles = hap_str.split("+")
    else:
        sub_alleles = [hap_str]
    expanded = []
    for s in sub_alleles:
        expanded.extend(_expand_multiplication(s))
    if len(expanded) == 1:
        return expanded[0]
    return "+".join(sorted(expanded, key=_sort_allele_key))


def normalize_genotype(gt_str):
    """
    Normalize a genotype string for comparison.
      1. Sort sub-alleles within each haplotype:  *36+*10  -> *10+*36
      2. Sort the two haplotypes:                 *10+*36/*1 -> *1/*10+*36
      3. Strip whitespace
    """
    if gt_str is None or gt_str in ("None", "", "Error"):
        return None

    gt_str = gt_str.strip()

    if "/" in gt_str:
        parts = gt_str.split("/")
        parts = [_normalize_haplotype(p) for p in parts]
        parts_sorted = sorted(parts, key=_sort_allele_key)
        return "/".join(parts_sorted)

    return _normalize_haplotype(gt_str)


def compare_genotype(called, gold):
    """
    Compare called genotype to gold standard.
    Returns: 'concordant', 'discordant', 'no_call', or 'no_gold'
    """
    norm_called = normalize_genotype(called)
    norm_gold = normalize_genotype(gold)

    if norm_called is None:
        return "no_call"
    if norm_gold is None:
        return "no_gold"
    if norm_called == norm_gold:
        return "concordant"
    return "discordant"


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [iter-%(process)d] %(message)s",
    )

    cohort_cfg = cfg.COHORTS[args.cohort]
    panel_size = args.panel_size
    iteration = args.iteration
    out_dir = args.out_dir
    work_dir = args.work_dir

    os.makedirs(out_dir, exist_ok=True)

    logging.info("=== Iteration %d | %s | panel_size=%d | seed=%d ===",
                 iteration, args.cohort, panel_size, args.seed)
    logging.info("  out-dir:  %s", out_dir)
    logging.info("  work-dir: %s", work_dir)

    # -- 1. Load sample lists --
    diploid_ids = sorted(read_sample_list(
        cohort_cfg["diploid_csv"], cfg.COL_SAMPLE_ID
    ))
    logging.info("  Loaded %d known diploid samples", len(diploid_ids))

    gold_map = read_csv_map(
        cohort_cfg["gold_csv"], cfg.COL_SAMPLE_ID, cfg.COL_GENOTYPE
    )
    logging.info("  Loaded %d gold-standard genotypes", len(gold_map))

    # -- 2. Discover BAMs --
    bam_dir = cohort_cfg["bam_dir"]
    suffix = cohort_cfg["bam_suffix"]

    all_bams = glob.glob(os.path.join(bam_dir, "*", "*" + suffix))
    if not all_bams:
        all_bams = glob.glob(os.path.join(bam_dir, "*" + suffix))

    all_sample_ids = {
        os.path.basename(b).replace(suffix, ""): b for b in all_bams
    }
    logging.info("  Found %d BAM files in %s", len(all_sample_ids), bam_dir)

    # Only diploids that have corresponding BAMs are eligible as references
    valid_diploid_ids = [s for s in diploid_ids if s in all_sample_ids]
    if len(valid_diploid_ids) < panel_size:
        logging.error(
            "  Only %d diploid BAMs found, need %d. Aborting.",
            len(valid_diploid_ids), panel_size,
        )
        sys.exit(1)

    # -- 3. Random reference panel selection --
    rng = random.Random(args.seed)
    ref_ids = sorted(rng.sample(valid_diploid_ids, panel_size))
    test_ids = sorted([s for s in all_sample_ids if s not in ref_ids])

    logging.info("  Selected %d references: %s%s",
                 len(ref_ids), ref_ids[:5], "..." if len(ref_ids) > 5 else "")
    logging.info("  Test set: %d samples", len(test_ids))

    # Save reference selection for reproducibility
    with open(os.path.join(out_dir, "reference_samples.json"), "w") as f:
        json.dump({"seed": args.seed, "references": ref_ids}, f, indent=2)

    # -- 4. Build isolated work dir (symlink CyriPanel tree) --
    setup_work_dir(work_dir, cfg.CYRIPANEL_SRC_DIR)

    # -- 5. Populate ref_dir/ inside the work dir --
    ref_bam_paths = []
    for sid in ref_ids:
        bam_path = all_sample_ids[sid]
        _, bai_path = find_bam(bam_dir, sid, suffix)
        if bai_path is None:
            candidate = bam_path + ".bai"
            if os.path.exists(candidate):
                bai_path = candidate
        ref_bam_paths.append((bam_path, bai_path))

    setup_ref_dir(work_dir, ref_bam_paths)

    # -- 6. Write test manifest --
    manifest_path = os.path.join(out_dir, "test_manifest.txt")
    test_bam_paths = [all_sample_ids[s] for s in test_ids]
    write_manifest(manifest_path, test_bam_paths)

    # -- 7. Run CyriPanel via the isolated work dir --
    prefix = "%s_n%d_iter%03d" % (args.cohort, panel_size, iteration)
    threads = getattr(cfg, "SLURM_CPUS_PER_TASK", 1)
    rc, stdout, stderr, elapsed = run_cyripanel(
        work_dir, manifest_path, out_dir, prefix, threads=threads
    )

    # -- 8. Parse CyriPanel JSON output --
    json_path = os.path.join(out_dir, "%s.json" % prefix)
    cyripanel_results = {}
    if os.path.exists(json_path):
        with open(json_path) as f:
            raw_results = json.load(f)
        # Remap keys from e.g. "HG00438.markdup" -> "HG00438"
        stem_to_strip = suffix.replace(".bam", "").replace(".cram", "")
        for key, val in raw_results.items():
            clean_key = key.replace(stem_to_strip, "") if stem_to_strip else key
            cyripanel_results[clean_key] = val
        logging.info("  Parsed %d results from JSON (stripped '%s' from keys)",
                     len(cyripanel_results), stem_to_strip)
    else:
        logging.error("  CyriPanel JSON not found: %s", json_path)

    # -- 9. Collect CNVPanelizer MeanRatios from iteration-specific dir --
    cnvp_dir = os.path.join(work_dir, "depth_calling", "cnv_panelizer_results")
    mean_ratios = collect_mean_ratios(cnvp_dir, test_ids)

    # -- 10. Compare to gold standard --
    per_sample = {}
    n_concordant = n_discordant = n_no_call = n_pass = n_tested = 0

    for sid in test_ids:
        gold_gt = gold_map.get(sid)
        result = cyripanel_results.get(sid, {})
        called_gt = result.get("Genotype")
        filt = result.get("Filter")
        total_cn = result.get("Total_CN")
        spacer_cn = result.get("Spacer_CN")
        cnv_group = result.get("CNV_group")

        comparison = compare_genotype(
            str(called_gt) if called_gt else None, gold_gt
        )

        if gold_gt is not None:
            n_tested += 1
            if comparison == "concordant":
                n_concordant += 1
            elif comparison == "discordant":
                n_discordant += 1
            elif comparison == "no_call":
                n_no_call += 1

        if filt == "PASS":
            n_pass += 1

        per_sample[sid] = {
            "gold_genotype": gold_gt,
            "called_genotype": str(called_gt) if called_gt else None,
            "filter": str(filt) if filt else None,
            "comparison": comparison,
            "total_cn": total_cn,
            "spacer_cn": spacer_cn,
            "cnv_group": cnv_group,
            "mean_ratios": mean_ratios.get(sid, {}),
        }

    # -- 11. Compute summary statistics --
    concordance = n_concordant / n_tested if n_tested > 0 else None
    pass_rate = n_pass / len(test_ids) if test_ids else None
    no_call_rate = n_no_call / n_tested if n_tested > 0 else None

    mr_stats = {}
    for region in ["CYP2D6plusREP6_hapcn2", "CYP2D7", "D7spacer_hapcn1"]:
        values = [
            per_sample[s]["mean_ratios"].get(region)
            for s in test_ids
            if per_sample[s]["mean_ratios"].get(region) is not None
        ]
        if values:
            arr = np.array(values)
            mr_stats[region] = {
                "mean":   float(np.mean(arr)),
                "std":    float(np.std(arr)),
                "median": float(np.median(arr)),
                "min":    float(np.min(arr)),
                "max":    float(np.max(arr)),
                "n":      len(values),
            }
        else:
            mr_stats[region] = {"mean": None, "std": None, "n": 0}

    summary = {
        "cohort": args.cohort,
        "panel_size": panel_size,
        "iteration": iteration,
        "seed": args.seed,
        "n_references": len(ref_ids),
        "n_test_samples": len(test_ids),
        "n_tested_with_gold": n_tested,
        "n_concordant": n_concordant,
        "n_discordant": n_discordant,
        "n_no_call": n_no_call,
        "concordance": concordance,
        "pass_rate": pass_rate,
        "no_call_rate": no_call_rate,
        "mean_ratio_stats": mr_stats,
        "elapsed_seconds": elapsed,
        "cyripanel_exit_code": rc,
        "reference_ids": ref_ids,
    }

    # -- 12. Write iteration_results.json --
    results_out = {"summary": summary, "per_sample": per_sample}
    results_path = os.path.join(out_dir, "iteration_results.json")
    with open(results_path, "w") as f:
        json.dump(results_out, f, indent=2)

    # -- 13. Write per-sample CSV (publication raw data) --
    csv_columns = [
        "Sample_ID",
        "Called_Genotype",
        "Gold_Genotype",
        "Filter",
        "Comparison",
        "Total_CN",
        "Spacer_CN",
        "CNV_Group",
        "MeanRatio_CYP2D6plusREP6",
        "MeanRatio_CYP2D7",
        "MeanRatio_D7spacer",
    ]
    csv_path = os.path.join(out_dir, "iteration_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns)
        writer.writeheader()
        for sid in test_ids:
            ps = per_sample[sid]
            mr = ps.get("mean_ratios", {})
            writer.writerow({
                "Sample_ID": sid,
                "Called_Genotype": ps["called_genotype"] or "",
                "Gold_Genotype": ps["gold_genotype"] or "",
                "Filter": ps["filter"] or "",
                "Comparison": ps["comparison"],
                "Total_CN": ps["total_cn"] if ps["total_cn"] is not None else "",
                "Spacer_CN": ps["spacer_cn"] if ps["spacer_cn"] is not None else "",
                "CNV_Group": ps["cnv_group"] or "",
                "MeanRatio_CYP2D6plusREP6": mr.get("CYP2D6plusREP6_hapcn2", ""),
                "MeanRatio_CYP2D7": mr.get("CYP2D7", ""),
                "MeanRatio_D7spacer": mr.get("D7spacer_hapcn1", ""),
            })
    logging.info("  CSV -> %s  (%d rows)", csv_path, len(test_ids))

    # -- 14. Summary log --
    logging.info("  -- SUMMARY --")
    if concordance is not None:
        logging.info("  Concordance: %d/%d = %.3f",
                     n_concordant, n_tested, concordance)
    if pass_rate is not None:
        logging.info("  PASS rate:   %d/%d = %.3f",
                     n_pass, len(test_ids), pass_rate)
    logging.info("  No-call:     %d/%d", n_no_call, n_tested)
    logging.info("  Discordant:  %d/%d", n_discordant, n_tested)
    logging.info("  Results -> %s", results_path)

    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
