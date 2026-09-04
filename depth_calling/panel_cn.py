#!/usr/bin/env python3
#
# CyriPanel: CYP2D6 genotyper for targeted sequencing panels
# Modified and integrated by Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# Original Cyrius Copyright (c) 2019-2020 Illumina, Inc.
# Original Author: Xiao Chen <xchen2@illumina.com>
# BCyrius (updated CYP2D6 star alleles) Copyright (c) 2024 Andreas Halman
#
# Modifications include:
# - Integration of CNVPanelizer to override Gaussian Mixture Model (GMM) CNV calculations.
# - Optimization specifically tailored for targeted sequencing panels.
# - Added fallback mechanisms to improve robustness against the noisier read depths of targeted sequencing data.
# - Updated CYP2D6 star allele definitions based on BCyrius.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


# --- NEW IN CYRIPANEL: CNVPanelizer integration module ---
# Replaces Cyrius GMM-based CNV calling with CNVPanelizer's reference-panel normalization approach.
# Executes external R script and parses output to determine total_cn and d7_spacer for targeted panels.
#
# --- MODIFIED: Added pysam-based fallback CN estimator ---
# When CNVPanelizer produces implausibly low MeanRatios (typically due to BWA-MEM assigning
# low MAPQ to reads in the CYP2D6/CYP2D7 paralogous region), a fallback estimator counts
# all passing reads (MQ>=0) across the three target regions and computes CN ratios directly.
# This addresses the systematic under-estimation seen with BWA-MEM aligned data where
# >80-90% of reads receive MQ=0 in the D6/D7 region.


import os
import subprocess
import glob
import math
import logging

import pandas as pd
import pysam


# ──────────────────────────────────────────────────────────────────────
# BED-defined target regions for CN estimation (GRCh38 coordinates)
# ──────────────────────────────────────────────────────────────────────
_CN_REGIONS = {
    "CYP2D6plusREP6_hapcn2": ("chr22", 42123192, 42132032),
    "D7spacer_hapcn1":       ("chr22", 42138124, 42139676),
    "CYP2D7":                ("chr22", 42139676, 42145745),
}

# MeanRatio threshold: if both CYP2D6 and CYP2D7 MeanRatios from
# CNVPanelizer fall below this value, trigger the pysam fallback.
_LOW_RATIO_THRESHOLD = 0.3


# ──────────────────────────────────────────────────────────────────────
# Pysam-based read counting (fallback CN estimator)
# ──────────────────────────────────────────────────────────────────────

def _count_region_reads(bam_path, chrom, start, end, min_mapq=0):
    """
    Count passing reads in a genomic region.

    Passing reads: mapped, non-duplicate, non-secondary, non-supplementary,
    with mapping quality >= min_mapq.

    Returns:
        int: Number of passing reads.
    """
    count = 0
    try:
        bf = pysam.AlignmentFile(bam_path, "rb")
        for read in bf.fetch(chrom, start, end):
            if (read.is_unmapped or read.is_duplicate
                    or read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality >= min_mapq:
                count += 1
        bf.close()
    except Exception as e:
        logging.warning("Failed to count reads in %s:%d-%d for %s: %s",
                        chrom, start, end, bam_path, e)
    return count


def _get_reference_mean_counts(reference_dir_path, min_mapq=0):
    """
    Compute mean read counts across all reference BAMs for each target region.

    Args:
        reference_dir_path: Directory containing reference BAM files.
        min_mapq: Minimum mapping quality for read counting.

    Returns:
        dict: {region_name: mean_count} for each region in _CN_REGIONS.
    """
    ref_bams = sorted(glob.glob(os.path.join(reference_dir_path, "*.bam")))
    if not ref_bams:
        logging.error("No reference BAMs found in %s for fallback CN estimation.",
                      reference_dir_path)
        return None

    region_totals = {name: 0.0 for name in _CN_REGIONS}
    n_refs = len(ref_bams)

    for bam_path in ref_bams:
        for region_name, (chrom, start, end) in _CN_REGIONS.items():
            count = _count_region_reads(bam_path, chrom, start, end, min_mapq)
            region_totals[region_name] += count

    ref_means = {name: total / n_refs for name, total in region_totals.items()}
    logging.info("  Fallback CN: Reference mean counts from %d BAMs (MQ>=%d):",
                 n_refs, min_mapq)
    for name, mean_val in ref_means.items():
        logging.info("    %s: %.1f", name, mean_val)

    return ref_means


def _estimate_cn_from_reads(bam_path, reference_dir_path, min_mapq=0):
    """
    Estimate total_cn and d7_spacer from raw read counts (pysam fallback).

    Uses the same formula as CNVPanelizer:
        MeanRatio = sample_count / reference_mean_count
        total_cn  = round(D6_ratio * 2) + round(D7_ratio * 2)
        d7_spacer = round(spacer_ratio * 2)

    Args:
        bam_path: Path to the test sample BAM file.
        reference_dir_path: Directory containing reference BAM files.
        min_mapq: Minimum mapping quality for counting. Default 0 to
                  capture BWA-MEM low-MAPQ reads in paralogous regions.

    Returns:
        tuple: (total_cn, d7_spacer, ratios_dict) or (None, None, None) on failure.
               ratios_dict maps region name to MeanRatio for logging.
    """
    ref_means = _get_reference_mean_counts(reference_dir_path, min_mapq)
    if ref_means is None:
        return None, None, None

    ratios = {}
    for region_name, (chrom, start, end) in _CN_REGIONS.items():
        sample_count = _count_region_reads(bam_path, chrom, start, end, min_mapq)
        mean_ref = ref_means[region_name]
        if mean_ref > 0:
            ratios[region_name] = sample_count / mean_ref
        else:
            logging.warning("  Fallback CN: Reference mean is 0 for %s", region_name)
            ratios[region_name] = 0.0

    d6_ratio = ratios.get("CYP2D6plusREP6_hapcn2", 0.0)
    d7_ratio = ratios.get("CYP2D7", 0.0)
    spacer_ratio = ratios.get("D7spacer_hapcn1", 0.0)

    total_cn = int(round(d6_ratio * 2) + round(d7_ratio * 2))
    d7_spacer = int(round(spacer_ratio * 2))

    logging.info("  Fallback CN ratios (MQ>=%d): D6=%.3f  D7=%.3f  Spacer=%.3f",
                 min_mapq, d6_ratio, d7_ratio, spacer_ratio)
    logging.info("  Fallback CN result: total_cn=%d, d7_spacer=%d", total_cn, d7_spacer)

    return total_cn, d7_spacer, ratios


# ──────────────────────────────────────────────────────────────────────
# Main entry point: CNVPanelizer with pysam fallback
# ──────────────────────────────────────────────────────────────────────

def get_cn_from_cnvpanelizer(bam_file, r_script_path, output_dir,
                              bed_file_path, reference_dir_path):
    """
    Runs the CNVPanelizer R script and calculates total_cn and d7_spacer.

    If CNVPanelizer produces implausibly low MeanRatios (both CYP2D6 and
    CYP2D7 below threshold), falls back to a pysam-based read count
    estimator that includes low-MAPQ reads. This addresses BWA-MEM's
    tendency to assign MQ=0 to reads in the CYP2D6/CYP2D7 paralogous
    region, which causes CNVPanelizer (at MQ>=20) to severely
    underestimate copy numbers.

    Args:
        bam_file (str): Path to the sample BAM file.
        r_script_path (str): Path to the run_CNVPanelizer.R script.
        output_dir (str): Directory where the CNVPanelizer report will be saved.
        bed_file_path (str): Path to the BED file for CNVPanelizer.
        reference_dir_path (str): Path to the directory of reference BAMs.

    Returns:
        tuple: A tuple containing (total_cn, d7_spacer), or (None, None) on failure.
    """
    # --- Validate input files and create output directory. ---
    if not os.path.exists(bam_file):
        logging.error("BAM file not found: %s", bam_file)
        return None, None

    if not os.path.exists(r_script_path):
        logging.error("R script not found: %s", r_script_path)
        return None, None

    os.makedirs(output_dir, exist_ok=True)

    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
    report_filename = os.path.join(
        output_dir, "{}_CNV_exon_level_report.csv".format(bam_basename)
    )

    # --- STEP 1: Execute CNVPanelizer R script via subprocess ---
    logging.info("Running CNVPanelizer for %s...", bam_basename)
    try:
        r_command = [
            "Rscript",
            r_script_path,
            bam_file,
            bed_file_path,
            reference_dir_path,
            output_dir
        ]
        subprocess.run(r_command, check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
        logging.info("CNVPanelizer R script executed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error("Error executing CNVPanelizer R script for %s:", bam_basename)
        logging.error(e.stderr)
        # CNVPanelizer failed entirely — try pysam fallback directly
        logging.warning("CNVPanelizer failed for %s. Attempting pysam fallback.",
                        bam_basename)
        result = _estimate_cn_from_reads(bam_file, reference_dir_path, min_mapq=0)
        if result[0] is not None:
            return result[0], result[1]
        return None, None

    # --- STEP 2: Load CNVPanelizer CSV output report ---
    if not os.path.exists(report_filename):
        logging.error("CNVPanelizer report not found: %s", report_filename)
        return None, None

    try:
        report_df = pd.read_csv(report_filename, index_col=0)
    except Exception as e:
        logging.error("Failed to read or parse the CNVPanelizer report: %s", e)
        return None, None

    # --- STEP 3: Extract MeanRatio values for key genomic regions. ---
    try:
        mean_ratio_2d6 = report_df.loc["CYP2D6plusREP6_hapcn2", "MeanRatio"]
        mean_ratio_2d7 = report_df.loc["CYP2D7", "MeanRatio"]
        mean_ratio_d7_spacer = report_df.loc["D7spacer_hapcn1", "MeanRatio"]
        logging.info("CNVPanelizer MeanRatios: D6=%.3f  D7=%.3f  Spacer=%.3f",
                     mean_ratio_2d6, mean_ratio_2d7, mean_ratio_d7_spacer)
    except KeyError as e:
        logging.error("Could not find required gene name in the report's index: %s", e)
        return None, None

    # --- STEP 3b: Check for implausibly low MeanRatios (BWA-MEM MAPQ issue) ---
    # When BWA-MEM assigns MQ=0 to most reads in the paralogous region,
    # CNVPanelizer (which filters at MQ>=20) sees near-zero depth and produces
    # MeanRatio ≈ 0. We detect this and fall back to pysam counting at MQ>=0.
    needs_fallback = False

    if mean_ratio_2d6 < _LOW_RATIO_THRESHOLD and mean_ratio_2d7 < _LOW_RATIO_THRESHOLD:
        logging.warning(
            "Both CYP2D6 (%.3f) and CYP2D7 (%.3f) MeanRatios below %.2f — "
            "likely BWA-MEM low-MAPQ issue. Triggering pysam fallback.",
            mean_ratio_2d6, mean_ratio_2d7, _LOW_RATIO_THRESHOLD
        )
        needs_fallback = True

    if needs_fallback:
        result = _estimate_cn_from_reads(bam_file, reference_dir_path, min_mapq=0)
        if result[0] is not None:
            return result[0], result[1]
        else:
            logging.error("Pysam fallback also failed for %s.", bam_basename)
            # Fall through to use CNVPanelizer values as last resort

    # --- STEP 4: Calculate total_cn and d7_spacer from MeanRatio values. ---
    cyp2d6_eps = 0
    cyp2d7_eps = 0
    total_cn = int(round((mean_ratio_2d6 + cyp2d6_eps) * 2)
                   + round((mean_ratio_2d7 + cyp2d7_eps) * 2))

    d7_spacer_eps = 0
    d7_spacer = int(round((mean_ratio_d7_spacer + d7_spacer_eps) * 2))

    logging.info("Calculated total_cn: %d", total_cn)
    logging.info("Calculated d7_spacer: %d", d7_spacer)

    return total_cn, d7_spacer
