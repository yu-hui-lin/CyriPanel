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


import os
import subprocess
import pandas as pd
import math
import logging


def get_cn_from_cnvpanelizer(bam_file, r_script_path, output_dir, bed_file_path, reference_dir_path):
    """
    Runs the CNVPanelizer R script and calculates total_cn and d7_spacer.

    Args:
        bam_file (str): Path to the sample BAM file.
        r_script_path (str): Path to the run_CNVPanelizer.R script.
        output_dir (str): Directory where the CNVPanelizer report will be saved.
        bed_file_path (str): Path to the BED file for CNVPanelizer.
        reference_dir_path (str): Path to the directory of reference BAMs.

    Returns:
        tuple: A tuple containing (total_cn, d7_spacer), or (None, None) on failure.
               Returns (None, None) if the process fails.
    """
    # --- Validate input files and create output directory. ---
    if not os.path.exists(bam_file):
        logging.error(f"BAM file not found: {bam_file}")
        return None, None

    if not os.path.exists(r_script_path):
        logging.error(f"R script not found: {r_script_path}")
        return None, None

    os.makedirs(output_dir, exist_ok=True)

    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
    report_filename = os.path.join(output_dir, f"{bam_basename}_CNV_exon_level_report.csv")

    # --- STEP 1: Execute CNVPanelizer R script via subprocess ---
    # Calls external R script with sample BAM, BED file, reference BAMs directory, and output directory.
    # CNVPanelizer normalizes sample depth against reference panel to calculate copy number ratios. 
    logging.info(f"Running CNVPanelizer for {bam_basename}...")
    try:
        r_command = [
            "Rscript", 
            r_script_path, 
            bam_file, 
            bed_file_path, 
            reference_dir_path, 
            output_dir
        ]
        subprocess.run(r_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        logging.info("CNVPanelizer R script executed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing CNVPanelizer R script for {bam_basename}:")
        logging.error(e.stderr)
        return None, None

    # --- STEP 2: Load CNVPanelizer CSV output report ---
    # Report contains MeanRatio values for each genomic region (gene-level CN estimates). 
    if not os.path.exists(report_filename):
        logging.error(f"CNVPanelizer report not found: {report_filename}")
        return None, None

    try:
        # First column (gene names) used as index
        report_df = pd.read_csv(report_filename, index_col=0)
    except Exception as e:
        logging.error(f"Failed to read or parse the CNVPanelizer report: {e}")
        return None, None

    # --- STEP 3: Extract MeanRatio values for key genomic regions. ---
    # Three regions extracted: CYP2D6+REP6, CYP2D7, and D7 spacer (each has expected haploid CN). 
    try:
        mean_ratio_2d6 = report_df.loc["CYP2D6plusREP6_hapcn2", "MeanRatio"]
        mean_ratio_2d7 = report_df.loc["CYP2D7", "MeanRatio"]
        mean_ratio_d7_spacer = report_df.loc["D7spacer_hapcn1", "MeanRatio"]
        logging.info("Successfully extracted MeanRatio values from the report.")
    except KeyError as e:
        logging.error(f"Could not find required gene name in the report's index: {e}")
        return None, None

    # --- STEP 4: Calculate total_cn and d7_spacer from MeanRatio values. ---
    # Formula: total_cn = round(2D6_ratio*2) + round(2D7_ratio*2)
    # Combines CYP2D6 and CYP2D7 CN estimates to get total gene cluster copy number.
    # d7_spacer = round(spacer_ratio*2) gives spacer region copy number.
    # Epsilon values (currently 0) allow tuning/hysteresis based on reference panel characteristics. 
    cyp2d6_eps = 0  # Hysteresis window: tune on reference samples if needed
    cyp2d7_eps = 0  # Hysteresis window: tune on reference samples if needed
    total_cn = int(round((mean_ratio_2d6 + cyp2d6_eps) * 2) + round((mean_ratio_2d7 + cyp2d7_eps) * 2))

    d7_spacer_eps = 0  # Hysteresis window: tune on reference samples if needed
    d7_spacer = int(round((mean_ratio_d7_spacer + d7_spacer_eps) * 2))

    logging.info(f"Calculated total_cn: {total_cn}")
    logging.info(f"Calculated d7_spacer: {d7_spacer}")

    return total_cn, d7_spacer
