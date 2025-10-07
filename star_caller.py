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


import os
import sys
import argparse
import json
import logging
import datetime
from collections import namedtuple, OrderedDict
import pysam

# --- MODIFIED: Import glob module for automated BED file discovery in data directory. ---
import glob
# --- END OF MODIFICATION ---

from caller.call_variants import FALLBACK_CN_COMBINATIONS

from depth_calling.snp_count import (
    get_supporting_reads,
    get_supporting_reads_single_region,
    get_fraction,
    get_snp_position,
)

# --- MODIFIED: Remove GMM-related imports as CNVPanelizer replaces Gaussian Mixture Model for CN determination. ---
# from depth_calling.gmm import Gmm
# --- END OF MODIFICATION ---

from depth_calling.utilities import (
    # --- MODIFIED: Remove parse_gmm_file import as GMM parameters are no longer needed. ---
    # parse_gmm_file,
    # --- END OF MODIFICATION ---
    parse_region_file,
    open_alignment_file,
)

# --- MODIFIED: Import CNVPanelizer integration function to obtain total_cn and d7_spacer from panel data. ---
from depth_calling.panel_cn import get_cn_from_cnvpanelizer
# --- END OF MODIFICATION ---

# --- MODIFIED: Remove bin_count module imports as read depth normalization is handled by CNVPanelizer. ---
# from depth_calling.bin_count import (
#     get_normed_depth,
#     get_normed_depth_from_count,
#     get_read_length,
# )
# --- END OF MODIFICATION ---

from caller.call_variants import (
    NOISY_VAR,
    call_cn_snp,
    call_cn_var,
    call_cn_var_homo,
    get_allele_counts_var42128936,
    update_var42128936,
    get_called_variants,
    call_exon9gc,
    call_var42126938,
    call_var42127526_var42127556,
    call_var42127803hap,
)
from caller.cnv_hybrid import get_cnvtag
from caller.construct_star_table import get_hap_table
from caller.match_star_allele import match_star

MAD_THRESHOLD = 0.11
EXON9_SITE1 = 7
EXON9_SITE2 = 8
# --- MODIFIED: Remove HIGH_CN_DEPTH_THRESHOLD as GMM-based confidence scoring is not used in CyriPanel. ---
# HIGH_CN_DEPTH_THRESHOLD = 7.5 
# --- END OF MODIFICATION ---
HAPLOTYPE_VAR = ["g.42126938C>T", "g.42127803C>T", "g.42127526C>T_g.42127556T>C"]

# --- MODIFIED: Remove gmm_parameter from resource_info namedtuple as GMM is not used in CyriPanel. ---
resource_info = namedtuple(
    "resource_info",
    "genome region_dic snp_db var_db var_homo_db haplotype_db var_list star_combinations",
)
# --- END OF MODIFICATION ---

exon9_values = namedtuple(
    "exon9_values", "exon9_cn exon9cn_in_consensus exon9_raw_site1 exon9_raw_site2"
)
# Below are the SV configurations that the caller is able to call
CNV_ACCEPTED = [
    "star5_star5",
    "star13_star13",
    "star13intron1_star13intron1",
    "star5",
    "star13",
    "star13intron1",
    "star5_star5_star68",
    "star5_star68",
    "cn2",
    "exon9hyb_star5",
    "dup_star13",
    "dup_star13intron1",
    "star13_star68",
    "cn3",
    "exon9hyb",
    "star68",
    "cn4",
    "exon9hyb_exon9hyb",
    "star68_star68",
    "dup_exon9hyb",
    "dup_star68",
    "exon9hyb_star68",
    "cn5",
    "exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68",
    "cn6",
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68_star68",
]


# --- MODIFIED: Add automated BED file discovery function to support flexible panel configurations.
# Users can provide their own panel BED file in the data directory without hardcoding filenames. ---
def find_bed_file(data_dir):
    """
    Find and validate .bed files in the data directory.
    
    Returns:
        str: Path to the .bed file to use
    
    Raises:
        Exception: If no .bed files found or validation fails
    """
    # Search for all .bed files in the data directory
    bed_pattern = os.path.join(data_dir, "*.bed")
    bed_files = glob.glob(bed_pattern)
    
    if not bed_files:
        raise Exception(f"No .bed files found in the data directory: {data_dir}")
    
    if len(bed_files) > 1:
        bed_filenames = [os.path.basename(f) for f in bed_files]
        logging.warning(
            f"Multiple .bed files found in data directory: {', '.join(bed_filenames)}. "
            f"Using the first one: {bed_filenames[0]}"
        )
    
    selected_bed_file = bed_files[0]
    logging.info(f"Using BED file for CNVPanelizer: {selected_bed_file}")
    
    # Check if the selected file exists and is readable
    if not os.path.exists(selected_bed_file):
        raise Exception(f"Selected BED file does not exist: {selected_bed_file}")
    
    if not os.access(selected_bed_file, os.R_OK):
        raise Exception(f"Selected BED file is not readable: {selected_bed_file}")
    
    return selected_bed_file
# --- END OF MODIFICATION ---


def load_parameters():
    """Return parameters."""
    parser = argparse.ArgumentParser(
        description="Call CYP2D6 genotypes from a targeted sequencing BAM file. Diploid references are required."
    )
    parser.add_argument(
        "-m",
        "--manifest",
        help="Manifest listing absolute paths to BAM/CRAM files",
        required=True,
    )
    
    # --- MODIFIED: Restrict to hg38 only as CyriPanel is validated for this reference build. Original Cyrius supported hg19/37/38. ---
    parser.add_argument('--genome', type=int, choices=[38], help='Genome build version')
    # --- END OF MODIFICATION ---
    
    parser.add_argument("-o", "--outDir", help="Output directory", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix to output file", required=True)
    parser.add_argument(
        "-t",
        "--threads",
        help="Optional, number of threads to use. Default is 1",
        type=int,
        required=False,
        default=1,
    )
    parser.add_argument(
        "--countFilePath",
        help="Optional path to count files (not used in targeted sequencing mode)",
        required=False
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Optional path to reference fasta file for CRAM",
        required=False,
    )

    args = parser.parse_args()
    if args.genome != 38:
        raise Exception("Genome not recognized. Only hg38 (genome=38) is supported.")

    return args


# --- MODIFIED: Add helper function to format count pairs for consistent output formatting.
# Handles D6/D7 allele count pairs with proper ordering and None value handling. ---
def _fmt_pair(pair, order="d7,d6", none_token="NA"):
    """
    Format a 2-element count pair (typically [d6_alt, d7_alt]) as 'x,y'.
    Emits none_token for missing (None) values.
    order='d7,d6' keeps legacy output convention.
    """
    if not pair or len(pair) != 2:
        return f"{none_token},{none_token}"
    d6, d7 = pair[0], pair[1]
    def tok(v): return str(int(v)) if v is not None else none_token
    return f"{tok(d7)},{tok(d6)}" if order == "d7,d6" else f"{tok(d6)},{tok(d7)}"
# --- END OF MODIFICATION ---


# --- MODIFIED: Add total_cn and d7_spacer_param parameters to function signature to receive CNVPanelizer outputs.
# These replace GMM-derived CN calls and enable direct integration of panel-based CN determination. ---
def d6_star_caller(
    bam, call_parameters, threads, total_cn, d7_spacer_param,
    count_file=None, reference_fasta=None, index_name=None
):
# --- END OF MODIFICATION ---
    """Return CYP2D6 star allele diplotype calls for each sample."""
    d6_call = namedtuple(
        "d6_call",
        "Coverage_MAD Median_depth Total_CN Spacer_CN Total_CN_raw \
        Spacer_CN_raw Variants_called CNV_group Genotype Filter Raw_star_allele \
        Call_info Exon9_CN CNV_consensus d67_snp_call d67_snp_raw \
        Variant_raw_count",
    )

    # --- MODIFIED: Replace GMM-based depth normalization and CN calling with predetermined CNVPanelizer values.
    # Original Cyrius Steps 1-2: (1) Read counting and normalization, (2) GMM and CN call.
    # CyriPanel: Bypass these steps entirely and use total_cn/d7_spacer from CNVPanelizer.
    # Set placeholder values (0.0) for Coverage_MAD and Median_depth as these metrics are GMM-specific.
    # The raw depth values (d67_depth_placeholder, spacer_depth_placeholder) are also initialized as placeholders. ---
    coverage_mad_value = 0.0
    median_depth_value = 0.0
    d67_depth_placeholder = float(total_cn) if total_cn is not None else None
    spacer_depth_placeholder = 0.0

    cn_call_tuple_def = namedtuple("cn_call", "d67_cn d67_depth spacer_cn spacer_depth")

    # Handle the d7_spacer value. If a specific value is provided via d7_spacer_param, use it. Otherwise, default to 4.
    spacer_cn_value = d7_spacer_param if d7_spacer_param is not None else 4
    if d7_spacer_param is not None:
         logging.info(f"Using CNVPanelizer input of d7_spacer_cn: {d7_spacer_param}")
    else:
         logging.info("Using default d7_spacer_cn: 4")

    # Directly construct the raw_cn_call tuple using the predetermined total_cn and spacer_cn_value, bypassing the GMM calculation.
    raw_cn_call = cn_call_tuple_def(
        d67_cn=total_cn,
        d67_depth=d67_depth_placeholder,
        spacer_cn=spacer_cn_value,
        spacer_depth=spacer_depth_placeholder
    )

    high_cn_low_confidence = False # This remains from the original but is less relevant without GMM.
    bamfile = open_alignment_file(bam, reference_fasta, index_filename=index_name)

    # Add a check for invalid total_cn input. If the provided total_cn is None or negative, the function returns a no-call with an "Invalid_Total_CN" filter.
    if raw_cn_call.d67_cn is None or raw_cn_call.d67_cn < 0:
        logging.warning(f"Total CN ({total_cn}) is invalid. Reporting no-call.")
        sample_call = d6_call(
            coverage_mad_value,
            median_depth_value,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None, None, None, "Invalid_Total_CN", None, None, None, None, None, None, OrderedDict()
        )
        if bamfile: bamfile.close()
        return sample_call       
    # --- END OF MODIFICATION ---

    # --- REST OF THE ORIGINAL CYRIUS LOGIC ---
    # 3. Get allele counts at D6/D7 SNP (base difference) sites and target variant sites
    # D6/D7 base difference sites. Get read counts at both D6/D7 positions.
    snp_db = call_parameters.snp_db
    snp_d6, snp_d7 = get_supporting_reads(
        bamfile, snp_db.dsnp1, snp_db.dsnp2, snp_db.nchr, snp_db.dindex
    )

    # Variants not in homology regions. Get read counts only at D6 positions.
    var_db = call_parameters.var_db
    var_alt, var_ref, var_alt_forward, var_alt_reverse = get_supporting_reads_single_region(
        bamfile, var_db.dsnp1, var_db.nchr, var_db.dindex
    )
    # Look more carefully for insertions at 42128936 from reads.
    var_list = call_parameters.var_list
    ref_read, long_ins_read, short_ins_read = get_allele_counts_var42128936(
        bamfile, call_parameters.genome
    )
    var_alt, var_ref = update_var42128936(
        var_list, var_alt, var_ref, ref_read, long_ins_read, short_ins_read
    )
    # Variants in homology regions. Get read counts at both D6/D7 positions.
    var_homo_db = call_parameters.var_homo_db
    var_homo_alt, var_homo_ref = get_supporting_reads(
        bamfile,
        var_homo_db.dsnp1,
        var_homo_db.dsnp2,
        var_homo_db.nchr,
        var_homo_db.dindex,
    )
    # This ordered dictionary is for final reporting.
    raw_count = OrderedDict()
    non_homology_variant_count = len(var_alt)
    for i in range(len(call_parameters.var_list)):
        if i < non_homology_variant_count:
            if var_list[i] in NOISY_VAR:
                raw_count.setdefault(
                    var_list[i],
                    "%i(%i:%i),%i"
                    % (var_alt[i], var_alt_forward[i], var_alt_reverse[i], var_ref[i]),
                )
            else:
                raw_count.setdefault(var_list[i], "%i,%i" % (var_alt[i], var_ref[i]))
        else:
            raw_count.setdefault(
                var_list[i],
                "%i,%i"
                % (
                    var_homo_alt[i - non_homology_variant_count],
                    var_homo_ref[i - non_homology_variant_count],
                ),
            )

    # --- MODIFIED: Original code returned no-call if GMM failed to determine CN. 
    # Remove GMM-based no-call check as total_cn is now predetermined from CNVPanelizer. ---
    # if raw_cn_call.d67_cn is None:
    #     sample_call = d6_call(...)
    #     return sample_call
    # --- END OF MODIFICATION ---

    # 4. Call CNV and hybrids
    d6_fraction = get_fraction(snp_d6, snp_d7)
    if raw_cn_call.d67_cn == 0:
        raw_d6_cn = [0.0] * len(d6_fraction)
    else:
        raw_d6_cn = [round(raw_cn_call.d67_cn * a, 3) for a in d6_fraction]

    cn_call_snp = call_cn_snp(raw_cn_call.d67_cn, snp_d6, snp_d7)

    # exon9gc
    exon9gc_call_stringent = call_exon9gc(
        snp_d6[EXON9_SITE1 : EXON9_SITE2 + 1],
        snp_d7[EXON9_SITE1 : EXON9_SITE2 + 1],
        raw_cn_call.d67_cn,
    )
    cnvtag, consensus = get_cnvtag(
        raw_cn_call.d67_cn,
        raw_d6_cn,
        cn_call_snp,
        exon9gc_call_stringent,
        raw_cn_call.spacer_cn, # This now carries predetermined spacer_cn
    )

    # no-call due to CNV group calling
    if cnvtag is None or cnvtag not in CNV_ACCEPTED:
        sample_call = d6_call(
            coverage_mad_value,
            median_depth_value,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None,
            cnvtag,
            None,
            "CNV_Group_Not_Accepted" if cnvtag not in CNV_ACCEPTED else "CNV_Group_Undetermined",
            None,
            None,
            exon9gc_call_stringent,
            ",".join(str(a) for a in consensus) if consensus else None,
            ",".join(str(a) for a in cn_call_snp),
            ",".join(str(a) for a in raw_d6_cn),
            raw_count,
        )
        if bamfile: bamfile.close()
        return sample_call

    # 5. Call variants
    # homology region
    cn_call_var_homo = call_cn_var_homo(raw_cn_call.d67_cn, var_homo_alt, var_homo_ref)
    # non-homology region
    cn_call_var = call_cn_var(
        cnvtag, var_alt, var_ref, var_alt_forward, var_alt_reverse, var_list, var_db
    )
    # call haplotypes
    haplotype_db = call_parameters.haplotype_db
    site42126938_count, var42126938, var42126938_G_haplotype = call_var42126938(
        bamfile, raw_cn_call.d67_cn, haplotype_db["g.42126938C>T"]
    )
    
    # --- MODIFIED: Use _fmt_pair() helper function for consistent count pair formatting. ---
    raw_count["g.42126938C>T"] = _fmt_pair(site42126938_count, order="d7,d6")
    # --- END OF MODIFICATION ---

    site42127526_count, site42127556_count, var42127526 = call_var42127526_var42127556(
        bamfile, cnvtag, haplotype_db["g.42127526C>T_g.42127556T>C"]
    )
    
    # --- MODIFIED: Use _fmt_pair() helper function for consistent count pair formatting. ---
    raw_count["g.42127526C>T"] = _fmt_pair(site42127526_count, order="d7,d6")
    raw_count["g.42127556T>C"] = _fmt_pair(site42127556_count, order="d7,d6")
    # --- END OF MODIFICATION ---

    var42127803_diff_haplotype = call_var42127803hap(
        bamfile, cnvtag, haplotype_db["g.42127803C>T"]
    )

    # 6. Call star allele
    total_callset = get_called_variants(var_list, cn_call_var)
    called_var_homo = get_called_variants(var_list, cn_call_var_homo, len(cn_call_var))
    total_callset += called_var_homo
    total_callset += var42126938
    total_callset += var42127526

    star_called = match_star(
        total_callset,
        cnvtag,
        raw_cn_call.spacer_cn,
        call_parameters.star_combinations,
        exon9_values(
            exon9gc_call_stringent,
            consensus.exon9_and_downstream if consensus else None,
            raw_d6_cn[EXON9_SITE1] if len(raw_d6_cn) > EXON9_SITE1 else None,
            raw_d6_cn[EXON9_SITE2] if len(raw_d6_cn) > EXON9_SITE2 else None,
        ),
        var42126938_G_haplotype,
        var42127803_diff_haplotype,
    )

    genotype_filter = None
    final_star_allele_call = None
    # no-call due to star allele matching
    if star_called.call_info and star_called.call_info != "no_match":
        final_star_allele_call = star_called.clean_call
        if final_star_allele_call:
            if ";" in final_star_allele_call:
                genotype_filter = "More_than_one_possible_genotype"
            elif "/" not in final_star_allele_call:
                genotype_filter = "Not_assigned_to_haplotypes"
            
            # --- MODIFIED: Remove GMM-based low confidence filter as it's not applicable to CNVPanelizer approach. ---
            # elif high_cn_low_confidence:
            #     genotype_filter = "LowQ_high_CN"
            # --- END OF MODIFICATION ---
            
            if genotype_filter is None : genotype_filter = "PASS"

    sample_call = d6_call(
        coverage_mad_value,
        median_depth_value,
        raw_cn_call.d67_cn,
        raw_cn_call.spacer_cn,
        raw_cn_call.d67_depth,
        raw_cn_call.spacer_depth,
        star_called.variants_called.split() if star_called.variants_called else None,
        cnvtag,
        final_star_allele_call,
        genotype_filter,
        star_called.raw_call,
        star_called.call_info,
        exon9gc_call_stringent,
        ",".join(str(a) for a in consensus) if consensus else None,
        ",".join(str(a) for a in cn_call_snp),
        ",".join(str(a) for a in raw_d6_cn),
        raw_count,
    )
    if bamfile: bamfile.close()
    return sample_call


# --- MODIFIED: Add bed_file_path parameter to support dynamic BED file selection for different panels. ---
def prepare_resource(datadir, parameters, bed_file_path):
# --- END OF MODIFICATION ---
    genome = parameters.genome
    
    # --- MODIFIED: Use dynamically detected BED file instead of hardcoded genome-specific filename.
    # This allows users to provide custom panel BED files without code modification. ---
    region_file = bed_file_path  # Use the dynamically detected BED file
    # --- END OF MODIFICATION ---
    
    snp_file = os.path.join(datadir, "CYP2D6_SNP_%s.txt" % genome)
    
    # --- MODIFIED: Remove GMM file reference as it's not used in CyriPanel. ---
    # gmm_file = os.path.join(datadir, "CYP2D6_gmm.txt")
    # --- END OF MODIFICATION ---
    
    star_table = os.path.join(datadir, "star_table.txt")
    variant_file = os.path.join(datadir, "CYP2D6_target_variant_%s.txt" % genome)
    variant_homology_file = os.path.join(
        datadir, "CYP2D6_target_variant_homology_region_%s.txt" % genome
    )
    haplotype_file = os.path.join(datadir, "CYP2D6_haplotype_%s.txt" % genome)
    star_combinations = get_hap_table(star_table)

    for required_file in [
        region_file,
        snp_file,
        variant_file,
        variant_homology_file,
        haplotype_file,
        # --- MODIFIED: Remove GMM file from validation check. ---
        # gmm_file,
        # --- END OF MODIFICATION ---
    ]:
        if not os.path.exists(required_file):
            raise Exception("File %s not found." % required_file)

    # --- MODIFIED: Pass genome parameter to get_snp_position() calls for proper coordinate handling. ---
    snp_db = get_snp_position(snp_file, genome)
    var_db = get_snp_position(variant_file, genome)
    var_homo_db = get_snp_position(variant_homology_file, genome)
    haplotype_db = {}
    for variant in HAPLOTYPE_VAR:
        haplotype_db.setdefault(variant, get_snp_position(haplotype_file, variant, genome))
    # --- END OF MODIFICATION ---
    
    var_list = []
    with open(variant_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    with open(variant_homology_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    
    # --- MODIFIED: Remove GMM parameter parsing and region_dic parsing without genome parameter. ---
    # gmm_parameter = parse_gmm_file(gmm_file)
    
    region_dic = parse_region_file(region_file, genome)

    # Remove gmm_parameter from resource_info initialization.
    call_parameters = resource_info(
        genome,
        region_dic,
        snp_db,
        var_db,
        var_homo_db,
        haplotype_db,
        var_list,
        star_combinations,
    )
    # --- END OF MODIFICATION ---
    
    return call_parameters


def main():
    parameters = load_parameters()
    manifest = parameters.manifest
    outdir = parameters.outDir
    prefix = parameters.prefix
    reference_fasta = parameters.reference
    threads = parameters.threads

    logging.basicConfig(level=logging.INFO)

    # --- MODIFIED: Set up CNVPanelizer integration paths automatically based on script location.
    # Determines paths for CNVPanelizer R script, output directory, and reference BAMs directory
    # without requiring manual configuration, assuming standard project structure. ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    depth_calling_dir = os.path.join(script_dir, "depth_calling")
    cnv_panelizer_rscript = os.path.join(depth_calling_dir, "run_CNVPanelizer.R")
    cnv_panelizer_output_dir = os.path.join(depth_calling_dir, "cnv_panelizer_results")

    # Ensure the output directory for CNVPanelizer reports exists
    os.makedirs(cnv_panelizer_output_dir, exist_ok=True)

    # Path to the reference BAMs directory
    reference_dir_path = os.path.join(script_dir, "ref_dir")

    # Count the number of BAM files in the reference directory for the log message
    try:
        num_bam_files = len([f for f in os.listdir(reference_dir_path) if f.endswith(('.bam', '.cram'))])
    except FileNotFoundError:
        num_bam_files = 0    
    logging.info(f"Expecting reference {num_bam_files} BAMs for CNVPanelizer in: {reference_dir_path}")

    if not os.path.isdir(reference_dir_path):
        logging.warning(f"Reference directory not found at: {reference_dir_path}. CNVPanelizer may fail.")

    # Automatically find and validate BED file instead of using hardcoded genome-specific filename.
    data_dir = os.path.join(script_dir, "data")
    try:
        bed_file_path = find_bed_file(data_dir)
    except Exception as e:
        logging.error(str(e))
        sys.exit(1)
    # --- END OF MODIFICATION ---

    datadir = os.path.join(os.path.dirname(__file__), "data")
    call_parameters = prepare_resource(datadir, parameters, bed_file_path)

    out_json = os.path.join(outdir, prefix + ".json")
    out_tsv = os.path.join(outdir, prefix + ".tsv")
    final_output = {}
    with open(manifest) as read_manifest:
        for line in read_manifest:
            bam_name = line.strip()
            index_name = None
            if '##idx##' in bam_name:
                bam_name, index_name = bam_name.split('##idx##')

            sample_id = os.path.splitext(os.path.basename(bam_name))[0]
            
            if "://" not in bam_name and not os.path.exists(bam_name):
                logging.warning("Input file for sample %s does not exist.", sample_id)
                continue

            logging.info(
                "Processing sample %s at %s", sample_id, datetime.datetime.now()
            )

            # --- MODIFIED: Call CNVPanelizer to determine total_cn and d7_spacer for each sample. ---
            logging.info(f"Determining CN values for {sample_id} using CNVPanelizer...")
            initial_total_cn, initial_d7_spacer = get_cn_from_cnvpanelizer(
                bam_file=bam_name,
                r_script_path=cnv_panelizer_rscript,
                output_dir=cnv_panelizer_output_dir,
                bed_file_path=bed_file_path,
                reference_dir_path=reference_dir_path
            )

            if initial_total_cn is None:
                logging.error(f"Failed to determine CN for {sample_id}. Skipping.")
                cyp2d6_call = {
                    "Genotype": "Error",
                    "Filter": "CN_Determination_Failed"
                }
                final_output[sample_id] = cyp2d6_call
                continue
            # --- END OF MODIFICATION ---
                        
            # --- MODIFIED: Implement comprehensive three-tier fallback mechanism for robust genotype calling in targeted sequencing.
            # Tier 1: Initial attempt with CNVPanelizer values
            # Tier 2: Strategic CN adjustment (±1-3 from initial, keeping spacer constant) to address CNVPanelizer under/overestimation
            # Tier 3: Comprehensive search through FALLBACK_CN_COMBINATIONS prioritized by proximity to initial CN
            # Accepts "PASS" or "More_than_one_possible_genotype" as success. Stores best "Not_assigned_to_haplotypes" as fallback. ---
            
            # 1. Initial attempt: First, try to call the genotype using the predetermined --total_cn and --d7_spacer values. 
            logging.info(f"Initial attempt for {sample_id} with total_cn={initial_total_cn} and d7_spacer={initial_d7_spacer}")
            cyp2d6_call = d6_star_caller(
                bam_name,
                call_parameters,
                threads,
                initial_total_cn,
                initial_d7_spacer,
                reference_fasta=reference_fasta,
                index_name=index_name
            )._asdict()

            # 2. Check for a high-confidence successful result: If the initial attempt results in a "PASS" or "More_than_one_possible_genotype" filter, consider it a success and store the result.
            if cyp2d6_call.get("Filter") in ["PASS", "More_than_one_possible_genotype"]:
                logging.info(f"SUCCESS: Initial parameters worked for {sample_id} with Filter: {cyp2d6_call.get('Filter')}.")
                final_output[sample_id] = cyp2d6_call
            else:
                # If the initial attempt fails, start the enhanced fallback process.
                logging.warning(
                    f"Initial attempt failed for {sample_id} with Filter: {cyp2d6_call.get('Filter')}. "
                    f"Starting enhanced fallback."
                )
                
                initial_failed_call = cyp2d6_call # Store the very first result.
                found_success = False
                
                # ENHANCED FALLBACK STRATEGY
                # Strategy 1: Try INCREMENTING total_cn while keeping spacer_cn constant
                # This addresses CNVPanelizer underestimation issues
                # Try both increasing AND decreasing total_cn
                if initial_d7_spacer is not None:
                    # First try incrementing (more common case)
                    for cn_increment in [1, 2, 3]:  # Try +1, +2, and +3
                        adjusted_total_cn = initial_total_cn + cn_increment
                        
                        # Skip if this would result in an impossible configuration
                        if adjusted_total_cn > 10:  # Reasonable upper limit
                            continue
                            
                        logging.info(
                            f"Trying INCREASED total_cn={adjusted_total_cn} with same spacer_cn={initial_d7_spacer}"
                        )
                        
                        test_call = d6_star_caller(
                            bam_name, call_parameters, threads,
                            adjusted_total_cn, initial_d7_spacer,
                            reference_fasta=reference_fasta, index_name=index_name
                        )._asdict()
                        
                        if test_call.get("Filter") in ["PASS", "More_than_one_possible_genotype"]:
                            # Verify this makes sense
                            expected_exon9_cn = adjusted_total_cn - initial_d7_spacer
                            logging.info(
                                f"SUCCESS with increased total_cn: total_cn={adjusted_total_cn}, "
                                f"spacer_cn={initial_d7_spacer}, expected_exon9_cn={expected_exon9_cn}"
                            )
                            final_output[sample_id] = test_call
                            found_success = True
                            break
                    
                    # If incrementing didn't work, try decrementing (less common but possible)
                    if not found_success and initial_total_cn > 2:
                        for cn_decrement in [1, 2]:  # Try -1 and -2
                            adjusted_total_cn = initial_total_cn - cn_decrement
                            
                            # Skip if this would result in impossible configuration
                            if adjusted_total_cn < 2 or adjusted_total_cn <= initial_d7_spacer:
                                continue
                            
                            logging.info(
                                f"Trying DECREASED total_cn={adjusted_total_cn} with same spacer_cn={initial_d7_spacer}"
                            )
                            
                            test_call = d6_star_caller(
                                bam_name, call_parameters, threads,
                                adjusted_total_cn, initial_d7_spacer,
                                reference_fasta=reference_fasta, index_name=index_name
                            )._asdict()
                            
                            if test_call.get("Filter") in ["PASS", "More_than_one_possible_genotype"]:
                                expected_exon9_cn = adjusted_total_cn - initial_d7_spacer
                                logging.info(
                                    f"SUCCESS with decreased total_cn: total_cn={adjusted_total_cn}, "
                                    f"spacer_cn={initial_d7_spacer}, expected_exon9_cn={expected_exon9_cn}"
                                )
                                final_output[sample_id] = test_call
                                found_success = True
                                break
                
                # Strategy 2: If Strategy 1 fails, fall back to original comprehensive search
                if not found_success:
                    logging.info("Same spacer_cn strategy failed. Trying comprehensive fallback combinations.")
                    
                    # These variables will track the best "Not_assigned_to_haplotypes" call.
                    best_not_assigned_call = None
                    min_distance = float('inf') # Initialize with a very large number.
                    
                    # Sort the fallback combinations with a smarter ordering:
                    # 1. First prioritize combinations close to initial total_cn
                    # 2. Within same distance, prefer larger total_cn (for underestimation cases)
                    fallback_order = sorted(
                        FALLBACK_CN_COMBINATIONS, 
                        key=lambda x: (abs(x[0] - initial_total_cn), -x[0])  # Distance first, then prefer larger CN
                    )

                    # 3. Begin primary fallback loop: Iterate through the sorted list of pre-defined CN combinations to find a best working set.
                    for total_cn_fallback, d7_spacer_fallback in fallback_order:
                        # Skip the combination if it's the same as the initial attempt.
                        if total_cn_fallback == initial_total_cn and d7_spacer_fallback == initial_d7_spacer:
                            continue
                        # Also skip combinations we already tried in Strategy 1
                        if d7_spacer_fallback == initial_d7_spacer and abs(total_cn_fallback - initial_total_cn) <= 3:
                            continue

                        logging.info(f"Comprehensive fallback for {sample_id}: total_cn={total_cn_fallback}, d7_spacer={d7_spacer_fallback}")
                        
                        fallback_call_dict = d6_star_caller(
                            bam_name, call_parameters, threads, total_cn_fallback,
                            d7_spacer_fallback, reference_fasta=reference_fasta, index_name=index_name
                        )._asdict()

                        # Check for high-confidence success.
                        if fallback_call_dict.get("Filter") in ["PASS", "More_than_one_possible_genotype"]:
                            logging.info(f"SUCCESS: Comprehensive fallback successful for {sample_id} with Filter: {fallback_call_dict.get('Filter')}.")
                            final_output[sample_id] = fallback_call_dict
                            found_success = True
                            break # Exit loop on first high-confidence success.

                        # Logic for the 'Not_assigned_to_haplotypes' fallback: If a high-confidence call isn't found, this looks for the "best" lower-confidence call.
                        if fallback_call_dict.get("Filter") == "Not_assigned_to_haplotypes":
                            distance = abs(total_cn_fallback - initial_total_cn)
                            # If this call's CN is closer to the initial guess than any previous one, store it.
                            if distance < min_distance:
                                min_distance = distance
                                best_not_assigned_call = fallback_call_dict
                                logging.info(f"Found a candidate 'Not_assigned_to_haplotypes' result with distance {distance}. Storing as best option so far.")

                # 4. Determine the final result after the loop is complete.
                if found_success:
                    pass # The successful result is already in final_output.
                elif best_not_assigned_call is not None:
                    # If no "PASS" call was found, use the best 'Not_assigned_to_haplotypes' call found during the loop.
                    logging.warning(f"FALLBACK TIER 2: Using best 'Not_assigned_to_haplotypes' result for {sample_id}.")
                    final_output[sample_id] = best_not_assigned_call
                else:
                    # If all fallbacks fail to produce a usable genotype, report a complete failure with None.
                    logging.error(f"COMPLETE FAILURE: No usable genotype found for {sample_id}.")
                    initial_failed_call['Genotype'] = None
                    initial_failed_call['Filter'] = None # Set Filter to None as requested
                    final_output[sample_id] = initial_failed_call
            # --- END OF MODIFICATION ---

    # Write to json
    logging.info("Writing to json at %s", datetime.datetime.now())
    with open(out_json, "w") as json_output:
        json.dump(final_output, json_output, indent=4)

    # Write to tsv
    logging.info("Writing to tsv at %s", datetime.datetime.now())
    header = ["Sample", "Genotype", "Filter"]
    with open(out_tsv, "w") as tsv_output:
        tsv_output.write("\t".join(header) + "\n")
        for sample_id_out, final_call in final_output.items():
            output_per_sample = [
                sample_id_out,
                str(final_call.get("Genotype", "None")),
                str(final_call.get("Filter", "None")),
            ]
            tsv_output.write("\t".join(output_per_sample) + "\n")

if __name__ == "__main__":
    main()
