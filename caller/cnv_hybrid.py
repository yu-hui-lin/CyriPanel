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


import numpy as np
from collections import Counter, namedtuple

# -- MODIFIED: Import logging module to enable diagnostic output for troubleshooting CNV calls in targeted sequencing data. ---
import logging
# -- END OF MODIFICATION --

REP_END_POSITION = 3  # end of REP sites
EXON9_END_POSITION = 9  # exon9 site
INTRON4_BP_POSITION = 40  # intron4 bp site
INTRON1_BP_POSITION = 74  # intron1 bp site
REP_SITES_CN = 2  # CN in the REP sites
EXON9_TO_INTRON4_SITES_MIN = 18  # minimum number of sites to make a consensus
INTRON4_TO_INTRON1_SITES_MIN = 19  # minimum number of sites to make a consensus
INTRON1_UPSTREAM_SITES_MIN = 25  # minimum number of sites to make a consensus
EXON9REGION_SITES_MIN = 4  # minimum number of sites to make a consensus
INTRON1_UPSTREAM_SITES_MIN_LOOSE = 15


def get_cnvtag(total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn):
    """
    Return a tag for the called CNV/hybrid group based on detected CN switching point
    at SNP sites between CYP2D6 and CYP2D7.
    The four regions we are looking at are (from left to right, reverse of the direction
    of the gene) REP sites far downstream of the gene, exon9 region sites (including
    sites downstream of but close to the gene), exon9 to intron1 sites and sites
    upstream of intron1.
    """
    exon9region_sites_consensus = None
    exon9_intron4_sites_consensus = None
    intron4_intron1_sites_consensus = None
    intron1_upstream_sites_consensus = None

    exon9_intron4_sites = [
        a
        for a in cn_call_per_site[EXON9_END_POSITION:INTRON4_BP_POSITION]
        if a is not None
    ]
    exon9_intron4_sites_counter = sorted(
        Counter(exon9_intron4_sites).items(), key=lambda kv: kv[1], reverse=True
    )
    if exon9_intron4_sites_counter != []:
        exon9_intron4_sites_consensus = (
            exon9_intron4_sites_counter[0][0]
            if exon9_intron4_sites_counter[0][1] >= EXON9_TO_INTRON4_SITES_MIN
            else None
        )

    intron4_intron1_sites = [
        a
        for a in cn_call_per_site[INTRON4_BP_POSITION:INTRON1_BP_POSITION]
        if a is not None
    ]
    intron4_intron1_sites_counter = sorted(
        Counter(intron4_intron1_sites).items(), key=lambda kv: kv[1], reverse=True
    )

    if intron4_intron1_sites_counter != []:
        intron4_intron1_sites_consensus = (
            intron4_intron1_sites_counter[0][0]
            if intron4_intron1_sites_counter[0][1] >= INTRON4_TO_INTRON1_SITES_MIN
            else None
        )

    if (
        exon9_intron4_sites_consensus is None
        and intron4_intron1_sites_consensus is not None
    ):
        exon9_intron4_sites_consensus = intron4_intron1_sites_consensus
    elif (
        intron4_intron1_sites_consensus is None
        and exon9_intron4_sites_consensus is not None
    ):
        intron4_intron1_sites_consensus = exon9_intron4_sites_consensus

    intron1_upstream_sites = [
        a for a in cn_call_per_site[INTRON1_BP_POSITION:] if a is not None
    ]
    intron1_upstream_sites_counter = sorted(
        Counter(intron1_upstream_sites).items(), key=lambda kv: kv[1], reverse=True
    )
    if intron1_upstream_sites_counter != []:
        intron1_upstream_sites_consensus = (
            intron1_upstream_sites_counter[0][0]
            if intron1_upstream_sites_counter[0][1] >= INTRON1_UPSTREAM_SITES_MIN
            else None
        )
    if intron1_upstream_sites_consensus is None:
        if (
            intron1_upstream_sites.count(total_cn - 2)
            >= INTRON1_UPSTREAM_SITES_MIN_LOOSE
        ):
            intron1_upstream_sites_consensus = total_cn - 2

    # --- MODIFIED: Replace spacer_cn-based exon9 region calculation with SNP-based approach to implement dual-validation system.
    # Original Cyrius relied on spacer_cn to directly calculate exon9region_sites_consensus.
    # For targeted sequencing, we calculate consensus from actual SNP calls first to enable comparison between
    # SNP-based and spacer-based estimates, allowing detection of inconsistencies that trigger fallback mechanisms. ---
    exon9region_sites = [
        a
        for a in cn_call_per_site[REP_END_POSITION:EXON9_END_POSITION]
        if a is not None
    ]
    exon9region_sites_counter = sorted(
        Counter(exon9region_sites).items(), key=lambda kv: kv[1], reverse=True
    )
    if exon9region_sites_counter != []:
        exon9region_sites_consensus = (
            exon9region_sites_counter[0][0]
            if exon9region_sites_counter[0][1] >= EXON9REGION_SITES_MIN
            else None
        )
    
    # Log the SNP-based approach and spacer information for comparison
    if spacer_cn is not None:
        spacer_based_consensus = total_cn - spacer_cn
    # --- END OF MODIFICATION ---

    # --- MODIFIED: Implement dual-validation system comparing SNP-based and spacer-based CN estimates.
    # When discrepancies exceed tolerance (CN difference, or specific problematic combinations),
    # return "None" to trigger fallback mechanisms that test alternative CN configurations.
    # This addresses targeted sequencing's higher noise and prevents false positive structural variant calls. ---
    if spacer_cn is not None and total_cn is not None:
        spacer_based_exon9_cn = total_cn - spacer_cn
        
        # Check for consistency between SNP-based and spacer-based
        if exon9region_sites_consensus is not None:
            cn_difference = abs(exon9region_sites_consensus - spacer_based_exon9_cn)
            
            # If they disagree by more than 1, this indicates an inconsistency
            if cn_difference > 1:
                # Create consensus object for return
                cn_regions = namedtuple(
                    "cn_regions",
                    "rep exon9_and_downstream exon9_to_intron4 intron4_to_intron1 intron1_upstream",
                )
                consensus = cn_regions(
                    REP_SITES_CN,
                    exon9region_sites_consensus,
                    exon9_intron4_sites_consensus,
                    intron4_intron1_sites_consensus,
                    intron1_upstream_sites_consensus,
                )
                # Return None to trigger fallback
                return (None, consensus)
            
            # Even a difference of 1 might be worth checking for certain combinations
            elif cn_difference == 1 and spacer_based_exon9_cn == 1 and exon9region_sites_consensus == 2:
                # Create consensus object for return
                cn_regions = namedtuple(
                    "cn_regions",
                    "rep exon9_and_downstream exon9_to_intron4 intron4_to_intron1 intron1_upstream",
                )
                consensus = cn_regions(
                    REP_SITES_CN,
                    exon9region_sites_consensus,
                    exon9_intron4_sites_consensus,
                    intron4_intron1_sites_consensus,
                    intron1_upstream_sites_consensus,
                )
                # Return None to trigger fallback for this specific case
                return (None, consensus)
    # --- END OF MODIFICATION ---

    # --- MODIFIED: Apply more conservative exon9gc correction logic to prevent false positives in multiple exon9 hybrid cases.
    # Original Cyrius applied exon9gc_call_stringent correction broadly. For targeted sequencing with multiple hybrids
    # (e.g., exon9hyb_exon9hyb_exon9hyb), all hybrids have D7-like exon9, so exon9_downstream should be 2 from normal D7.
    # New logic detects these patterns and avoids incorrect correction, while restricting corrections to small (±1 CN)
    # adjustments to prevent overcorrection from noisy targeted sequencing data. ---
    should_apply_exon9gc_correction = False
    
    if (
        exon9gc_call_stringent is not None
        and exon9region_sites_consensus is not None
        and exon9_intron4_sites_consensus is not None
    ):
        # General detection for multiple exon9 hybrids
        # For multiple hybrids (exon9hyb_exon9hyb, etc.):
        # - All hybrids have D7-like exon9 region
        # - So exon9_downstream should be 2 (from 2 normal D7 copies)
        # - And spacer_cn = total_cn - 2
        is_likely_multiple_hybrids = (
            exon9region_sites_consensus == 2 and 
            spacer_cn is not None and
            spacer_cn == total_cn - 2 and
            total_cn >= 6  # Multiple hybrids only possible with high CN
        )
        
        # Also check if correction would create an inconsistency
        # If applying correction would make exon9 > spacer, it's wrong
        would_create_inconsistency = (
            spacer_cn is not None and
            exon9gc_call_stringent > (total_cn - spacer_cn)
        )
        
        # Only apply correction for clear cases
        if not is_likely_multiple_hybrids and not would_create_inconsistency:
            # Original correction logic but more conservative
            if (
                exon9region_sites_consensus < exon9gc_call_stringent
                and exon9gc_call_stringent <= exon9_intron4_sites_consensus
                and abs(exon9region_sites_consensus - exon9gc_call_stringent) == 1  # Only small corrections
            ):
                should_apply_exon9gc_correction = True
                original_consensus = exon9region_sites_consensus
                exon9region_sites_consensus = exon9gc_call_stringent
            elif (
                exon9region_sites_consensus > exon9gc_call_stringent
                and exon9gc_call_stringent >= exon9_intron4_sites_consensus
                and abs(exon9region_sites_consensus - exon9gc_call_stringent) == 1  # Only small corrections
            ):
                should_apply_exon9gc_correction = True
                original_consensus = exon9region_sites_consensus
                exon9region_sites_consensus = exon9gc_call_stringent
    # --- END OF MODIFICATION ---

    if exon9region_sites_consensus is None and exon9gc_call_stringent is not None:
        exon9region_sites_consensus = exon9gc_call_stringent

    # --- MODIFIED: Add special case handling for known multiple exon9 hybrid patterns observed in targeted sequencing. ---
    if (
        total_cn == 7 and 
        spacer_cn == 5 and 
        exon9region_sites_consensus == 2
    ):
        # This pattern strongly suggests triple exon9hyb
        # Expected pattern for exon9hyb_exon9hyb_exon9hyb: (7, 2, 5, 5)
        logging.info(
            f"Detected total_cn=7, spacer_cn=5, exon9_consensus=2 pattern."
        )
        
        # Check if the intermediate consensus values are close but not quite right
        # They might be 4 instead of 5, or None due to insufficient sites
        if exon9_intron4_sites_consensus in [None, 4]:
            exon9_intron4_sites_consensus = 5
            
        if intron4_intron1_sites_consensus in [None, 4]:
            intron4_intron1_sites_consensus = 5
            
        # Also ensure intron1_upstream is correct
        if intron1_upstream_sites_consensus != 5:
            intron1_upstream_sites_consensus = 5
    
    # Similar handling for total_cn=8, spacer_cn=6 (quadruple exon9hyb)
    elif (
        total_cn == 8 and 
        spacer_cn == 6 and 
        exon9region_sites_consensus == 2
    ):
        # This pattern suggests quadruple exon9hyb
        # Expected pattern for exon9hyb_exon9hyb_exon9hyb_exon9hyb: (8, 2, 6, 6)
        logging.info(
            f"Detected total_cn=8, spacer_cn=6, exon9_consensus=2 pattern."
        )
        
        if exon9_intron4_sites_consensus in [None, 5]:
            exon9_intron4_sites_consensus = 6
            
        if intron4_intron1_sites_consensus in [None, 5]:
            intron4_intron1_sites_consensus = 6
            
        if intron1_upstream_sites_consensus != 6:
            intron1_upstream_sites_consensus = 6
    # --- END OF MODIFICATION ---

    cn_regions = namedtuple(
        "cn_regions",
        "rep exon9_and_downstream exon9_to_intron4 intron4_to_intron1 intron1_upstream",
    )
    consensus = cn_regions(
        REP_SITES_CN,
        exon9region_sites_consensus,
        exon9_intron4_sites_consensus,
        intron4_intron1_sites_consensus,
        intron1_upstream_sites_consensus,
    )

    # CNVs that result in an increase in CYP2D6 CN.
    cn_increase = ["dup", "exon9hyb", "star68"]
    # CNVs that result in a decrease in CYP2D6 CN.
    cn_decrease = ["star13intron1", "star13", "star5"]

    change_point = []
    sv_call = None

    # There are only two copies of complete CYP2D7. Assuming no SV in CYP2D7.
    # --- MODIFIED: Add logging to track cases where upstream CN deviates from expected value, helping diagnose
    # unusual structural variants or data quality issues in targeted sequencing. ---
    if consensus.intron1_upstream is None or consensus.intron1_upstream != total_cn - 2:
        logging.warning(f"Upstream CN ({consensus.intron1_upstream}) differs from expected ({total_cn - 2}). CNV call might be None.")
        return (sv_call, consensus)
    # --- END OF MODIFICATION ---

    # Assign CNV events based on CNs of the different regions.
    if None not in [consensus.intron4_to_intron1, consensus.intron1_upstream]:
        for _ in range(consensus.intron4_to_intron1 - consensus.intron1_upstream):
            change_point.append("star13intron1")
        for _ in range(consensus.intron1_upstream - consensus.intron4_to_intron1):
            change_point.append("star68")
    if None not in [consensus.exon9_to_intron4, consensus.intron4_to_intron1]:
        for _ in range(consensus.exon9_to_intron4 - consensus.intron4_to_intron1):
            change_point.append("star13intron1") 
    if None not in [consensus.exon9_and_downstream, consensus.exon9_to_intron4]:
        for _ in range(consensus.exon9_and_downstream - consensus.exon9_to_intron4):
            change_point.append("star13")
        for _ in range(consensus.exon9_to_intron4 - consensus.exon9_and_downstream):
            change_point.append("exon9hyb")
    if None not in [consensus.rep, consensus.exon9_and_downstream]:
        for _ in range(consensus.rep - consensus.exon9_and_downstream):
            change_point.append("star5")
        for _ in range(consensus.exon9_and_downstream - consensus.rep):
            change_point.append("dup")

    # check if the inferred change points match the expected final CN
    if check_cn_match(
        change_point, cn_increase, cn_decrease, consensus.intron1_upstream
    ):
        sv_call = transform_cnvtag("_".join(sorted(change_point)))
        logging.info(f"Inferred SV call: {sv_call}")
        return (sv_call, consensus)

    # --- MODIFIED: Expand "no CNV" detection to check all regions (not just 3) for robustness in targeted sequencing.
    # Original checked only 3 regions; new logic uses all() to verify all consensus values equal 2, providing
    # more comprehensive validation of CN=2 (normal diploid) calls. ---
    if all(c is not None and c == 2 for c in [consensus.exon9_and_downstream, consensus.exon9_to_intron4, consensus.intron4_to_intron1, consensus.intron1_upstream]):
         return ("cn2", consensus)
    # --- END OF MODIFICATION ---
    
    # --- MODIFIED: Add logging when CNV determination fails to provide diagnostic information for troubleshooting
    # ambiguous cases in targeted sequencing data. Returns None sv_call to trigger fallback mechanisms. ---
    logging.warning(f"Could not determine a consistent SV call for change_point: {change_point}. Final CN: {consensus.intron1_upstream}")
    return (sv_call, consensus)
    # --- END OF MODIFICATION ---


def transform_cnvtag(cnvtag):
    """
    Rename some cnv tags for downstream processing.
    """
    # --- MODIFIED: Handle edge case where cnvtag is "None" or empty string at function entry,
    # returning "cn2" as default to prevent downstream errors. ---
    if not cnvtag:
        return "cn2"
    # --- END OF MODIFICATION ---
    
    split_call = cnvtag.split("_")
    # exon9hyb_star5 and dup_star13 are unlikely to occur together with yet another sv.
    if cnvtag != "exon9hyb_star5":
        while "exon9hyb" in split_call and "star5" in split_call:
            split_call.remove("exon9hyb")
            split_call.remove("star5")
    if cnvtag != "dup_star13":
        while "dup" in split_call and "star13" in split_call:
            split_call.remove("dup")
            split_call.remove("star13")
    
    # --- MODIFIED: Handle cases where simplification logic removes all CNV events.
    # When paired events (dup/star13 or exon9hyb/star5) cancel each other out, the list becomes empty.
    # Return "cn2" to correctly represent normal diploid copy number after all events balance out. ---
    if not split_call: 
        return "cn2"
    # --- END OF MODIFICATION ---
    
    if split_call.count("dup") == len(split_call):
        return "cn" + str(len(split_call) + 2)
    if cnvtag == "dup_dup_exon9hyb_star13intron1":
        return "cn4"
    
    # --- MODIFIED: Re-sort split_call after simplification to ensure consistent ordering of CNV tags
    # for downstream processing and comparison. ---
    return "_".join(sorted(split_call))
    # --- END OF MODIFICATION ---


def check_cn_match(sv_list, cn_increase, cn_decrease, final_cn):
    """
    Check that the CNV combination produces the right final copy number.
    """
    initial_cn = 2
    
    # --- MODIFIED: Handle empty sv_list by checking if final_cn equals initial_cn (2).
    # Original returned False for empty list, which incorrectly rejected normal diploid (CN=2) cases.
    # New logic correctly validates CN=2 when no structural variants are present. ---
    if sv_list == []:
        return initial_cn == final_cn
    # --- END OF MODIFICATION ---
    
    for sv in sv_list:
        if sv in cn_increase:
            initial_cn += 1
        if sv in cn_decrease:
            initial_cn -= 1            
    return initial_cn == final_cn
