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


from scipy.stats import poisson

POSTERIOR_CUTOFF_STRINGENT = 0.9
ERROR_RATE = 0.1 # --- BCyirus: Modified from 0.01 to 0.1 ---
# --- MODIFIED: For handling low allele fractions for CN=1 in targeted panels.
# Add a new variable to check if the total copy number (full_cn) is 1.
# If the total copy number is 1, the allele fraction is checked against a threshold (MIN_ALLELE_FRACTION_FOR_CN1_TARGETED).
# This helps to avoid incorrect calls due to low-level noise that might otherwise be misinterpreted. ---
MIN_ALLELE_FRACTION_FOR_CN1_TARGETED = 0.15
# --- MODIFIED: For handling ambiguous CN=2-vs-3 calls in targeted panels.
# Add a new variable to check the confidence when calling CN=3. If the confidence is below this threshold, the call is preferentially set to CN=2. 
# This helps resolve ambiguity in targeted panel data where signals can be noisy. ---
AMBIGUOUS_CALL_THRESHOLD = 0.85
# --- END OF MODIFICATION ---


def call_reg1_cn(full_cn, count_reg1, count_reg2, min_read=0):
    """
    Return the reg1 copy number call at each site based on Poisson likelihood,
    with a minimum read support cutoff
    """
    if full_cn is None:
        return [None]
    if full_cn == 0:
        return [0]
    # --- MODIFIED: To more confidently call a copy number of 1, which is particularly useful for targeted sequencing data where read counts can be variable. 
    # Before this change, the function relied solely on a Poisson distribution model to determine the most likely copy number. 
    # While effective, this Poisson distribution model can sometimes be ambiguous when read counts are low or noisy, especially in targeted sequencing data.
    # Now, the function also checks the allele fraction to ensure that the call is not due to low-level noise. ---
    total_reads_at_site = count_reg1 + count_reg2

    if full_cn == 1:
        if total_reads_at_site > 20:  # Only apply this logic with sufficient read depth
            allele_fraction = count_reg1 / total_reads_at_site
            
            # CONFLICT DETECTION: If the fraction suggests heterozygosity,
            # it contradicts the full_cn=1 assumption. Return None to indicate a conflict.
            if 0.25 < allele_fraction < 0.75:
                return [None]

            # CONFIDENT CALLS: If the fraction is overwhelmingly high or low, make a confident call.
            # This is much stricter than the original logic.
            if allele_fraction >= 0.85 and count_reg1 > min_read:
                return [1]  # Confidently a variant on the single copy
            if allele_fraction <= 0.15 and count_reg2 > min_read:
                return [0]  # Confidently not a variant (reference) on the single copy
    # --- END OF MODIFICATION ---
    prob = []
    nsum = total_reads_at_site
    if nsum == 0:
        return [None]
    for i in range(full_cn + 1):
        depthexpected = float(nsum) * float(i) / float(full_cn)
        if i == 0:
            depthexpected = (ERROR_RATE / 3) * float(nsum)
        if i == full_cn:
            depthexpected = float(nsum) - ERROR_RATE * float(nsum)
        if count_reg1 <= count_reg2:
            prob.append(poisson.pmf(int(count_reg1), depthexpected))
        else:
            prob.append(poisson.pmf(int(count_reg2), depthexpected))
    sum_prob = sum(prob)
    if sum_prob == 0:
        return [None]
    post_prob = [float(a) / float(sum_prob) for a in prob]
    if count_reg2 < count_reg1:
        post_prob = post_prob[::-1]
    post_prob_sorted = sorted(post_prob, reverse=True)
    if (
        post_prob.index(post_prob_sorted[0]) != 0
        and count_reg1 <= min_read
        and count_reg2 >= min_read
    ):
        return [0]
    if post_prob_sorted[0] >= POSTERIOR_CUTOFF_STRINGENT:
        return [post_prob.index(post_prob_sorted[0])]
    # output the two most likely scenarios
    cn_prob_filtered = [
        post_prob.index(post_prob_sorted[0]),
        round(post_prob_sorted[0], 3),
        post_prob.index(post_prob_sorted[1]),
        round(post_prob_sorted[1], 3),
    ]
    return cn_prob_filtered


def process_raw_call_gc(cn_prob, post_cutoff, keep_none=True):
    """
    Filter raw CN calls based on posterior probablity cutoff.
    For gene conversion cases, i.e. SNVs between paralogs.
    -- MODIFIED to handle ambiguous 2-vs-3 calls for targeted panels --
    """
    cn_prob_filtered = []
    for cn_call in cn_prob:
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        # --- MODIFIED: For handling ambiguous 2-vs-3 calls in targeted panels.
        # Add a new variable to check the confidence when calling CN=3. If the confidence is below this threshold, the call is preferentially set to CN=2. 
        # This helps resolve ambiguity in targeted panel data where signals can be noisy ---
        elif (len(cn_call) > 2 and cn_call[0] == 3 and cn_call[2] == 2 and cn_call[1] < AMBIGUOUS_CALL_THRESHOLD):
            cn_prob_filtered.append(2)  # Preferentially call CN=2
        # --- END OF MODIFICATION ---
        elif cn_call[1] > post_cutoff:  # Original logic
            cn_prob_filtered.append(cn_call[0])
        elif keep_none:
            cn_prob_filtered.append(None)
    return cn_prob_filtered


def process_raw_call_denovo(
    cn_prob, post_cutoff1, post_cutoff2, list_total_cn=None, keep_none=True
):
    """
    Filter raw CN calls based on posterior probablity cutoff.
    For de novel variant calling, i.e. non-gene-conversion cases.
    For less confident calls that are not copy number zero,
    return the smaller CN call.
    Also keep the variant if called CN is equal to total CN at the site.
    This makes sure a variant is always kept if CN>=1
    but we can't distinguish 1 vs 2.
    """
    cn_prob_filtered = []
    for i, cn_call in enumerate(cn_prob):
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        else:
            current_total_cn_at_site = None
            if list_total_cn is not None:
                current_total_cn_at_site = list_total_cn[i]
            keep_var_conditionally = False
            if current_total_cn_at_site is not None:
                keep_var_conditionally = (cn_call[0] > 0 and cn_call[2] > 0) or \
                                         (cn_call[0] == current_total_cn_at_site or cn_call[2] == current_total_cn_at_site)
            else:
                keep_var_conditionally = cn_call[0] > 0 and cn_call[2] > 0
            if cn_call[1] > post_cutoff1:
                cn_prob_filtered.append(cn_call[0])
            elif keep_var_conditionally:
                if cn_call[1] > post_cutoff2:
                    cn_prob_filtered.append(cn_call[0])
                else:
                    cn_prob_filtered.append(min(cn_call[0], cn_call[2]))
            elif keep_none:
                cn_prob_filtered.append(None)
    return cn_prob_filtered
