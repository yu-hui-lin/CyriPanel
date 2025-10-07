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
from collections import namedtuple
from scipy.stats import poisson, fisher_exact
import pysam
import numpy as np

dir_name = os.path.join(os.path.dirname(os.path.dirname(__file__)), "depth_calling")
if os.path.exists(dir_name):
    sys.path.append(dir_name)
from depth_calling.copy_number_call import (
    call_reg1_cn,
    process_raw_call_gc,
    process_raw_call_denovo,
)
from depth_calling.haplotype import (
    get_haplotypes_from_bam,
    get_haplotypes_from_bam_single_region,
    extract_hap,
)
from depth_calling.snp_count import (
    get_supporting_reads,
    get_supporting_reads_single_region,
)


INTRON1_BP_APPROX = 42130500
EXON9_BP_APPROX = 42126611

# --- MODIFIED: Relax strand bias filtering thresholds for targeted sequencing with high depth.
# Original P_CUTOFF=0.05 was too strict for panels generating very high depth (>1000x),
# where statistically significant p-values can arise from minor proportional deviations that don't indicate true bias.
# New P_CUTOFF=0.001 with MIN_ACCEPTABLE_STRAND_RATIO allows variants to pass
# if minor strand has ≥15-20% of major strand reads, even when Fisher's exact test is significant. ---
P_CUTOFF = 0.001
MIN_ACCEPTABLE_STRAND_RATIO_FOR_NOISY = 0.15
# --- END OF MODIFICATION ---

cn_regions = namedtuple(
    "cn_regions", "total_cn exon9_and_downstream exon9_to_intron1 intron1_upstream"
)
CNVTAG_LOOKUP_TABLE = {
    "star5_star5": cn_regions(2, 0, 0, 0),
    "star13_star13": cn_regions(2, 2, 0, 0),
    "star13intron1_star13intron1": cn_regions(2, 2, 2, 0),
    "star5": cn_regions(3, 1, 1, 1),
    "star13": cn_regions(3, 2, 1, 1),
    "star13intron1": cn_regions(3, 2, 2, 1),
    "star5_star5_star68": cn_regions(3, 0, 0, 1),
    "star5_star68": cn_regions(4, 1, 1, 2),
    "cn2": cn_regions(4, 2, 2, 2),    
    # --- MODIFIED: Process exon9hyb_star5 as cn2 for improved compatibility, following BCyrius approach. ---
    #'exon9hyb_star5': cn_regions(4, 1, 2, 2),
    "exon9hyb_star5": cn_regions(4, 2, 2, 2), 
    # --- END OF MODIFICATION ---    
    "dup_star13": cn_regions(4, 3, 2, 2),
    "dup_star13intron1": cn_regions(4, 3, 3, 2),
    "star13_star68": cn_regions(4, 2, 1, 2),
    "cn3": cn_regions(5, 3, 3, 3),
    "exon9hyb": cn_regions(5, 2, 3, 3),
    "star68": cn_regions(5, 2, 2, 3),
    "cn4": cn_regions(6, 4, 4, 4),
    "exon9hyb_exon9hyb": cn_regions(6, 2, 4, 4),
    "star68_star68": cn_regions(6, 2, 2, 4),
    "dup_exon9hyb": cn_regions(6, 3, 4, 4),
    "dup_star68": cn_regions(6, 3, 3, 4),
    "exon9hyb_star68": cn_regions(6, 2, 3, 4),
    "cn5": cn_regions(7, 5, 5, 5),
    "exon9hyb_exon9hyb_exon9hyb": cn_regions(7, 2, 5, 5),
    "star68_star68_star68": cn_regions(7, 2, 2, 5),
    "cn6": cn_regions(8, 6, 6, 6),
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb": cn_regions(8, 2, 6, 6),
    "star68_star68_star68_star68": cn_regions(8, 2, 2, 6),
}

# For these variants the region is clean and we use a less stringent minimum read cutoff
CLEAN_VAR = [
    "g.42129809T>C",
    "g.42129819G>T",
    "g.42128945C>T",
    "g.42126611C>G",
    "g.42130692G>A",
    "g.42127941G>A",
]

# --- MODIFIED: Remove g.42128181A>T and g.42128185C>T from NOISY_VAR list following BCyrius update.
# These variants were found to not require special noisy handling in updated star allele definitions. ---
# These are noisy (mostly gene conversion) variants that may have misalignments
# resulting in strand bias
NOISY_VAR = [
    "g.42127473C>T",
    # "g.42128181A>T",  # BCyrius: Remove this variant from NOISY_VAR
    # "g.42128185C>T",  # BCyrius: Remove this variant from NOISY_VAR
    "g.42129042T>C",
    "g.42129174C>A",
    "g.42129180A>T",
    "g.42127526C>T",
    "g.42128325A>G",
    "g.42126877G>A",
    "g.42127973T>C",
    "g.42127556T>C",
]
# --- END OF MODIFICATION ---


def get_total_cn_per_site(cnvtag, var_db, var_list):
    """
    For variants in non-homology regions, get the total expected CN
    at each site based on the CNV configuration.
    """
    num_var_sites = len(var_db.dsnp1)
    variant_names = var_list[:num_var_sites]
    num_var_sites_before_intron1bp = 0
    num_var_sites_after_exon9bp = 0
    for var_name in variant_names:
        var_pos = int(var_name[2:10])
        if var_pos >= INTRON1_BP_APPROX:
            num_var_sites_before_intron1bp += 1
        if var_pos <= EXON9_BP_APPROX:
            num_var_sites_after_exon9bp += 1

    if cnvtag not in CNVTAG_LOOKUP_TABLE:
        return None
    cn_pattern = CNVTAG_LOOKUP_TABLE[cnvtag]
    cn_list = []
    for _ in range(num_var_sites_after_exon9bp):
        cn_list.append(cn_pattern.exon9_and_downstream)
    for _ in range(
        num_var_sites - num_var_sites_before_intron1bp - num_var_sites_after_exon9bp
    ):
        cn_list.append(cn_pattern.exon9_to_intron1)
    for _ in range(num_var_sites_before_intron1bp):
        cn_list.append(cn_pattern.intron1_upstream)
    return cn_list


def call_cn_snp(total_cn, lsnp1, lsnp2, threshold=0.6):
    """
    Call CN for SNP sites between CYP2D6 and CYP2D7.
    Use a loose cutoff as this is for CNV/hybrid group calling.
    """
    cn_prob = []
    for i, count1 in enumerate(lsnp1):
        count2 = lsnp2[i]
        cn_prob.append(call_reg1_cn(total_cn, count1, count2))
    cn_call_initial = process_raw_call_gc(cn_prob, threshold)
    
    # --- MODIFIED: Add post-processing correction to override statistical model when raw read fractions contradict calls.
    # Critical for targeted sequencing where noise can cause model to call CN=3 when D6 fraction <0.48 indicates CN=2.
    # This raw-data-based override prevents false positive CN3 calls that would propagate through downstream analysis. ---
    corrected_cn_calls = []
    for i, called_cn in enumerate(cn_call_initial):
        count1 = lsnp1[i]
        count2 = lsnp2[i]
        d6_fraction = 0.0
        if (count1 + count2) > 0:
            d6_fraction = count1 / (count1 + count2)
        if called_cn == 3 and d6_fraction < 0.48:
            corrected_cn_calls.append(2)
        else:
            corrected_cn_calls.append(called_cn)            
    return corrected_cn_calls
    # --- END OF MODIFICATION ---


def call_cn_var_homo(total_cn, lsnp1, lsnp2):
    """
    Call CN for variant sites in homology regions.
    """
    cn_prob = []
    for i, count1 in enumerate(lsnp1):
        count2 = lsnp2[i]
        cn_prob.append(call_reg1_cn(total_cn, count1, count2, 4))
    cn_call = []
    for site_call in process_raw_call_denovo(cn_prob, 0.8, 0.65):
        if site_call is None:
            cn_call.append(None)
        else:
            cn_call.append(min(site_call, total_cn - 2))
    return cn_call


def call_cn_var(cnvtag, var_alt, var_ref, alt_forward, alt_reverse, var_list, var_db):
    """
    Call CN for variant sites in non-homology regions.
    Use different minimum read cutoffs for clean variant sites and other sites.
    Total CN at each site is also considered during filtering.
    """
    total_cn_at_site_list = get_total_cn_per_site(cnvtag, var_db, var_list)
    assert total_cn_at_site_list is not None
    cn_prob = []

    # --- MODIFIED: Implement relaxed strand bias filtering for targeted sequencing high-depth data.
    # Original strict P_CUTOFF=0.05 filtered too many true variants in panels with >1000x depth where
    # minor deviations from 50/50 strand distribution become statistically significant despite being biologically real.
    # New logic: (1) Use stricter P_CUTOFF=0.001, (2) Calculate strand ratio, (3) Allow variants with
    # minor strand ≥15% of major strand even if p-value is significant, preventing over-filtering of real variants. ---
    for i, forward_reads in enumerate(alt_forward):
        reverse_reads = alt_reverse[i]
        current_total_ref = var_ref[i]
        current_total_var = var_alt[i]

        if current_total_var > 0 and var_list[i] in NOISY_VAR:
            ntotal_alt = forward_reads + reverse_reads
            table_for_fisher = [[forward_reads, reverse_reads], [ntotal_alt // 2, ntotal_alt - (ntotal_alt // 2)]]
            oddsratio, pvalue = fisher_exact(table_for_fisher)
            
            filter_this_variant_due_to_strand_issue = False
            if forward_reads <= 1 or reverse_reads <= 1:
                filter_this_variant_due_to_strand_issue = True
            elif pvalue < P_CUTOFF:
                if forward_reads == 0 or reverse_reads == 0:
                    strand_ratio = 0.0
                else:
                    strand_ratio = min(forward_reads, reverse_reads) / max(forward_reads, reverse_reads)
                
                if strand_ratio < MIN_ACCEPTABLE_STRAND_RATIO_FOR_NOISY:
                    filter_this_variant_due_to_strand_issue = True
            
            if filter_this_variant_due_to_strand_issue:
                current_total_var = 0

        current_site_expected_cn = total_cn_at_site_list[i]
        if var_list[i] in CLEAN_VAR:
            cn_prob.append(call_reg1_cn(current_site_expected_cn, current_total_var, current_total_ref, 2))
        elif var_list[i] in NOISY_VAR:
            cn_prob.append(call_reg1_cn(current_site_expected_cn, current_total_var, current_total_ref, 7))
        else:
            cn_prob.append(call_reg1_cn(current_site_expected_cn, current_total_var, current_total_ref, 4))
    # --- END OF MODIFICATION ---
    
    cn_call = process_raw_call_denovo(cn_prob, 0.8, 0.65, total_cn_at_site_list)
    return cn_call


def good_read(read):
    """
    Define read filters
    """
    return read.is_secondary == 0 and read.is_supplementary == 0


def get_allele_counts_var42128936(bamfile_handle, genome):
    """
    Search for the inserstions at 42128936 defining
    *30/*40/*58 in read sequences
    """
    long_ins_read = 0
    short_ins_read = 0
    ref_read = 0
    
    # --- MODIFIED: Remove hg19/37 reference coordinates, restricting to hg38 only for CyriPanel. ---
    dregion = {
        # "19": ("chr22", 42524850, 42524980),
        # "37": ("22", 42524850, 42524980),
        38: ("chr22", 42128848, 42128978),
    }
    # --- END OF MODIFICATION ---
    
    region = dregion[genome]
    for read in bamfile_handle.fetch(region[0], region[1], region[2]):
        seq = read.query_sequence
        if good_read(read):
            if "TGGGGCGAAAGGGGCGAAAGGGGCGAAAGGGGCGT" in seq:
                long_ins_read += 1
            elif "TTGGGGCGAAAGGGGCGAAAGGGGCGTC" in seq:
                short_ins_read += 1
            elif "TTGGGGCGAAAGGGGCGTC" in seq:
                ref_read += 1
    return (ref_read, long_ins_read, short_ins_read)


def update_var42128936(
    var_list, var_alt, var_ref, ref_read, long_ins_read, short_ins_read
):
    """
    Update variant read counts for g42128936.
    """
    if "g.42128936-42128937insGGGGCGAAAGGGGCGAAA" in var_list:
        long_ins_index = var_list.index("g.42128936-42128937insGGGGCGAAAGGGGCGAAA")
        var_alt[long_ins_index] = long_ins_read
        var_ref[long_ins_index] = short_ins_read + ref_read
    if "g.42128936-42128937insGGGGCGAAA" in var_list:
        short_ins_index = var_list.index("g.42128936-42128937insGGGGCGAAA")
        var_alt[short_ins_index] = short_ins_read
        var_ref[short_ins_index] = long_ins_read + ref_read
    return var_alt, var_ref


def call_exon9gc(d6_count, d7_count, full_length_cn):
    """
    Call exon 9 conversion
    """
    lsnp1 = d6_count
    lsnp2 = d7_count

    if full_length_cn is not None:
        full_length_cn = int(full_length_cn)
        d6_values = []
        cn_prob = []
        for i, count1 in enumerate(lsnp1):
            count2 = lsnp2[i]
            if (count1 + count2) > 0:
                 d6_values.append(full_length_cn * count1 / (count1 + count2))
            else:
                 d6_values.append(0)
            cn_prob.append(call_reg1_cn(full_length_cn, count1, count2, 3))
        cn_prob_processed_stringent = process_raw_call_gc(cn_prob, 0.88)

        cn_calls = list(set([a for a in cn_prob_processed_stringent if a is not None]))
        if len(cn_calls) == 1:
            count1 = np.mean(d6_count)
            count2 = np.mean(d7_count)
            ave_call = process_raw_call_gc(
                [call_reg1_cn(full_length_cn, count1, count2, 3)], 0.75
            )
            if ave_call[0] is not None and ave_call[0] == cn_calls[0]:
                if ave_call[0] == 1:
                    if d6_values and min(d6_values) < 1.2 and max(d6_values) < 1.3:
                        return ave_call[0]
                else:
                    return ave_call[0]
    return None


# --- MODIFIED: Completely rewrite call_var42126938() for panel-safe, direct pileup-based variant calling.
# Original Cyrius used complex haplotype phasing suitable for WGS but unreliable in targeted panels with irregular coverage.
# New implementation: (1) Direct pileup counting matching samtools, (2) Minimal filtering to avoid over-filtering in high-depth panels,
# (3) Conservative allele fraction thresholds (≥20% for het, ≥75% for hom) with D7 contamination guard (<10%),
# ensuring accurate calls without haplotype dependencies. ---
def call_var42126938(bamfile, full_length_cn, base_db):
    """
    Simple, direct implementation that matches samtools counting
    """
    import pysam
    
    var_called = []
    G_haplotype = False
    
    # Open BAM if path provided
    if isinstance(bamfile, str):
        bam = pysam.AlignmentFile(bamfile)
        should_close = True
    else:
        bam = bamfile
        should_close = False
    
    try:
        # Direct pileup counting - match samtools exactly
        def count_at_position(chrom, pos, ref_base, alt_base):
            ref_count = alt_count = 0
            
            # Use same parameters as samtools depth (minimal filtering)
            for col in bam.pileup(chrom, pos-1, pos, truncate=True,
                                stepper="nofilter",  # Don't filter reads
                                min_base_quality=0,   # Don't filter by base quality
                                ignore_overlaps=False):
                if col.pos != pos-1:
                    continue
                    
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip:
                        continue
                    
                    aln = pr.alignment
                    # Minimal filtering - only skip obvious junk
                    if aln.mapping_quality == 0:
                        continue
                        
                    try:
                        base = aln.query_sequence[pr.query_position].upper()
                        if base == alt_base:
                            alt_count += 1
                        elif base == ref_base:
                            ref_count += 1
                    except (IndexError, TypeError):
                        continue
            
            return ref_count, alt_count
        
        # Count at both positions
        d6_ref, d6_alt = count_at_position("chr22", 42126938, "C", "T")
        d7_ref, d7_alt = count_at_position("chr22", 42140642, "T", "C")        
        d6_d7_base_count = [d6_alt, d7_alt]
        
        # Simple calling logic
        total_d6 = d6_ref + d6_alt
        total_d7 = d7_ref + d7_alt
        
        if total_d6 > 50:  # Require decent depth
            alt_fraction = d6_alt / total_d6
            d7_contamination = (d7_alt / total_d7) if total_d7 > 0 else 0
            
            # Conservative calling - avoid D7 cross-mapping
            if d7_contamination < 0.1 and alt_fraction >= 0.20:
                if alt_fraction >= 0.75:
                    var_called = ["g.42126938C>T", "g.42126938C>T"]  # Homozygous
                else:
                    var_called = ["g.42126938C>T"]  # Heterozygous
    
    finally:
        if should_close:
            bam.close()
    
    return d6_d7_base_count, var_called, G_haplotype
# --- END OF MODIFICATION ---


# --- MODIFIED: Completely rewrite call_var42127526_var42127556() for panel-safe calling without CN-dependence.
# Original relied on CNV tag and haplotype phasing, both unreliable in targeted panels. New CN-neutral implementation:
# (1) Direct pileup counting with flexible chromosome naming, (2) Conservative thresholds: ≥50x depth, ≥5 alt reads,
# ≥25% alt fraction on D6, (3) D7 cross-mapping guard (≤5% alt on D7), (4) No haplotype phasing required,
# enabling accurate calls across all CN configurations without complex dependencies. ---
def call_var42127526_var42127556(
    bamfile,                # str path or pysam.AlignmentFile
    cnvtag,                 # kept for API compatibility; logic below is CN-neutral
    base_db,
    mapq=10,
    baseq=10,
    min_alt_frac=0.25,      # heterozygote threshold on D6
    min_alt_count=5,        # minimal alt reads on D6
    min_total_d6=50,        # minimal total depth on D6
    max_d7_alt_frac=0.05,   # guard for D7 cross-mapping
):
    """
    Panel-safe caller for the paired homology-block sites:
      g.42127526C>T  and  g.42127556T>C   (GRCh38)
    Returns:
      (site42127526_count, site42127556_count, var_called)
      where each *_count is [d6_alt or None, d7_alt or None],
      and var_called is a list like ["g.42127526C>T", "g.42127556T>C"] when supported.
    """
    import pysam

    # ---- local helpers (path or open handle safe) ---------------------------
    def _as_bam(bam_in):
        if hasattr(bam_in, "pileup"):   # AlignmentFile handle
            return bam_in, False
        return pysam.AlignmentFile(bam_in), True

    def _bam_refs(bam):
        try: return set(bam.references)
        except Exception: return set()

    def _norm_chr(contig, bam_in):
        bam, need_close = _as_bam(bam_in)
        try:
            refs = _bam_refs(bam)
            if contig in refs: return contig
            if contig.startswith("chr") and contig[3:] in refs: return contig[3:]
            if ("chr"+contig) in refs: return "chr"+contig
            return contig
        finally:
            if need_close: bam.close()

    def _pileup_one(bam_in, chrom, pos, ref, alt):
        bam, need_close = _as_bam(bam_in)
        ref_n = alt_n = 0
        try:
            for col in bam.pileup(chrom, pos-1, pos, truncate=True,
                                  stepper="all", min_base_quality=baseq,
                                  ignore_overlaps=False):
                if col.pos != pos-1: continue
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip: continue
                    aln = pr.alignment
                    if aln.mapping_quality < mapq: continue
                    b = aln.query_sequence[pr.query_position].upper()
                    if b == alt: alt_n += 1
                    elif b == ref: ref_n += 1
        finally:
            if need_close: bam.close()
        return ref_n, alt_n

    def _find_idx(dsnp, target):
        if isinstance(dsnp, dict):
            if not dsnp: return None
            for k, v in dsnp.items():
                if v == target or (isinstance(v, (list, tuple)) and target in v):
                    return k
            return None
        try: return dsnp.index(target)
        except Exception:
            for i, v in enumerate(dsnp or []):
                if v == target or (isinstance(v, (list, tuple)) and target in v):
                    return i
            return None

    def _coerce_pos(val, prefer):
        if isinstance(val, (list, tuple)):
            return prefer if prefer in val else val[0]
        return val

    # ---- constants & resolve exact rows ------------------------------------
    # GRCh38 canonical coordinates (your haplotype file matches these):
    D6_26, D7_26 = 42127526, 42141231  # C>T on D6; T>C on D7
    D6_56, D7_56 = 42127556, 42141261  # T>C on D6; C>T on D7

    dsnp1 = getattr(base_db, "dsnp1", None)
    dsnp2 = getattr(base_db, "dsnp2", None)

    i26 = _find_idx(dsnp1, D6_26) if dsnp1 is not None else None
    i56 = _find_idx(dsnp1, D6_56) if dsnp1 is not None else None

    if i26 is None or dsnp2 is None:
        d6_26, d7_26 = D6_26, D7_26
    else:
        if isinstance(dsnp1, dict):
            d6_26 = _coerce_pos(dsnp1.get(i26), D6_26)
            d7_26 = _coerce_pos(dsnp2.get(i26), D7_26)
        else:
            d6_26 = _coerce_pos(dsnp1[i26], D6_26)
            d7_26 = _coerce_pos(dsnp2[i26], D7_26)

    if i56 is None or dsnp2 is None:
        d6_56, d7_56 = D6_56, D7_56
    else:
        if isinstance(dsnp1, dict):
            d6_56 = _coerce_pos(dsnp1.get(i56), D6_56)
            d7_56 = _coerce_pos(dsnp2.get(i56), D7_56)
        else:
            d6_56 = _coerce_pos(dsnp1[i56], D6_56)
            d7_56 = _coerce_pos(dsnp2[i56], D7_56)

    nchr = _norm_chr(getattr(base_db, "nchr", "chr22"), bamfile)

    # ---- counts via helper; fallback to direct pileup ----------------------
    try:
        s26_d6, s26_d7 = get_supporting_reads(bamfile, [d6_26], [d7_26], nchr, 0)
    except Exception:
        s26_d6, s26_d7 = [], []
    if not s26_d6 or not s26_d7:
        d6_26_ref, d6_26_alt = _pileup_one(bamfile, nchr, d6_26, "C", "T")
        d7_26_ref, d7_26_alt = _pileup_one(bamfile, nchr, d7_26, "T", "C")
    else:
        d6_26_ref = s26_d6[-2] if len(s26_d6) >= 2 else None
        d6_26_alt = s26_d6[-1] if len(s26_d6) >= 1 else None
        d7_26_ref = s26_d7[-2] if len(s26_d7) >= 2 else None
        d7_26_alt = s26_d7[-1] if len(s26_d7) >= 1 else None

    try:
        s56_d6, s56_d7 = get_supporting_reads(bamfile, [d6_56], [d7_56], nchr, 0)
    except Exception:
        s56_d6, s56_d7 = [], []
    if not s56_d6 or not s56_d7:
        d6_56_ref, d6_56_alt = _pileup_one(bamfile, nchr, d6_56, "T", "C")
        d7_56_ref, d7_56_alt = _pileup_one(bamfile, nchr, d7_56, "C", "T")
    else:
        d6_56_ref = s56_d6[-2] if len(s56_d6) >= 2 else None
        d6_56_alt = s56_d6[-1] if len(s56_d6) >= 1 else None
        d7_56_ref = s56_d7[-2] if len(s56_d7) >= 2 else None
        d7_56_alt = s56_d7[-1] if len(s56_d7) >= 1 else None

    site42127526_count = [d6_26_alt, d7_26_alt]
    site42127556_count = [d6_56_alt, d7_56_alt]

    # ---- conservative panel-safe calling (CN-neutral) ----------------------
    var_called = []

    # 42127526 C>T on D6; T>C on D7
    tot26_d6 = (d6_26_ref or 0) + (d6_26_alt or 0)
    tot26_d7 = (d7_26_ref or 0) + (d7_26_alt or 0)
    d7_26_frac = (d7_26_alt / tot26_d7) if tot26_d7 else 0.0
    if tot26_d6 >= min_total_d6 and (d6_26_alt or 0) >= min_alt_count and d7_26_frac <= max_d7_alt_frac:
        if (d6_26_alt / tot26_d6) >= min_alt_frac:
            var_called.append("g.42127526C>T")

    # 42127556 T>C on D6; C>T on D7
    tot56_d6 = (d6_56_ref or 0) + (d6_56_alt or 0)
    tot56_d7 = (d7_56_ref or 0) + (d7_56_alt or 0)
    d7_56_frac = (d7_56_alt / tot56_d7) if tot56_d7 else 0.0
    if tot56_d6 >= min_total_d6 and (d6_56_alt or 0) >= min_alt_count and d7_56_frac <= max_d7_alt_frac:
        if (d6_56_alt / tot56_d6) >= min_alt_frac:
            var_called.append("g.42127556T>C")

    return site42127526_count, site42127556_count, var_called
# --- END OF MODIFICATION ---


def call_var42127803hap(bamfile, cnvtag, base_db):
    """
    Call haplotype with regard to g.42127803C>T and g.42127941G>A
    """
    diff_haplotype = False
    if cnvtag == "cn2":
        haplotype_per_read = get_haplotypes_from_bam_single_region(
            bamfile, base_db, range(len(base_db.dsnp1))
        )
        recombinant_read_count = extract_hap(haplotype_per_read, [0, 1])
        if (
            "12" in recombinant_read_count
            and sum(recombinant_read_count["12"]) > 1
            and "21" in recombinant_read_count
            and sum(recombinant_read_count["21"]) > 1
        ):
            diff_haplotype = True
    return diff_haplotype


# --- MODIFIED: Add BCyrius function call_var42130655insA() for *15.003 allele identification.
# This insertion is relevant for identifying the *15.003 star allele variant. ---
def call_var42130655insA(bamfile, full_length_cn, base_db):
    """
    Call haplotype with regard to g.42130655-42130656insA (only for calling *15.003)
    """
    var_called = []
    has_long_insert_size_reads_d67 = False
    length_between_d67_regions = 13000 # Distance in bases (at least) between reads in one pair (one is aligned on d6 and the other one on d7)
    no_AA_on_d6 = 0
    no_A_on_d6 = 0

    var_AA_d6, var_AA_d7, var_ref_forward_AA, var_ref_reverse_AA = get_supporting_reads_single_region(bamfile, base_db.dsnp2, base_db.nchr, base_db.dindex, length_between_d67_regions)
    no_AA_on_d6 = var_AA_d7[-1]

    # Set threshold at least 5 read on d7 where the other pair is aligned on d6 and at the same time the number of reads with the insertion on d6 has to be 2 or less
    if int(var_AA_d7[-1]) >= 5 and var_AA_d6[-1] <= 2:
        has_long_insert_size_reads_d67 = True

    if has_long_insert_size_reads_d67:
        # Also count how many AA on d6 (don't use long insert size)
        var_AA_d6, var_AA_d7, var_ref_forward_AA, var_ref_reverse_AA = get_supporting_reads_single_region(bamfile, base_db.dsnp2, base_db.nchr, base_db.dindex)
        no_AA_on_d6 += var_AA_d6[-1]

        var_A_d6, var_A_d7, var_ref_forward_A, var_ref_reverse_A = get_supporting_reads_single_region(bamfile, base_db.dsnp1, base_db.nchr, base_db.dindex)

        no_A_on_d6 = var_A_d6[-1] # This is to determine if heterozygous

        if no_A_on_d6 >= no_AA_on_d6:
            var_called.append("g.42130655-42130656insA")
        else:
            var_called.append("g.42130655-42130656insA")
            var_called.append("g.42130655-42130656insA")

    return [no_A_on_d6, no_AA_on_d6], var_called
# --- END OF MODIFICATION ---


def get_called_variants(var_list, cn_prob_processed, starting_index=0):
    """
    Return called variants based on called copy number and list of variant names
    """
    total_callset = []
    
    # --- MODIFIED: Add bounds checking to prevent index errors when processing variant lists.
    # Ensures starting_index and cn_prob_processed length are compatible with var_list to avoid crashes. ---
    if starting_index > len(var_list) or len(cn_prob_processed) > (len(var_list) - starting_index):
        pass

    for i, cn_called in enumerate(cn_prob_processed):
        current_var_index = i + starting_index
        if current_var_index < len(var_list):
            if cn_called is not None and cn_called != 0:
                for _ in range(cn_called):
                    total_callset.append(var_list[current_var_index])
    return total_callset
    # --- END OF MODIFICATION ---


# --- MODIFIED: Add get_fallback_combinations() to generate FALLBACK_CN_COMBINATIONS for fallback mechanism.
# Extracts all valid (total_cn, d7_spacer) pairs from CNVTAG_LOOKUP_TABLE, creating a comprehensive list
# of alternative CN configurations to test when initial CNVPanelizer values fail. Sorted by total_cn for
# systematic exploration starting from simpler configurations, enabling robust genotyping in ambiguous cases. ---
def get_fallback_combinations():
    """
    Generates a list of (total_cn, d7_spacer) tuples from the CNVTAG_LOOKUP_TABLE
    to be used as fallback options for genotyping.
    """
    fallback_combinations = set()
    for cnv_tag, cn_profile in CNVTAG_LOOKUP_TABLE.items():
        total_cn = cn_profile.total_cn
        # exon9_cn is derived from total_cn - spacer_cn
        # So, spacer_cn = total_cn - exon9_and_downstream_cn
        d7_spacer = total_cn - cn_profile.exon9_and_downstream
        if d7_spacer >= 0:
            fallback_combinations.add((total_cn, d7_spacer))
        
    return sorted(list(fallback_combinations)) # Sort the list to have a consistent order, starting with lower total_cn

FALLBACK_CN_COMBINATIONS = get_fallback_combinations()
# --- END OF MODIFICATION ---
