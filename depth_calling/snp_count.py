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


from collections import namedtuple
import pysam


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def reverse_complement(sequence):
    """Return the reverse complement of a sequence."""
    return "".join(COMPLEMENT[b] for b in sequence[::-1])


def get_nm(ltag):
    """Return the value of the NM tag."""
    for tag in ltag:
        if tag[0] == "NM":
            return tag[1]
    return None


def get_snp_position(pos_file, genome, group=None):
    """Get all base differences listed in the SNP location file."""
    dsnp1 = {}
    dsnp2 = {}
    dindex = {}
    nchr = None

    with open(pos_file) as read_pos:
        counter = -1
        for line in read_pos:
            if line[0] != "#" and line[0] != "\n":
                split_line = line.strip().split()

                # --- MODIFIED: Extract chromosome name from file data rather than using genome parameter.
                # Original assigned nchr from split_line[0] only after loop completion. 
                # New approach assigns it on first valid line, ensuring nchr is available earlier and preserving 'chr' prefix consistency
                # (e.g., 'chr22' instead of '22') for better BAM compatibility. ---
                if nchr is None:
                    nchr = split_line[0]
                # --- END OF MODIFICATION ---

                if group is None or split_line[-1] == group:
                    counter += 1
                    reg1_name = split_line[1] + "_" + str(counter)
                    reg2_name = split_line[3] + "_" + str(counter)
                    reg1_base = split_line[2].upper()
                    reg2_base = split_line[4].upper()
                    if split_line[-1] != "-":
                        dsnp1.setdefault(reg1_name, "_".join([reg1_base, reg2_base]))
                        dsnp2.setdefault(reg2_name, "_".join([reg1_base, reg2_base]))
                    else:
                        dsnp1.setdefault(
                            reg1_name,
                            "_".join([reg1_base, reverse_complement(reg2_base)]),
                        )
                        dsnp2.setdefault(
                            reg2_name,
                            "_".join([reverse_complement(reg1_base), reg2_base]),
                        )
                    dindex.setdefault(reg1_name, counter)
                    dindex.setdefault(reg2_name, counter)

    snp_lookup = namedtuple("snp_lookup", "dsnp1 dsnp2 nchr dindex")
    dbsnp = snp_lookup(dsnp1, dsnp2, nchr, dindex)
    return dbsnp


def passing_read(pileupread):
    """Return whether a read passes filter."""
    return (
        not pileupread.is_del
        and not pileupread.is_refskip
        and pileupread.alignment.is_secondary == 0
        and pileupread.alignment.is_supplementary == 0
        and pileupread.alignment.is_duplicate == 0
    )


def passing_read_stringent(pileupread):
    """Return whether a read passes more stringent filter."""
    number_mismatch = get_nm(pileupread.alignment.tags)
    align_len = pileupread.alignment.query_alignment_length
    read_len = len(pileupread.alignment.query_sequence)
    insert_size = pileupread.alignment.template_length
    # number_mismatch <= float(align_len) * 0.08
    # and pileupread.query_position > 0
    # and pileupread.query_position < read_len - 1
    return abs(insert_size) < 1000


def get_reads_by_region(
    bamfile_handle, nchr, dsnp, dindex, min_mapq=0, min_insert_size=None, stringent=False
):
    """
    Return the number of reads supporting region1 and region2, forward and reverse.
    """
    lsnp1_forward = []
    lsnp1_reverse = []
    lsnp2_forward = []
    lsnp2_reverse = []
    for _ in dsnp:
        lsnp1_forward.append(set())
        lsnp1_reverse.append(set())
        lsnp2_forward.append(set())
        lsnp2_reverse.append(set())

    for snp_position_ori in dsnp:
        snp_position = int(snp_position_ori.split("_")[0])
        for pileupcolumn in bamfile_handle.pileup(
            nchr,
            snp_position - 1,
            snp_position + 1,
            truncate=True,
            stepper="nofilter",
            ignore_overlaps=False,
            ignore_orphan=False,
        ):
            site_position = pileupcolumn.pos + 1
            if site_position == snp_position:
                reg1_allele, reg2_allele = dsnp[snp_position_ori].split("_")
                for read in pileupcolumn.pileups:
                    
                    # --- MODIFIED: Add min_insert_size parameter for insert size filtering following BCyrius update.
                    # Enables filtering reads by template length, useful for distinguishing between CYP2D6 and CYP2D7 in specific scenarios 
                    # where read pairs spanning large distances indicate cross-paralog mapping. Used specifically in call_var42130655insA() for *15.003 detection. ---
                    if passing_read(read) and read.alignment.mapping_quality >= min_mapq:
                        if min_insert_size is not None and abs(read.alignment.template_length) < min_insert_size:
                            continue
                    # --- END OF MODIFICATION ---
                        
                        dsnp_index = dindex[snp_position_ori]
                        read_name = read.alignment.query_name
                        read_seq = read.alignment.query_sequence
                        if stringent is False or passing_read_stringent(read):
                            reg1_allele_split = reg1_allele.split(",")
                            reg2_allele_split = reg2_allele.split(",")
                            start_pos = read.query_position
                            for allele in reg1_allele_split:
                                end_pos = start_pos + len(allele)
                                if read_seq[start_pos:end_pos] == allele:
                                    if read.alignment.is_reverse:
                                        lsnp1_reverse[dsnp_index].add(read_name)
                                    else:
                                        lsnp1_forward[dsnp_index].add(read_name)
                            for allele in reg2_allele_split:
                                end_pos = start_pos + len(allele)
                                if read_seq[start_pos:end_pos] == allele:
                                    if read.alignment.is_reverse:
                                        lsnp2_reverse[dsnp_index].add(read_name)
                                    else:
                                        lsnp2_forward[dsnp_index].add(read_name)
    return lsnp1_forward, lsnp1_reverse, lsnp2_forward, lsnp2_reverse


def get_fraction(lsnp1, lsnp2):
    """Return the fraction of reads supporting region1."""
    reg1_fraction = []
    for index in range(len(lsnp1)):
        sumdepth = lsnp1[index] + lsnp2[index]
        if sumdepth == 0:
            reg1_fraction.append(0)
        else:
            reg1_fraction.append(float(lsnp1[index]) / float(sumdepth))
    return reg1_fraction


def merge_reads(list_to_merge):
    """
    Merge sets of reads (forward/reverse, region1/region2)
    """
    merged_reads = []
    for i in range(len(list_to_merge[0])):
        reads_per_site = set()
        for lsnp in list_to_merge:
            reads_per_site = reads_per_site.union(lsnp[i])
        merged_reads.append(reads_per_site)
    return merged_reads


def get_supporting_reads(bamfile_handle, dsnp1, dsnp2, nchr, dindex):
    """
    Return the number of supporting reads at each position in
    both region1 and region2.
    """
    assert len(dsnp1) == len(dsnp2)
    # Go through SNP sites in both regions,
    # and count the number of reads supporting each gene.
    lsnp1_reg1_for, lsnp1_reg1_rev, lsnp2_reg1_for, lsnp2_reg1_rev = get_reads_by_region(
        bamfile_handle, nchr, dsnp1, dindex
    )
    lsnp1_reg2_for, lsnp1_reg2_rev, lsnp2_reg2_for, lsnp2_reg2_rev = get_reads_by_region(
        bamfile_handle, nchr, dsnp2, dindex
    )
    lsnp1 = merge_reads(
        [lsnp1_reg1_for, lsnp1_reg1_rev, lsnp1_reg2_for, lsnp1_reg2_rev]
    )
    lsnp2 = merge_reads(
        [lsnp2_reg1_for, lsnp2_reg1_rev, lsnp2_reg2_for, lsnp2_reg2_rev]
    )
    return [len(a) for a in lsnp1], [len(a) for a in lsnp2]


def get_supporting_reads_single_region(bamfile_handle, dsnp1, nchr, dindex, min_ins_size=None):
    """
    Return the number of supporting reads at each position only in region1, as well 
    as the number of alt reads in forward and reverse.
    """
    
    # --- MODIFIED: Add min_ins_size parameter to pass through to get_reads_by_region() following BCyrius update.
    # Enables insert size filtering for single-region read counting, specifically used in call_var42130655insA() to detect *15.003 allele 
    # by identifying read pairs with one mate in D6 and the other in D7 (spanning ≥13kb), indicating the insertion variant. ---
    lsnp1_for, lsnp1_rev, lsnp2_for, lsnp2_rev = get_reads_by_region(
        bamfile_handle, nchr, dsnp1, dindex, 10, min_ins_size
    )
    # --- END OF MODIFICATION ---
    
    lsnp1 = merge_reads([lsnp1_for, lsnp1_rev])
    lsnp2 = merge_reads([lsnp2_for, lsnp2_rev])
    return (
        [len(a) for a in lsnp1],
        [len(a) for a in lsnp2],
        [len(a) for a in lsnp1_for],
        [len(a) for a in lsnp1_rev],
    )
