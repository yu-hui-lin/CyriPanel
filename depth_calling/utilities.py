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


def parse_region_file(region_file, genome):
    """Return the set of regions for counting from a bed file."""
    region_dic = {}
    with open(region_file) as read_region:
        for line in read_region:
            # Skip any header or malformed lines
            if line.startswith("#") or len(line.strip().split()) != 4:
                continue

            # --- MODIFIED: Adapt BED file parsing for standard 4-column format instead of 6-column Cyrius format.
            # Original Cyrius expected: chr, start, end, name, type, GC_content (6 columns)
            # CyriPanel uses standard BED: chr, start, end, name (4 columns)
            # Modification unpacks only 4 columns and assigns default placeholder values for missing region_type
            # ("target") and region_gc ("0.0"), enabling compatibility with standard BED files from panel designs
            # while maintaining downstream code structure that expects these fields. ---
            nchr, region_start, region_end, region_name = line.strip().split()
            region_type = "target"  # A generic type for all regions
            region_gc = "0.0"       # A default GC content value
            # --- END OF MODIFICATION ---

            region_start = int(region_start)
            region_end = int(region_end)
            region = (nchr, region_start, region_end, region_name)
            
            # The rest of the logic remains the same
            region_dic.setdefault(region_type, []).append((region, region_gc))            
    return region_dic


# --- MODIFIED: Remove parse_gmm_file() function (not used in CyriPanel).
# Original Cyrius relied on Gaussian Mixture Model (GMM) parameters stored in separate files for CNV calling. 
# CyriPanel replaces GMM-based CNV detection with CNVPanelizer integration. ---
# def parse_gmm_file(gmm_file):
#     """Return the gmm parameters stored in input file."""
#     dpar_tmp = {}
#     with open(gmm_file) as read_gmm:
#         for line in read_gmm:
#             split_line = line.strip().split()
#             dpar_tmp.setdefault(split_line[0], {})
#             list_value = [a.split(":")[-1] for a in split_line[2:]]
#             dpar_tmp[split_line[0]].setdefault(split_line[1], list_value)
#     return dpar_tmp
# --- END OF MODIFICATION ---


def open_alignment_file(alignment_file, reference_fasta=None, index_filename=None):
    if alignment_file.endswith("cram"):
        return pysam.AlignmentFile(
            alignment_file, "rc", reference_filename=reference_fasta, index_filename=index_filename
        )
    return pysam.AlignmentFile(alignment_file, "rb", index_filename=index_filename)
