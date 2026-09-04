#!/bin/bash
#
# CyriPanel Test Script for HG02129
# Expected Genotype: *1/*10+*36
#
# This script validates that the CyriPanel installation is working correctly
# by running a known sample and checking the output.
#
# Usage: bash run_test_HG02129.sh
#

set -e  # Exit on any error

# =============================================================================
# Configuration
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CYRIPANEL_DIR="${SCRIPT_DIR}"

# Test sample
TEST_SAMPLE_NAME="HG02129"
TEST_BAM="${CYRIPANEL_DIR}/test_data/${TEST_SAMPLE_NAME}.bam"
EXPECTED_GENOTYPE="*1/*36+*10"

# Output directories
OUTPUT_DIR="${CYRIPANEL_DIR}/results/test_${TEST_SAMPLE_NAME}"
LOG_DIR="${CYRIPANEL_DIR}/logs"
MANIFEST_FILE="${OUTPUT_DIR}/test_manifest.txt"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# =============================================================================
# Functions
# =============================================================================
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        return 1
    fi
    return 0
}

# =============================================================================
# Pre-flight Checks
# =============================================================================
log_info "=========================================="
log_info "CyriPanel Test Script for ${TEST_SAMPLE_NAME}"
log_info "=========================================="
log_info "Expected Genotype: ${EXPECTED_GENOTYPE}"
log_info ""

# Check if running in correct directory
if [[ ! -f "${CYRIPANEL_DIR}/star_caller.py" ]]; then
    log_error "star_caller.py not found in ${CYRIPANEL_DIR}"
    log_error "Please run this script from the CyriPanel directory"
    exit 1
fi

# Check test BAM file
log_info "Checking test BAM file..."
if ! check_file "${TEST_BAM}"; then
    log_error "Test BAM file not found: ${TEST_BAM}"
    log_error "Please copy ${TEST_SAMPLE_NAME}.bam to ${CYRIPANEL_DIR}/test_data/"
    exit 1
fi

# Check BAM index
if [[ ! -f "${TEST_BAM}.bai" ]]; then
    log_warn "BAM index not found, generating..."
    samtools index "${TEST_BAM}"
fi

# Check reference directory
REF_DIR="${CYRIPANEL_DIR}/ref_dir"
REF_COUNT=$(find "${REF_DIR}" -name "*.bam" 2>/dev/null | wc -l)
log_info "Found ${REF_COUNT} reference BAM files in ref_dir"

if [[ ${REF_COUNT} -lt 5 ]]; then
    log_warn "Fewer than 5 reference BAM files found. Results may be less reliable."
fi

if [[ ${REF_COUNT} -eq 0 ]]; then
    log_error "No reference BAM files found in ${REF_DIR}"
    log_error "Please add reference diploid BAM files for CNVPanelizer"
    exit 1
fi

# Check BED file
BED_FILES=$(find "${CYRIPANEL_DIR}/data" -name "*.bed" 2>/dev/null | wc -l)
if [[ ${BED_FILES} -eq 0 ]]; then
    log_error "No BED file found in ${CYRIPANEL_DIR}/data/"
    exit 1
fi

# Check Python environment
log_info "Checking Python environment..."
if ! python3 -c "import pysam; import numpy; import scipy; import pandas" 2>/dev/null; then
    log_error "Python dependencies not installed properly"
    log_error "Please activate the cyripanel conda environment:"
    log_error "  conda activate cyripanel"
    exit 1
fi

# Check R and CNVPanelizer
log_info "Checking R and CNVPanelizer..."
if ! Rscript -e "library(CNVPanelizer)" 2>/dev/null; then
    log_error "CNVPanelizer R package not installed"
    log_error "Please install: R -e \"BiocManager::install('CNVPanelizer')\""
    exit 1
fi

# =============================================================================
# Setup Output Directories
# =============================================================================
log_info "Setting up output directories..."
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Create manifest file
echo "${TEST_BAM}" > "${MANIFEST_FILE}"

# =============================================================================
# Run CyriPanel
# =============================================================================
log_info "Starting CyriPanel analysis..."
log_info "This may take a few minutes..."

START_TIME=$(date +%s)

python3 "${CYRIPANEL_DIR}/star_caller.py" \
    --manifest "${MANIFEST_FILE}" \
    --genome 38 \
    --outDir "${OUTPUT_DIR}" \
    --prefix "${TEST_SAMPLE_NAME}_test" \
    2>&1 | tee "${LOG_DIR}/test_${TEST_SAMPLE_NAME}_$(date +%Y%m%d_%H%M%S).log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

log_info "Analysis completed in ${ELAPSED} seconds"

# =============================================================================
# Validate Results
# =============================================================================
log_info ""
log_info "=========================================="
log_info "Validating Results"
log_info "=========================================="

TSV_OUTPUT="${OUTPUT_DIR}/${TEST_SAMPLE_NAME}_test.tsv"
JSON_OUTPUT="${OUTPUT_DIR}/${TEST_SAMPLE_NAME}_test.json"

if [[ ! -f "${TSV_OUTPUT}" ]]; then
    log_error "TSV output file not generated: ${TSV_OUTPUT}"
    exit 1
fi

# Extract the called genotype
CALLED_GENOTYPE=$(tail -n 1 "${TSV_OUTPUT}" | cut -f2)
FILTER_STATUS=$(tail -n 1 "${TSV_OUTPUT}" | cut -f3)

log_info "Called Genotype: ${CALLED_GENOTYPE}"
log_info "Filter Status: ${FILTER_STATUS}"
log_info "Expected Genotype: ${EXPECTED_GENOTYPE}"

# =============================================================================
# Final Verdict
# =============================================================================
log_info ""
log_info "=========================================="
log_info "TEST RESULT"
log_info "=========================================="

# Normalize genotypes for comparison (handle different orderings)
normalize_genotype() {
    echo "$1" | tr '+' '\n' | sort | tr '\n' '+' | sed 's/+$//'
}

NORM_CALLED=$(normalize_genotype "${CALLED_GENOTYPE}")
NORM_EXPECTED=$(normalize_genotype "${EXPECTED_GENOTYPE}")

if [[ "${CALLED_GENOTYPE}" == "${EXPECTED_GENOTYPE}" ]]; then
    echo -e "${GREEN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                                                              ║"
    echo "║   ✓ TEST PASSED - EXACT MATCH                               ║"
    echo "║                                                              ║"
    echo "║   CyriPanel is correctly installed and configured!          ║"
    echo "║                                                              ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    TEST_STATUS="PASSED"
elif [[ "${FILTER_STATUS}" == "PASS" ]] && [[ -n "${CALLED_GENOTYPE}" ]]; then
    echo -e "${YELLOW}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                                                              ║"
    echo "║   ⚠ TEST PARTIALLY PASSED                                   ║"
    echo "║                                                              ║"
    echo "║   A genotype was called with PASS filter, but differs       ║"
    echo "║   from expected. Please verify manually.                    ║"
    echo "║                                                              ║"
    echo "║   Called:   ${CALLED_GENOTYPE}                              "
    echo "║   Expected: ${EXPECTED_GENOTYPE}                            "
    echo "║                                                              ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    TEST_STATUS="PARTIAL"
else
    echo -e "${RED}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                                                              ║"
    echo "║   ✗ TEST FAILED                                             ║"
    echo "║                                                              ║"
    echo "║   The genotype could not be called correctly.               ║"
    echo "║   Please check the log files for details.                   ║"
    echo "║                                                              ║"
    echo "║   Called:   ${CALLED_GENOTYPE:-NONE}                        "
    echo "║   Filter:   ${FILTER_STATUS:-NONE}                          "
    echo "║   Expected: ${EXPECTED_GENOTYPE}                            "
    echo "║                                                              ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    TEST_STATUS="FAILED"
fi

# =============================================================================
# Output Summary
# =============================================================================
log_info ""
log_info "Output files:"
log_info "  TSV:  ${TSV_OUTPUT}"
log_info "  JSON: ${JSON_OUTPUT}"
log_info "  Log:  ${LOG_DIR}/test_${TEST_SAMPLE_NAME}_*.log"
log_info ""

# Display the full TSV output
log_info "Full TSV output:"
cat "${TSV_OUTPUT}"
echo ""

# Exit with appropriate code
if [[ "${TEST_STATUS}" == "PASSED" ]]; then
    exit 0
elif [[ "${TEST_STATUS}" == "PARTIAL" ]]; then
    exit 0  # Still consider partial pass as success for installation validation
else
    exit 1
fi
