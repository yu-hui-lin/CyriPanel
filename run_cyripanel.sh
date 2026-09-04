#!/bin/bash
#
# CyriPanel Production Run Script
# CYP2D6 Genotyper for Targeted Sequencing Panels
#
# Usage: 
#   bash run_cyripanel.sh <input_bam> [output_prefix]
#   bash run_cyripanel.sh --manifest <manifest_file> [output_prefix]
#
# Examples:
#   bash run_cyripanel.sh /path/to/sample.bam
#   bash run_cyripanel.sh /path/to/sample.bam sample001
#   bash run_cyripanel.sh --manifest samples.manifest batch001
#

set -e

# =============================================================================
# Configuration
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CYRIPANEL_DIR="${SCRIPT_DIR}"

# Default output directory
DEFAULT_OUTPUT_DIR="${CYRIPANEL_DIR}/results"

# =============================================================================
# Usage
# =============================================================================
usage() {
    echo "CyriPanel - CYP2D6 Genotyper for Targeted Sequencing Panels"
    echo ""
    echo "Usage:"
    echo "  $(basename "$0") <input_bam> [output_prefix]"
    echo "  $(basename "$0") --manifest <manifest_file> [output_prefix]"
    echo ""
    echo "Arguments:"
    echo "  <input_bam>       Path to input BAM file"
    echo "  --manifest        Use manifest file listing multiple BAM paths"
    echo "  [output_prefix]   Optional prefix for output files (default: sample name)"
    echo ""
    echo "Options:"
    echo "  -h, --help        Show this help message"
    echo "  -o, --outdir DIR  Output directory (default: ${DEFAULT_OUTPUT_DIR})"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") /data/sample001.bam"
    echo "  $(basename "$0") /data/sample001.bam my_sample"
    echo "  $(basename "$0") --manifest samples.txt batch_analysis"
    echo "  $(basename "$0") -o /custom/output /data/sample001.bam"
    echo ""
    exit 1
}

# =============================================================================
# Parse Arguments
# =============================================================================
OUTPUT_DIR="${DEFAULT_OUTPUT_DIR}"
MANIFEST_MODE=false
MANIFEST_FILE=""
INPUT_BAM=""
OUTPUT_PREFIX=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -o|--outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --manifest)
            MANIFEST_MODE=true
            MANIFEST_FILE="$2"
            shift 2
            ;;
        *)
            if [[ -z "${INPUT_BAM}" ]] && [[ "${MANIFEST_MODE}" == "false" ]]; then
                INPUT_BAM="$1"
            elif [[ -z "${OUTPUT_PREFIX}" ]]; then
                OUTPUT_PREFIX="$1"
            fi
            shift
            ;;
    esac
done

# =============================================================================
# Validation
# =============================================================================
if [[ "${MANIFEST_MODE}" == "true" ]]; then
    if [[ -z "${MANIFEST_FILE}" ]] || [[ ! -f "${MANIFEST_FILE}" ]]; then
        echo "ERROR: Manifest file not found: ${MANIFEST_FILE}"
        exit 1
    fi
    # Set default prefix from manifest filename
    if [[ -z "${OUTPUT_PREFIX}" ]]; then
        OUTPUT_PREFIX=$(basename "${MANIFEST_FILE}" | sed 's/\.[^.]*$//')
    fi
else
    if [[ -z "${INPUT_BAM}" ]]; then
        usage
    fi
    if [[ ! -f "${INPUT_BAM}" ]]; then
        echo "ERROR: Input BAM file not found: ${INPUT_BAM}"
        exit 1
    fi
    # Set default prefix from BAM filename
    if [[ -z "${OUTPUT_PREFIX}" ]]; then
        OUTPUT_PREFIX=$(basename "${INPUT_BAM}" | sed 's/\.[^.]*$//')
    fi
fi

# Check star_caller.py exists
if [[ ! -f "${CYRIPANEL_DIR}/star_caller.py" ]]; then
    echo "ERROR: star_caller.py not found in ${CYRIPANEL_DIR}"
    exit 1
fi

# =============================================================================
# Setup
# =============================================================================
mkdir -p "${OUTPUT_DIR}"
LOG_DIR="${CYRIPANEL_DIR}/logs"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/${OUTPUT_PREFIX}_${TIMESTAMP}.log"

# Create temporary manifest if single BAM mode
if [[ "${MANIFEST_MODE}" == "false" ]]; then
    TEMP_MANIFEST=$(mktemp)
    echo "${INPUT_BAM}" > "${TEMP_MANIFEST}"
    MANIFEST_FILE="${TEMP_MANIFEST}"
    trap "rm -f ${TEMP_MANIFEST}" EXIT
fi

# =============================================================================
# Run CyriPanel
# =============================================================================
echo "=========================================="
echo "CyriPanel - CYP2D6 Genotyper"
echo "=========================================="
echo "Start Time:    $(date)"
echo "Input:         ${MANIFEST_FILE}"
echo "Output Dir:    ${OUTPUT_DIR}"
echo "Output Prefix: ${OUTPUT_PREFIX}"
echo "Log File:      ${LOG_FILE}"
echo "=========================================="
echo ""

python3 "${CYRIPANEL_DIR}/star_caller.py" \
    --manifest "${MANIFEST_FILE}" \
    --genome 38 \
    --outDir "${OUTPUT_DIR}" \
    --prefix "${OUTPUT_PREFIX}" \
    2>&1 | tee "${LOG_FILE}"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================="
echo "Analysis Complete"
echo "=========================================="
echo "End Time: $(date)"
echo ""
echo "Output files:"
echo "  TSV:  ${OUTPUT_DIR}/${OUTPUT_PREFIX}.tsv"
echo "  JSON: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.json"
echo "  Log:  ${LOG_FILE}"
echo ""

# Display results
if [[ -f "${OUTPUT_DIR}/${OUTPUT_PREFIX}.tsv" ]]; then
    echo "Results:"
    cat "${OUTPUT_DIR}/${OUTPUT_PREFIX}.tsv"
fi
