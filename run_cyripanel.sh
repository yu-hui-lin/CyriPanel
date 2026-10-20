#!/bin/bash
#
# CyriPanel general run wrapper
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

set -euo pipefail

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
    exit "${1:-1}"
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
            usage 0
            ;;
        -o|--outdir)
            [[ $# -ge 2 && -n "$2" && "$2" != -* ]] || { echo "ERROR: --outdir requires a directory." >&2; exit 2; }
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --manifest)
            [[ $# -ge 2 && -n "$2" && "$2" != -* ]] || { echo "ERROR: --manifest requires a file." >&2; exit 2; }
            [[ -z "$MANIFEST_FILE" && -z "$INPUT_BAM" ]] || { echo "ERROR: Specify one BAM or one manifest." >&2; exit 2; }
            MANIFEST_MODE=true
            MANIFEST_FILE="$2"
            shift 2
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            exit 2
            ;;
        *)
            if [[ -z "${INPUT_BAM}" ]] && [[ "${MANIFEST_MODE}" == "false" ]]; then
                INPUT_BAM="$1"
            elif [[ -z "${OUTPUT_PREFIX}" ]]; then
                OUTPUT_PREFIX="$1"
            else
                echo "ERROR: Too many arguments." >&2
                exit 2
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

# Preserve BAM symlink paths so adjacent indexes remain selectable.
for value in "$INPUT_BAM" "$MANIFEST_FILE" "$OUTPUT_DIR" "$OUTPUT_PREFIX"; do
    [[ "$value" != *$'\n'* && "$value" != *$'\r'* ]] || {
        echo "ERROR: Newlines are not supported in paths or prefixes." >&2; exit 2;
    }
done
[[ -n "$OUTPUT_PREFIX" && "$OUTPUT_PREFIX" != */* &&
   "$OUTPUT_PREFIX" != "." && "$OUTPUT_PREFIX" != ".." ]] || {
    echo "ERROR: Prefix must be a filename, without directory components." >&2; exit 2;
}
absolute_path() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *) printf '%s/%s\n' "$PWD" "$1" ;;
    esac
}
OUTPUT_DIR=$(absolute_path "$OUTPUT_DIR")
if [[ "$MANIFEST_MODE" == "true" ]]; then
    MANIFEST_FILE=$(absolute_path "$MANIFEST_FILE")
else
    INPUT_BAM=$(absolute_path "$INPUT_BAM")
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
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/${OUTPUT_PREFIX}_${TIMESTAMP}.log"

# Create temporary manifest if single BAM mode
if [[ "${MANIFEST_MODE}" == "false" ]]; then
    TEMP_MANIFEST=$(mktemp)
    trap 'rm -f -- "$TEMP_MANIFEST"' EXIT
    printf '%s\n' "$INPUT_BAM" > "$TEMP_MANIFEST"
    MANIFEST_FILE="$TEMP_MANIFEST"
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
    --prefix="${OUTPUT_PREFIX}" \
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
