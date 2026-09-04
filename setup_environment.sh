#!/bin/bash
#
# CyriPanel Environment Setup Script for NCHC
#
# This script creates the conda environment and installs all dependencies
# for running CyriPanel on NCHC.
#
# Usage: bash setup_environment.sh
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="cyripanel"

echo "=========================================="
echo "CyriPanel Environment Setup"
echo "=========================================="
echo ""

# =============================================================================
# Step 1: Load Miniconda
# =============================================================================
log_info "Step 1: Loading Miniconda module..."

# Try different module names commonly used on NCHC
if module load miniconda3 2>/dev/null; then
    log_info "Loaded miniconda3 module"
elif module load anaconda3 2>/dev/null; then
    log_info "Loaded anaconda3 module"
elif module load conda 2>/dev/null; then
    log_info "Loaded conda module"
else
    log_warn "Could not load conda module. Checking if conda is available..."
    if ! command -v conda &> /dev/null; then
        log_error "Conda not found. Please install Miniconda or load the appropriate module."
        exit 1
    fi
fi

# Initialize conda for this shell
eval "$(conda shell.bash hook)"

# =============================================================================
# Step 2: Create Conda Environment
# =============================================================================
log_info "Step 2: Creating conda environment '${ENV_NAME}'..."

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    log_warn "Environment '${ENV_NAME}' already exists."
    read -p "Do you want to remove and recreate it? (y/N): " REPLY
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n ${ENV_NAME} -y
    else
        log_info "Keeping existing environment. Skipping creation."
        conda activate ${ENV_NAME}
    fi
fi

if ! conda env list | grep -q "^${ENV_NAME} "; then
    conda create -n ${ENV_NAME} python=3.10 -y
fi

conda activate ${ENV_NAME}
log_info "Activated environment: ${ENV_NAME}"

# =============================================================================
# Step 3: Install Python Dependencies
# =============================================================================
log_info "Step 3: Installing Python dependencies..."

pip install --upgrade pip
pip install numpy>=1.16 scipy>=1.2 pysam>=0.15.3 statsmodels>=0.9 pandas>=1.0.0

log_info "Python dependencies installed successfully"

# Verify Python packages
log_info "Verifying Python packages..."
python3 -c "
import numpy as np
import scipy
import pysam
import statsmodels
import pandas as pd
print(f'  numpy:       {np.__version__}')
print(f'  scipy:       {scipy.__version__}')
print(f'  pysam:       {pysam.__version__}')
print(f'  statsmodels: {statsmodels.__version__}')
print(f'  pandas:      {pd.__version__}')
"

# =============================================================================
# Step 4: Install R and CNVPanelizer
# =============================================================================
log_info "Step 4: Installing R and CNVPanelizer..."

# Install R through conda
conda install -c conda-forge r-base=4.2 -y

log_info "Installing R packages (this may take several minutes)..."

# Set R library path for user packages
R_LIB_PATH="${HOME}/R/library"
mkdir -p "${R_LIB_PATH}"

# Install BiocManager and CNVPanelizer
Rscript -e "
# Set library path
.libPaths(c('${R_LIB_PATH}', .libPaths()))

# Install BiocManager
if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager', repos='https://cloud.r-project.org', lib='${R_LIB_PATH}')
}

# Install CNVPanelizer
BiocManager::install('CNVPanelizer', update=FALSE, ask=FALSE, lib='${R_LIB_PATH}')

# Verify installation
library(CNVPanelizer)
cat('CNVPanelizer installed successfully!\n')
"

log_info "R and CNVPanelizer installed successfully"

# =============================================================================
# Step 5: Verify Complete Installation
# =============================================================================
log_info "Step 5: Verifying complete installation..."

echo ""
echo "=========================================="
echo "Installation Verification"
echo "=========================================="

# Check Python
echo ""
echo "Python Environment:"
which python3
python3 --version

# Check R
echo ""
echo "R Environment:"
which Rscript
Rscript --version

# Check CNVPanelizer
echo ""
echo "CNVPanelizer Check:"
Rscript -e "library(CNVPanelizer); cat('CNVPanelizer version:', as.character(packageVersion('CNVPanelizer')), '\n')"

# =============================================================================
# Create activation script
# =============================================================================
ACTIVATE_SCRIPT="${SCRIPT_DIR}/activate_cyripanel.sh"
cat > "${ACTIVATE_SCRIPT}" << 'EOF'
#!/bin/bash
# Activate CyriPanel environment
# Usage: source activate_cyripanel.sh

# Load conda module (adjust as needed for your NCHC setup)
module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || module load conda 2>/dev/null

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate environment
conda activate cyripanel

# Set R library path
export R_LIBS_USER="${HOME}/R/library"

echo "CyriPanel environment activated!"
echo "Python: $(which python3)"
echo "R: $(which Rscript)"
EOF
chmod +x "${ACTIVATE_SCRIPT}"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================="
echo -e "${GREEN}Installation Complete!${NC}"
echo "=========================================="
echo ""
echo "To use CyriPanel in future sessions:"
echo ""
echo "  1. Source the activation script:"
echo "     source ${ACTIVATE_SCRIPT}"
echo ""
echo "  2. Or manually activate:"
echo "     module load miniconda3"
echo "     conda activate ${ENV_NAME}"
echo ""
echo "  3. Run CyriPanel:"
echo "     bash run_cyripanel.sh <input.bam>"
echo ""
echo "  4. Run test to verify installation:"
echo "     bash run_test_HG02129.sh"
echo ""
