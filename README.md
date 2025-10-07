# CyriPanel
CYP2D6 genotyper for targeted panels with modified Cyrius and CNVPanelizer

## Introduction
Cyrius, an open-source tool by Illumina, excels at CYP2D6 genotyping from whole-genome sequencing (WGS) data. **However, its reliance on a Gaussian Mixture Model (GMM) trained on WGS coverage patterns from 2,504 samples in the 1000 Genomes Project creates significant limitations for laboratory-developed targeted sequencing panels, whose coverage patterns differ substantially from WGS data.**

To overcome this limitation, we developed CyriPanel, an integrated pipeline that combines CNVPanelizer with a modified Cyrius framework. **CyriPanel replaces the WGS-dependent GMM with CNVPanelizer's bootstrap-based approach, generating synthetic reference sets directly from panel data. We further optimized variant calling algorithms for high-depth sequencing noise and implemented fallback mechanisms to resolve ambiguous CNV or SV calls.** CyriPanel enables accurate CYP2D6 genotyping from targeted sequencing without requiring extensive reference cohorts, making comprehensive pharmacogenomic testing more accessible for clinical laboratories using custom panels.

---

## Installation Guide
CyriPanel requires both Python (≥3.6) and R (≥4.0) environments with specific packages.

### Step 1: Python Dependencies
#### Option A: Using pip (Recommended)
```bash
# Install Python dependencies
pip install -r requirements.txt
```

#### Option B: Using conda
```bash
# Create a new conda environment
conda create -n cyripanel python=3.9
conda activate cyripanel

# Install dependencies
conda install numpy scipy pysam statsmodels pandas
# or use pip within conda:
pip install -r requirements.txt
```

### Step 2: R Dependencies
#### Prerequisites
- R version ≥4.0.0 (check with `R --version`)
- Bioconductor package manager

#### Installation via R Console
```R
# Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install CNVPanelizer and its dependencies
BiocManager::install("CNVPanelizer")

# Verify installation
library(CNVPanelizer)
packageVersion("CNVPanelizer")  # Should be ≥1.24.0
```

#### Alternative: Command Line Installation
```bash
# Run R commands from shell
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
Rscript -e 'BiocManager::install("CNVPanelizer")'
```

### Step 3: Verify Installation
#### Test Python Environment
```python
import numpy
import scipy
import pysam
import statsmodels
import pandas

print("All Python packages imported successfully!")
print(f"pandas version: {pandas.__version__}")
```

#### Test R Environment
```bash
# Test CNVPanelizer availability
Rscript -e 'library(CNVPanelizer); cat("CNVPanelizer loaded successfully\n")'
```

### Step 4: Prepare Reference Panel
CyriPanel requires a reference panel of diploid samples (CN=2 at CYP2D6):
1. **Collect reference BAM files** from samples with known diploid CYP2D6 (recommend to have around 20 samples)
2. **Ensure all BAM files have index files** (.bai)
3. **Place these BAM and BAI files in the ref_dir:**
   ```
   CyriPanel/
   ├── star_caller.py
   ├── depth_calling
   ├── caller
   └── ref_dir/
       ├── sample01.bam
       ├── sample01.bam.bai
       ├── sample02.bam
       ├── sample02.bam.bai
       └── ...
   ```

---

## Common Installation Issues

### Issue 1: R Package Installation Fails
**Symptoms:** BiocManager::install() fails with compilation errors

**Solutions:**
```bash
# On Ubuntu/Debian
sudo apt-get install build-essential r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev

# On macOS
xcode-select --install
brew install libxml2
```

### Issue 2: pysam Installation Fails
**Symptoms:** Error compiling pysam C extensions

**Solutions:**
```bash
# On Ubuntu/Debian
sudo apt-get install libbz2-dev liblzma-dev

# On macOS
brew install xz bzip2

# Then retry pip install
pip install --upgrade pysam
```

### Issue 3: CNVPanelizer Not Found in Subprocess
**Symptoms:** Python script can't find Rscript or CNVPanelizer

**Solution:**
```bash
# Verify R is in PATH
which Rscript
R --version

# Set R library path explicitly (add to ~/.bashrc or ~/.zshrc)
export R_LIBS_USER="$HOME/R/library"

# Or specify full Rscript path in panel_cn.py
```

### Issue 4: Personal Library Path Issues
**Symptoms:** Permission denied when installing R packages

**Solution:** The run_CNVPanelizer.R script automatically handles this by setting:
```R
personal_lib_path <- file.path(Sys.getenv("HOME"), "R", "library")
.libPaths(c(personal_lib_path, .libPaths()))
```

---

## Version Compatibility Matrix

| Python | numpy  | scipy | pysam  | pandas | R    | CNVPanelizer |
|--------|--------|-------|--------|--------|------|--------------|
| 3.6    | 1.16+  | 1.2+  | 0.15.3+| 1.0+   | 4.0+ | 1.24.0+      |
| 3.7    | 1.17+  | 1.3+  | 0.16+  | 1.1+   | 4.1+ | 1.24.0+      |
| 3.8    | 1.18+  | 1.4+  | 0.16+  | 1.2+   | 4.1+ | 1.26.0+      |
| 3.9    | 1.19+  | 1.5+  | 0.17+  | 1.3+   | 4.2+ | 1.26.0+      |
| 3.10   | 1.21+  | 1.7+  | 0.18+  | 1.4+   | 4.2+ | 1.28.0+      |
| 3.11   | 1.23+  | 1.9+  | 0.19+  | 1.5+   | 4.3+ | 1.28.0+      |

**Note:** These are tested combinations. Other version combinations may work but are not guaranteed.

---

## Quick Start Verification Script

Save as `test_installation.sh`:
```bash
#!/bin/bash
echo "Testing CyriPanel Installation..."
echo ""

echo "1. Testing Python environment..."
python3 -c "import numpy, scipy, pysam, statsmodels, pandas; print('✓ Python packages OK')"

echo "2. Testing R environment..."
Rscript -e 'library(CNVPanelizer); cat("✓ R packages OK\n")'

echo "3. Testing Rscript in PATH..."
which Rscript && echo "✓ Rscript found in PATH"

echo ""
echo "Installation verification complete!"
```

Run with: `bash test_installation.sh`

---

## Support

For installation issues:
1. Check the [CyriPanel GitHub Issues](your-github-url)
2. Verify your environment matches the compatibility matrix
3. Review common issues section above
4. Contact: yhlin.md05@nycu.edu.tw
