#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J cyri_versions
#SBATCH -p ngs7G
#SBATCH --time=00:10:00
#SBATCH -c 1
#SBATCH --mem=7g
#SBATCH -o versions_out.log
#SBATCH -e versions_err.log
#SBATCH --mail-user=yhlin.md05@nycu.edu.tw
#SBATCH --mail-type=BEGIN,END,FAIL
#
# capture_versions.sh
# ===================
# Records the exact toolchain behind the CyriPanel results, for the Methods
# section and for JTM's "Availability of data and materials" software block
# (Other requirements: is a required field).
#
# This takes seconds. Run it on the login node — do not queue it:
#
#     module load biology
#     module load Samtools/1.15.1
#     source /work/u7715055/staging/biology/u7715055/analysis_env/bin/activate
#     bash capture_versions.sh
#
# The most important thing it captures is the CyriPanel commit hash and whether
# the working tree is clean, because the Zenodo archive must point at the code
# that produced the published numbers — not at whatever is in the folder today.
#
# Optional environment overrides:
#   CYRIPANEL_DIR   path to the CyriPanel git checkout  (default below)
#   OUT             output file                         (default CyriPanel_versions.txt)

set -u
CYRIPANEL_DIR=${CYRIPANEL_DIR:-/work/u7715055/staging/biology/u7715055/CyriPanel}
OUT=${OUT:-CyriPanel_versions.txt}
: > "$OUT"

say() { printf '%s\n' "$*" >> "$OUT"; }
hdr() { say ""; say "=== $* ==="; }

say "captured $(date -Iseconds) on $(hostname)"
say "CyriPanel directory: $CYRIPANEL_DIR"

# ─── 1. CyriPanel itself — the part that matters most ────────────────────────
hdr "1. CyriPanel source"
if [ -d "$CYRIPANEL_DIR/.git" ]; then
  say "HEAD commit  : $(git -C "$CYRIPANEL_DIR" rev-parse HEAD)"
  say "HEAD short   : $(git -C "$CYRIPANEL_DIR" rev-parse --short HEAD)"
  say "HEAD date    : $(git -C "$CYRIPANEL_DIR" log -1 --format=%cI)"
  say "branch       : $(git -C "$CYRIPANEL_DIR" rev-parse --abbrev-ref HEAD)"
  say "remote       : $(git -C "$CYRIPANEL_DIR" remote get-url origin 2>/dev/null || echo none)"
  say "tags at HEAD : $(git -C "$CYRIPANEL_DIR" tag --points-at HEAD | tr '\n' ' ')"
  dirty=$(git -C "$CYRIPANEL_DIR" status --porcelain | wc -l)
  if [ "$dirty" -eq 0 ]; then
    say "working tree : CLEAN"
  else
    say "working tree : *** $dirty UNCOMMITTED CHANGES — the hash above does not"
    say "               describe what is on disk. Commit or stash before quoting it. ***"
    git -C "$CYRIPANEL_DIR" status --porcelain | sed 's/^/                 /' >> "$OUT"
  fi
  say ""
  say "recent commits — identify the one that produced the benchmark results:"
  git -C "$CYRIPANEL_DIR" log -15 --date=short --format='  %h  %cd  %s' >> "$OUT" 2>&1
else
  say "  NOT a git checkout at $CYRIPANEL_DIR"
  say "  Record the release tag or Zenodo DOI manually — without it the archived"
  say "  version cannot be tied to the published numbers."
fi

# ─── 2. The dependency pin that decides reproducibility ──────────────────────
hdr "2. Is CNVPanelizer version-pinned in the code?"
RSCRIPT_R="$CYRIPANEL_DIR/depth_calling/run_CNVPanelizer.R"
if [ -f "$RSCRIPT_R" ]; then
  grep -n 'BiocManager::install\|install.packages\|version *=' "$RSCRIPT_R" >> "$OUT" 2>&1
  if grep -q 'BiocManager::install' "$RSCRIPT_R" && ! grep -q 'version *=' "$RSCRIPT_R"; then
    say ""
    say "  *** NOT PINNED. Anyone installing from this script gets the current"
    say "      Bioconductor release, not the version validated here. Fix with:"
    say "        BiocManager::install(\"CNVPanelizer\", version = \"3.16\", update = FALSE)"
  fi
else
  say "  $RSCRIPT_R not found"
fi

# ─── 3. R / Bioconductor / CNVPanelizer ──────────────────────────────────────
hdr "3. R, Bioconductor and CNVPanelizer"
if command -v Rscript >/dev/null 2>&1; then
  Rscript -e '
    .libPaths(c(file.path(Sys.getenv("HOME"), "R", "library"), .libPaths()))
    cat("R version         :", R.version.string, "\n")
    cat("Platform          :", R.version$platform, "\n")
    if (requireNamespace("BiocManager", quietly = TRUE))
      cat("Bioconductor      :", as.character(BiocManager::version()), "\n")
    if (requireNamespace("CNVPanelizer", quietly = TRUE)) {
      d <- packageDescription("CNVPanelizer")
      cat("CNVPanelizer      :", d$Version, "\n")
      cat("  Date/Publication:", d$"Date/Publication", "  <-- year for the citation\n")
      cat("  Built           :", d$Built, "\n")
      cat("  Installed at    :", find.package("CNVPanelizer"), "\n")
    } else cat("CNVPanelizer      : NOT INSTALLED in these .libPaths()\n")
    cat("\nlibrary paths:\n"); for (p in .libPaths()) cat("  ", p, "\n")
    cat("\nCNVPanelizer dependencies:\n")
    for (p in c("Rsamtools","GenomicRanges","IRanges","S4Vectors","BiocGenerics",
                "foreach","doParallel","exomeCopy","GenomeInfoDb","BiocParallel","edgeR"))
      if (requireNamespace(p, quietly = TRUE))
        cat(sprintf("  %-16s %s\n", p, as.character(packageVersion(p))))
  ' >> "$OUT" 2>&1
  say ""
  say "--- DESCRIPTION, verbatim cross-check ---"
  DESC="$HOME/R/library/CNVPanelizer/DESCRIPTION"
  if [ -f "$DESC" ]; then
    grep -E '^(Package|Version|Date/Publication|git_url|git_last_commit|Built):' "$DESC" \
      | sed 's/^/  /' >> "$OUT"
  else
    say "  not at $DESC — see the library paths above"
  fi
else
  say "  Rscript not on PATH. Run:  module load biology  (and an R module) first."
fi

# ─── 4. Python ───────────────────────────────────────────────────────────────
hdr "4. Python toolchain"
say "python3 resolves to: $(command -v python3 2>/dev/null || echo 'not found')"
say "VIRTUAL_ENV        : ${VIRTUAL_ENV:-<none active>}"
say ""
python3 - >> "$OUT" 2>&1 << 'PYEOF'
import sys
print("Python            :", sys.version.split()[0])
try:
    import importlib.metadata as md
    ver = md.version
except Exception:
    import pkg_resources
    ver = lambda p: pkg_resources.get_distribution(p).version
for p in ("pysam", "numpy", "scipy", "pandas", "statsmodels"):
    try:
        print("  %-16s %s" % (p, ver(p)))
    except Exception:
        print("  %-16s NOT INSTALLED" % p)
PYEOF

# ─── 5. Alignment and BAM tools ──────────────────────────────────────────────
hdr "5. Alignment and BAM tools"
for t in samtools bwa bwa-mem2 bcftools picard java; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '  %-10s %s\n' "$t" "$("$t" --version 2>&1 | head -1)" >> "$OUT"
  else
    printf '  %-10s not on PATH\n' "$t" >> "$OUT"
  fi
done

hdr "6. Loaded modules"
if command -v module >/dev/null 2>&1; then
  module list 2>&1 | sed 's/^/  /' >> "$OUT"
else
  say "  no module command in this shell"
fi

say ""
say "=== done ==="
say ""
say "Still to record by hand (nothing on the cluster knows these):"
say "  - Zenodo DOI, minted from the commit above after the pending fixes"
say "  - DRAGEN version and hash-table build"
say "  - deCYPher version/commit; PharmVar release behind data/star_table.txt"
cat "$OUT"
