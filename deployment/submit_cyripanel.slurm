#!/bin/bash
#SBATCH --job-name=CyriPanel
#SBATCH --output=cyripanel_%j.out
#SBATCH --error=cyripanel_%j.err
#SBATCH --export=ALL

# Before sbatch: activate the desired Conda environment.
# Supply account, partition, CPU, memory and time options to sbatch.
# Supply an absolute repository path followed by run_cyripanel.sh arguments.
# Scheduler logs default to the submission directory.
set -euo pipefail
usage() {
    echo "Usage: sbatch [cluster options] deployment/submit_cyripanel.slurm /absolute/repository [run_cyripanel.sh arguments]"
}
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage; exit 0
fi
[[ $# -ge 2 ]] || { usage >&2; exit 2; }
repo_dir=$1
shift
[[ "$repo_dir" == /* && -f "$repo_dir/run_cyripanel.sh" ]] || {
    echo "ERROR: First argument must be an absolute CyriPanel repository path." >&2; exit 2;
}
[[ -n "${CONDA_PREFIX:-}" ]] || {
    echo "ERROR: Activate the intended Conda environment before sbatch." >&2; exit 2;
}
for tool in python3 Rscript; do
    tool_path=$(command -v "$tool") || {
        echo "ERROR: Missing tool: $tool" >&2; exit 2;
    }
    [[ "$tool_path" == "$CONDA_PREFIX/bin/"* ]] || {
        echo "ERROR: $tool is outside the active Conda environment." >&2; exit 2;
    }
done
bash "$repo_dir/run_cyripanel.sh" "$@"
echo "CyriPanel job completed."
