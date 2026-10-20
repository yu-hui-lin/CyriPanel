#!/usr/bin/env bash
# Create the Linux example environment from the recorded dependency file.

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    printf '%s\n' "Run this script with bash; do not source it." >&2
    return 2
fi
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash setup_environment.sh [NEW_ENV_NAME]

Create a Linux x86_64 environment from environment.yml.
Default name: cyripanel-example
Use a different name when that environment already exists.

Example:
  bash setup_environment.sh cyripanel-rebuild
  conda activate cyripanel-rebuild

For the recorded package builds, see the explicit-file instructions
in README.md. This script does not activate the parent shell.
EOF
}

if [[ $# -eq 1 && ( "$1" == "--help" || "$1" == "-h" ) ]]; then
    usage
    exit 0
fi
if [[ $# -gt 1 ]]; then
    usage >&2
    exit 2
fi
env_name="${1:-cyripanel-example}"
if [[ ! "$env_name" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ||
      "$env_name" == "base" || "$env_name" == "root" ]]; then
    printf '%s\n' "ERROR: Choose a new environment name using letters, digits, _ or -." >&2
    exit 2
fi
if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    printf '%s\n' "ERROR: This environment specification targets Linux x86_64." >&2
    exit 1
fi
if ! command -v conda >/dev/null 2>&1; then
    printf '%s\n' "ERROR: Make Conda available before running this script." >&2
    exit 1
fi
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
spec="$script_dir/environment.yml"
if [[ ! -r "$spec" ]]; then
    printf '%s\n' "ERROR: environment.yml is missing or unreadable." >&2
    exit 1
fi

# Do not pass --yes: Conda 24.7.1 env create can remove an existing
# environment when that option is used.
CONDA_ALWAYS_YES=false CONDA_CHANNEL_PRIORITY=strict CONDA_SUBDIR=linux-64 \
    conda env create --solver libmamba --no-default-packages \
    --name "$env_name" --file "$spec"

printf 'Environment created. Activate it with: conda activate %s\n' "$env_name"
