#!/usr/bin/env bash
# Source from a Bash shell in which Conda has already been initialized.
# Usage: source activate_cyripanel.sh [ENV_NAME]

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf '%s\n' "Use: source activate_cyripanel.sh [ENV_NAME]" >&2
    exit 2
fi
if [[ $# -gt 1 ]]; then
    printf '%s\n' "Use: source activate_cyripanel.sh [ENV_NAME]" >&2
    return 2
fi
conda activate "${1:-cyripanel-example}" || return $?
printf 'Activated environment: %s\n' "${CONDA_DEFAULT_ENV}"
