#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_name="${1:-analysis-forward.yml}"

exec Rscript "${project_dir}/run-local.R" "${config_name}"
