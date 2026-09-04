#!/usr/bin/env bash
#SBATCH --job-name=phewasflow-finish
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

: "${PHEWASFLOW_CONFIG:?Set PHEWASFLOW_CONFIG to the YAML configuration path}"

phewasflow combine --config "${PHEWASFLOW_CONFIG}"
phewasflow plot --config "${PHEWASFLOW_CONFIG}"
