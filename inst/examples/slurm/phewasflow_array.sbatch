#!/usr/bin/env bash
#SBATCH --job-name=phewasflow
# This explicit 1-10 range matches execution.shards: 10 in the example YAML.
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

: "${PHEWASFLOW_CONFIG:?Set PHEWASFLOW_CONFIG to the YAML configuration path}"

phewasflow run-shard \
  --config "${PHEWASFLOW_CONFIG}" \
  --shard-id "${SLURM_ARRAY_TASK_ID}"
