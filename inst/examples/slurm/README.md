# Run a sharded PheWAS with SLURM

This directory contains two reusable batch templates:

- `phewasflow_array.sbatch` runs one phenotype shard per array task.
- `phewasflow_finish.sbatch` combines every shard, applies multiple-testing
  correction, and creates Manhattan and volcano plots in PNG and PDF.

Use these templates after the analysis has succeeded locally on simulated data
or a small representative phenotype set.

## 1. Install and check the command

From the package source directory:

```bash
Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_local(".", dependencies = NA, upgrade = "never")'
export PATH="$PWD/exec:$PATH"
```

The optional conda environment in the package README is also supported.
Reinstall the package after changing package code. Confirm that the submission
shell sees the intended installation:

```bash
command -v Rscript
command -v phewasflow
Rscript -e 'library(phewasFlow); packageVersion("phewasFlow")'
```

The batch jobs must see the same R environment and `phewasflow` executable.
The examples below export the current environment to SLURM. If the cluster does
not preserve it, add the cluster's normal environment initialization to both
batch templates.

## 2. Prepare a configuration

Start from a generated demo project or copy
`inst/examples/config/phewasflow.yml` into the analysis project. Set the input
files, analysis ID, model, eligibility thresholds, shard count, output
directory, and plot settings.

Keep paths relative to the YAML file when the project layout permits it. This
makes the same configuration usable in an interactive shell and in a batch
job.

Before a large run, make a separate smoke-test configuration with:

- a few representative phenotypes;
- a small shard count;
- its own `analysis_id` and output directory; and
- the same variable coding, covariates, and model settings planned for the
  larger analysis.

Run the smoke test through `combine` and `plot`, then review all result
statuses and sample counts.

## 3. Validate and create the manifest

From the directory containing the YAML file:

```bash
export PHEWASFLOW_CONFIG=analysis.yml

phewasflow validate --config "$PHEWASFLOW_CONFIG"
phewasflow manifest --config "$PHEWASFLOW_CONFIG"
```

Validation checks the YAML, participant data, phenotype metadata, and model
specification. The manifest assigns every requested phenotype to exactly one
shard and is written under:

```text
<output.directory>/.phewasflow/manifest.tsv
```

Read the `shard_count` printed by `manifest`. It must equal
`execution.shards` in the YAML and the array range submitted in the next
step.

## 4. Submit the array and finish job

The supplied array template contains example CPU, memory, and time settings.
Adjust those directives to the needs of the models and the scheduler's local
policy.

For a configuration with ten shards:

```bash
export PHEWASFLOW_CONFIG=analysis.yml
N_SHARDS=10
PHEWASFLOW_SLURM=$(Rscript -e 'cat(system.file("examples", "slurm", package = "phewasFlow"))')
mkdir -p logs

ARRAY_JOB_ID=$(sbatch --parsable \
  --array="1-${N_SHARDS}" \
  --output="logs/%x-%A_%a.out" \
  --error="logs/%x-%A_%a.err" \
  --export=ALL \
  "$PHEWASFLOW_SLURM/phewasflow_array.sbatch")

sbatch \
  --dependency="afterok:${ARRAY_JOB_ID}" \
  --output="logs/%x-%j.out" \
  --error="logs/%x-%j.err" \
  --export=ALL \
  "$PHEWASFLOW_SLURM/phewasflow_finish.sbatch"
```

Set `N_SHARDS` to the value in the YAML. The command-line `--array` option
overrides the example range in the template. `SLURM_ARRAY_TASK_ID` is passed
to `phewasflow run-shard` as the shard ID.

The finish job uses an `afterok` dependency, so it starts only after every
array task exits successfully. It then runs:

```bash
phewasflow combine --config "$PHEWASFLOW_CONFIG"
phewasflow plot --config "$PHEWASFLOW_CONFIG"
```

## 5. Monitor the jobs

Use the job ID printed during submission:

```bash
squeue -j "$ARRAY_JOB_ID"
sacct -j "$ARRAY_JOB_ID" --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

The example submission writes one output and error log per array task under
`logs/`. If an array task fails, read that task's error log before submitting
again.

## 6. Review the output

The configured output directory contains:

```text
results.rds
results.tsv
phewas.png
phewas.pdf
phewas-volcano.png
phewas-volcano.pdf
.phewasflow/
  manifest.tsv
  shards/
    shard-0001.rds
    shard-0002.rds
    ...
```

`plot.output` sets the plot stem. The finish job writes Manhattan and volcano
plots as both PNG and PDF; the names above use the stem `phewas`. Bonferroni is
the default displayed threshold. Start with the status counts in `results.tsv`:

- `ok` rows were included in multiple-testing correction;
- `skipped` rows did not meet an eligibility or estimability requirement;
  and
- `error` rows need review before the analysis is interpreted.

For non-`ok` rows, inspect `reason_code`, `message`, `warnings`, and the
sample counts. The combined result applies BH and Bonferroni correction once
across all successful rows, not separately within each shard.

## Restart an interrupted run

Submit the same array again with the same YAML after a temporary failure.
Every valid completed shard is checked and reused; only missing or invalid
work prevents the finish step from completing.

The combine step requires the exact shard set described by the manifest. It
rejects missing, unexpected, stale, or mismatched artifacts instead of mixing
runs.

When an input file, model setting, threshold, shard count, or package version
changes, use a new `analysis_id` and output directory. Do not place the new
run in a directory containing shards from an earlier specification.
