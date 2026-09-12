# Use PheWAS and phewasFlow together

PheWAS can construct clinical phenotypes and provide their descriptions and
groups. phewasFlow can then run a declared analysis, resume shards, and compare
adjustment models on matched participants. Both packages already support Firth
regression and local parallel execution.

These examples use synthetic data only. They neither validate clinical-code
mappings nor establish clinical findings.

## Install the optional comparison dependency

Install phewasFlow using the repository README. PheWAS is optional for normal
phewasFlow use. For these examples, install the upstream revision pinned in
`DESCRIPTION`:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github(
  "PheWAS/PheWAS@55dd1c24e228851922400cfba8d7db474565ccc7",
  dependencies = NA, upgrade = "never"
)
```

The scripts call the installed packages. Reinstall phewasFlow after editing its
source. Run the commands below from the repository root. Installed copies of
the scripts are also available under
`system.file("examples", "phewas", package = "phewasFlow")`.

## Run the phenotype bridge

```bash
Rscript inst/examples/phewas/use-phewas-phenotypes.R
```

The example creates 400 simulated participants and dated ICD9CM records,
constructs phenotypes with `PheWAS::createPhenotypes()`, gets metadata with
`PheWAS::addPhecodeInfo()`, joins participant tables with
`phewasFlow::assemble_phewas_data()`, and runs a PRS scan. It prints results
without saving participant data or fitted models. Source the file and call
`run_phewas_bridge()` to inspect the returned inputs and results in R.

The bridge attaches PheWAS with `library(PheWAS)`: the pinned version's
`addPhecodeInfo()` helper looks up its `pheinfo` dataset on the search path.

The bridge preserves exact column names, including numeric-looking phecodes.
Logical `FALSE` is declared explicitly as the control reference. `NA` stays
missing, including phenotype exclusions and records below the case threshold.
The metadata join must cover every selected phenotype exactly once; unknown
codes require explicit metadata rather than silently disappearing. For real
data, review code mappings, observation windows, population eligibility, and
exclusion definitions before this stage.

Exclusion mappings also introduce phenotype columns with no observed cases in
this simulation. The bridge keeps their `skipped` results with the reason
`insufficient_cases`; they do not enter the successful testing family.

## Compare numerical results and runtime

```bash
Rscript inst/examples/phewas/compare-packages.R output/phewas-comparison
```

Optional positional arguments set participant count, phenotype count **per
outcome type**, and measured repetitions:

```bash
Rscript inst/examples/phewas/compare-packages.R output/phewas-comparison-large 5000 10 5
```

With the pinned versions, the larger example includes very small p-values that
differ beyond the comparison tolerance, while coefficients and standard errors
agree. These differences arise from loss of precision when calculating a small
tail probability by subtracting from one. The script saves its results before
exiting nonzero. See the
[recorded comparison](https://fritschelab.github.io/phewasFlow/articles/phewas-interoperability.html#recorded-comparison)
for the tables and a brief explanation.

Use a new or empty output directory. The default is 1,000 participants, four
binary outcomes, four continuous outcomes, and three measured repetitions.
The seed is 20260912. Missing anchor, covariate, and phenotype values exercise
the complete-case selection in both packages.

| Outcome | PheWAS call | phewasFlow call |
|---|---|---|
| Binary | `phewas_ext(method = "logistf")` | `run_phewas()` with binary metadata |
| Continuous | `phewas()` | `run_phewas()` with continuous metadata |

Both packages receive the same participant table, predictor, covariates,
reference orientation, and eligibility thresholds. PRS standardization is
disabled in both. Each uses one worker. The comparison disables genetic
allele-frequency/HWE calculations in PheWAS because the predictor is a PRS.
It compares unadjusted association p-values; multiple-testing correction is
outside the timed workload.

The binary wrappers both use `logistf` with profile likelihood. phewasFlow also
explicitly computes the predictor likelihood-ratio test and allows up to 1,000
iterations; PheWAS uses the underlying defaults. These are comparisons of the
actual wrappers, not identical internal call sequences. Confidence intervals
are not compared. Continuous models use the same Gaussian linear model, through
PheWAS's GLM wrapper and phewasFlow's `lm` wrapper.

Each package and outcome type receives one unmeasured warmup. Measured repeats
alternate package order and run garbage collection before starting the clock.
Timings include scan validation and result assembly. They exclude simulation,
specification construction, package loading, plotting, disk writes, and
multiple-testing correction.

The output contains:

- `agreement.tsv`: coefficients, standard errors, natural-log p-values, sample
  counts, upstream notes, and acceptance flags for the warmup and each repeat;
- `timings.tsv` and `runtime-summary.tsv`: elapsed seconds and median/range;
- `metadata.json`: parameters, seed, package versions, upstream revision, and
  acceptance tolerances;
- `session-info.txt`: the R session and dependency versions; and
- `README.md`: a summary of agreement, runtime, and limits.

Acceptance requires successful phewasFlow fits, identical sample counts
(including cases and controls for binary outcomes), absolute coefficient and
standard-error differences at most `1e-6`, and absolute natural-log p-value
differences at most `1e-4`. A failed check produces a nonzero exit after writing
the comparison tables. Zero or nonfinite reported p-values fail the comparison;
`p_comparison` identifies zero upstream values separately from other tolerance
failures. The default fixture checks moderate associations; larger fixtures can
include more extreme p-values. There is no timing-based pass criterion.

This is a reproducible small-workload comparison, not a scalability benchmark.
It does not test count or ordinal models, reverse scans, rare or separated
outcomes, cluster recovery, or phenotype validity. Package tests separately
exercise model engines, failed/skipped associations, matched contrasts, and
sharding. Agreement does not establish greater accuracy, and timing ratios
depend on the data, hardware, dependency versions, and wrapper behavior.

The optional integration tests in `tests/testthat/test-phewas-interoperability.R`
check the bridge and shared-model agreement without asserting timing thresholds.
