# Run a PheWAS project

This directory is a working template for a file-based PheWAS. It includes
forward and reverse configurations, phenotype metadata, and one runner that
validates, fits, combines, and plots the results.

## Try the example

Install `phewasFlow`, then run these commands from the package source directory:

```bash
Rscript inst/examples/real-data/create-demo-project.R my-phewas
cd my-phewas
Rscript run-local.R analysis-forward.yml
```

The generated data are simulated. A successful run confirms that the package,
input files, model engines, and plotting packages work together.

If you are already inside the generated project, the run command is simply:

```bash
Rscript run-local.R analysis-forward.yml
```

The project contains:

```text
my-phewas/
  analysis-forward.yml
  analysis-reverse.yml
  phenotype-metadata-outcomes.tsv
  phenotype-metadata-predictors.tsv
  analysis-data-dictionary.tsv
  run-local.R
  inputs/
    analysis.rds
```

## Use your data

1. Replace `inputs/analysis.rds` with a data frame containing one row per
   participant. It must include the anchor, covariates, every phenotype named
   in the metadata, and any declared offset or exclusion columns.
2. Edit `phenotype-metadata-outcomes.tsv` for a PRS, exposure, or biomarker
   tested against many outcomes. Edit `phenotype-metadata-predictors.tsv` when
   many phenotypes are tested against one clinical endpoint.
3. Edit the corresponding YAML file. Set a new `analysis_id`, variable names,
   eligibility thresholds, output directory, and plot path.
4. Run the analysis from the project directory.

```bash
Rscript run-local.R analysis-forward.yml
```

For the reverse direction:

```bash
Rscript run-local.R analysis-reverse.yml
```

Relative input and output paths are resolved from the YAML file's directory,
so the runner also works when called by its full path from another directory.

## Participant table

The participant table may be RDS, CSV, TSV, TAB, or TXT. RDS is recommended
when a covariate is categorical because it preserves factor levels. Character
covariates are not model-ready; for tabular input, use K-1 numeric indicator
columns.

Before saving the table:

- require one unique, nonmissing participant ID;
- encode binary phenotypes consistently, such as `0` for controls and `1` for
  cases;
- keep an unavailable phenotype as `NA`, not as a control;
- convert missing-value codes such as `-9` or `999` to `NA`;
- make counts nonnegative integers and numeric variables finite;
- set factor reference levels explicitly.

If tables need to be joined by ID, use `assemble_phewas_data()` and review its
`join_report` before writing `inputs/analysis.rds`.

## Phenotype metadata

Every metadata row requires:

| Column | Meaning |
|---|---|
| `phenotype` | Exact participant-table column name and result identifier |
| `description` | Readable label |
| `group` | Clinical or phenotype group |
| `groupnum` | Numeric group order |
| `color` | Valid R color used for the group |

`variable_type` may be `numeric`, `binary`, `count`, or `ordinal`. If the
column is absent or a row value is `NA` or blank, the scanned phenotype
defaults to `binary`. For a scanned binary phenotype, a missing, `NA`, or
blank `reference` defaults to character `0`; declare both fields when those
defaults are not correct.

Forward metadata still requires `outcome_type` for every row: `continuous`,
`binary`, `count`, or `ordinal`. An ordinal phenotype requires pipe-separated
`levels`, for example `mild|moderate|severe`. A count outcome can name an
already-log-transformed `offset` column.

Rows in the same group must use the same `groupnum` and `color`. Different
groups must use different `groupnum` values.

## Configuration

Use `direction: phenotypes_as_outcomes` when the fixed `anchor` is the
predictor, such as a PRS. Use `direction: phenotypes_as_predictors` when the
fixed `anchor` is the outcome, such as a continuous symptom score or binary
disease endpoint.

Keep each multiple-testing family separate: one PRS or fixed endpoint per
configuration, with its own `analysis_id` and output directory. Set the
eligibility thresholds and reference values from the analysis plan.
The scanned-phenotype reference default does not apply to the fixed anchor;
always provide `anchor_reference` when that anchor is binary.

## Check the output

When at least one model has a finite p-value, the configured output directory
contains:

```text
results.rds
results.tsv
phewas.png
phewas.pdf
phewas-volcano.png
phewas-volcano.pdf
```

`plot.output` sets the filename stem. A value ending in `.png`, ending in
`.pdf`, or having no extension still creates all four files. PNG is convenient
for review; PDF keeps vector graphics for publication and post-editing. The
Bonferroni significance line is displayed by default.

The Manhattan plot shows associations across phenotype groups. The volcano
plot shows `-log10(p)` against the fitted model or link-scale coefficient. It
uses separate panels when effect measures differ; coefficient magnitudes from
different panels or variables with different units are not directly
comparable.

Start by checking `status` in `results.tsv`:

- `ok`: the association was fit and included in multiple-testing correction;
- `skipped`: it did not meet an eligibility or estimability requirement;
- `error`: validation, fitting, or convergence failed.

Review `reason_code`, `message`, `warnings`, sample counts, case/control counts,
reference values, and `testing_family_size` before interpreting effect sizes.
`native_effect` is a coefficient for continuous outcomes, an odds ratio for
binary outcomes, an incidence-rate ratio for counts, and a common odds ratio
for ordinal outcomes.

For full input examples and troubleshooting, see
`vignette("real-data-workflow", package = "phewasFlow")`.
