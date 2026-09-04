# Contributing

Keep statistical behavior documented and tested. Preserve caller-owned inputs.
Update deterministic tests and documentation when model behavior or result
fields change.

Use the installation and dependency commands in the
[README](https://github.com/FritscheLab/phewasFlow#development), reinstall the
current source so parallel workers use it, then run:

```r
testthat::test_local()
rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "warning")
```

Do not commit individual-level health data, fitted model objects, site-specific
paths, credentials, or generated analysis results.

Regenerate function documentation with `roxygen2::roxygenise()` after changing
documented interfaces. Preview the package website with
`pkgdown::build_site(preview = FALSE)`; it uses the Fritsche Lab palette and links
back to the lab website. The site is generated into ignored `docs/`.

Before a public release:

1. Confirm the approved license and copyright statement in `LICENSE`.
2. Run the package checks, including vignettes and the simulated workflows.
3. Review `git ls-files` and the contents of the source tarball. Publish from
   Git or the checked tarball; a copy of a working directory can include ignored
   private inputs and local workflows.
4. Confirm repository and support URLs, update the version and citation metadata,
   and review generated outputs before sharing them.

Only the source, documentation, and synthetic examples belong in this repository.
Git ignore rules do not remove previously tracked files or clean Git history.
