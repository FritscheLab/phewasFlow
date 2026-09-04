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

1. Keep the MIT license declaration and copyright statement consistent across
   `DESCRIPTION`, `LICENSE`, and `LICENSE.md`.
2. Run the package checks, including vignettes and the simulated workflows.
3. Review `git ls-files` and the contents of the source tarball. Publish from
   Git or the checked tarball; a copy of a working directory can include ignored
   private inputs and local workflows.
4. Confirm repository and support URLs, update the version and citation metadata,
   and review generated outputs before sharing them.

## Archive releases and maintain citations

Use `CITATION.cff` as the metadata source for GitHub and Zenodo. Keep its title,
authors, version, license, and project URLs consistent with `DESCRIPTION` and
`inst/CITATION`. Add authors and ORCIDs only when verified.

For the first archived release:

1. Connect the maintainer's GitHub account to Zenodo and
   [enable this repository](https://help.zenodo.org/docs/github/enable-repository/).
   Organization repositories may require organization access approval.
2. Set the release version in `DESCRIPTION` and `CITATION.cff`, add the actual
   `date-released` to `CITATION.cff`, and check the citation returned by
   `citation("phewasFlow")` after reinstalling the package.
3. Publish the checked GitHub release and wait for
   [Zenodo to archive it](https://help.zenodo.org/docs/github/archive-software/github-upload/).
4. Record the assigned DOI in `CITATION.cff` (`doi`) and `inst/CITATION`
   (the `doi` argument to `bibentry`), and link it from the README citation
   section. Use the DOI for the matching release and update the citation year
   to match its release date. Do not pair a previous release's DOI with a new
   version; obtain the matching DOI when archiving each release.

Keep scholarly citation requests in the documentation and citation files. Do
not add citation conditions to the standard MIT license text.

Only the source, documentation, and synthetic examples belong in this repository.
Git ignore rules do not remove previously tracked files or clean Git history.
