#' phewasFlow: reproducible bidirectional PheWAS
#'
#' `phewasFlow` runs an explicitly specified association model for every
#' phenotype in an analysis-ready phenome. The phenotype can be the response
#' (for example, a PGS-to-phenome scan) or the scanned predictor (a reverse
#' PheWAS). All supported outcome families share one result contract.
#'
#' The package complements the broader PheWAS toolkit with validated analysis
#' specifications, restartable sharded execution, detailed result records,
#' matched model contrasts, and count and ordinal outcome models. PheWAS can
#' construct the input phenotypes; it is not required to run the model wrappers
#' in this package. Firth regression and local parallel execution are supported
#' by both packages and are not unique contributions of phewasFlow.
#'
#' Upper-tail probabilities are calculated in log space and retained through
#' multiple-testing correction. Text and negative-log representations preserve
#' very small p-values for reporting and plotting when ordinary probabilities
#' lose precision or underflow.
#'
#' See `vignette("phewas-interoperability", package = "phewasFlow")` for the
#' input bridge and an optional numerical and runtime comparison. The package
#' does not claim greater accuracy or speed than PheWAS.
#'
#' The package never constructs clinical phenotypes from raw diagnosis codes
#' and never treats missing phenotypes as controls. Individual-level data and
#' fitted model objects are not written by the sharded workflow.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats coef confint complete.cases lm model.frame model.matrix
#'   p.adjust pchisq pnorm pt quantile reformulate sd setNames terms vcov
#' @importFrom utils read.csv read.delim write.table
NULL
