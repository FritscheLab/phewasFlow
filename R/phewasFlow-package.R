#' phewasFlow: reproducible bidirectional PheWAS
#'
#' `phewasFlow` runs an explicitly specified association model for every
#' phenotype in an analysis-ready phenome. The phenotype can be the response
#' (for example, a PGS-to-phenome scan) or the scanned predictor (a reverse
#' PheWAS). All supported outcome families share one result contract.
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

