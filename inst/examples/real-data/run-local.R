#!/usr/bin/env Rscript

# Run every file-based PheWAS stage for one configuration. Relative
# configuration paths are resolved beside this script so the command works
# from any working directory, including when the project path contains spaces.

.phewasflow_runner_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    stop("Cannot locate run-local.R; pass `project_dir` explicitly.",
         call. = FALSE)
  }
  script <- sub("^--file=", "", file_arg[[1L]])
  if (!file.exists(script)) {
    script <- gsub("~+~", " ", script, fixed = TRUE)
  }
  dirname(normalizePath(
    script,
    winslash = "/",
    mustWork = TRUE
  ))
}

.phewasflow_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}

run_phewasflow_local <- function(config_path = "analysis-forward.yml",
                                 project_dir = NULL) {
  if (!requireNamespace("phewasFlow", quietly = TRUE)) {
    stop("Install phewasFlow before running this project.", call. = FALSE)
  }
  if (!is.character(config_path) || length(config_path) != 1L ||
      is.na(config_path) || !nzchar(config_path)) {
    stop("`config_path` must be one non-empty path.", call. = FALSE)
  }
  if (is.null(project_dir)) {
    project_dir <- .phewasflow_runner_directory()
  }
  project_dir <- normalizePath(
    path.expand(project_dir), winslash = "/", mustWork = TRUE
  )
  config_path <- path.expand(config_path)
  if (!.phewasflow_absolute_path(config_path)) {
    config_path <- file.path(project_dir, config_path)
  }
  config_path <- normalizePath(config_path, winslash = "/", mustWork = TRUE)

  run_cli <- function(arguments) {
    status <- phewasFlow::phewasflow_cli(arguments)
    if (!identical(as.integer(status), 0L)) {
      stop("phewasFlow command failed: ", paste(arguments, collapse = " "),
           call. = FALSE)
    }
    invisible(status)
  }

  run_cli(c("validate", "--config", config_path))
  config <- phewasFlow::read_phewasflow_config(
    config_path, check_data = FALSE
  )
  run_cli(c("manifest", "--config", config_path))
  for (shard_id in seq_len(config$execution$shards)) {
    run_cli(c(
      "run-shard", "--config", config_path,
      "--shard-id", as.character(shard_id)
    ))
  }
  run_cli(c("combine", "--config", config_path))

  result <- phewasFlow::read_phewas_bundle(config$output$directory)
  plottable <- result$status == "ok" & is.finite(result$neg_log10_p)
  if (any(plottable, na.rm = TRUE)) {
    run_cli(c("plot", "--config", config_path))
  } else {
    message(
      "No successful association has a finite plotting value; ",
      "the Manhattan and volcano plots were not created."
    )
  }

  status_levels <- c("ok", "skipped", "error")
  status_counts <- table(factor(
    as.character(result$status), levels = status_levels
  ))
  message(sprintf(
    "Association status: ok=%d, skipped=%d, error=%d.",
    status_counts[["ok"]],
    status_counts[["skipped"]],
    status_counts[["error"]]
  ))
  if (any(result$status != "ok", na.rm = TRUE)) {
    message(
      "Inspect status, reason_code, message, and warnings in results.tsv ",
      "before interpreting the analysis."
    )
  }
  message("Results: ", config$output$directory)
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) {
    stop("Usage: Rscript run-local.R [CONFIG.yml]", call. = FALSE)
  }
  config_path <- if (length(args)) args[[1L]] else "analysis-forward.yml"
  run_phewasflow_local(config_path)
}
