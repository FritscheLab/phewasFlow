# Compare shared, synthetic forward scans. See README.md for methods and limits.
# No individual-level data or fitted models are written to disk.

make_phewas_comparison <- function(n = 1000L, per_type = 4L, seed = 20260912L) {
  if (length(n) != 1L || !is.finite(n) || n < 200L || n != as.integer(n) ||
      length(per_type) != 1L || !is.finite(per_type) || per_type < 1L ||
      per_type != as.integer(per_type)) {
    stop("Use integer n >= 200 and per_type >= 1.", call. = FALSE)
  }
  set.seed(seed)
  data <- data.frame(
    prs = stats::rnorm(n), age = stats::rnorm(n),
    sex = factor(sample(c("Female", "Male"), n, replace = TRUE),
                 levels = c("Female", "Male"))
  )
  metadata <- data.frame(
    phenotype = c(sprintf("binary_%03d", seq_len(per_type)),
                  sprintf("continuous_%03d", seq_len(per_type))),
    description = paste("Simulated outcome", seq_len(2L * per_type)),
    group = "Synthetic", groupnum = 1L, color = "#006f78",
    variable_type = rep(c("binary", "numeric"), each = per_type),
    outcome_type = rep(c("binary", "continuous"), each = per_type),
    reference = rep(c("FALSE", NA_character_), each = per_type)
  )
  for (j in seq_len(nrow(metadata))) {
    linear <- (0.15 + 0.05 * ((j - 1L) %% 4L)) * data$prs +
      0.15 * data$age + 0.1 * (data$sex == "Male")
    value <- if (metadata$outcome_type[[j]] == "binary") {
      stats::runif(n) < stats::plogis(-0.5 + linear)
    } else {
      linear + stats::rnorm(n)
    }
    value[sample.int(n, floor(0.03 * n))] <- NA
    data[[metadata$phenotype[[j]]]] <- value
  }
  data$prs[sample.int(n, floor(0.01 * n))] <- NA
  data$age[sample.int(n, floor(0.01 * n))] <- NA
  specs <- lapply(c("binary", "continuous"), function(type) {
    phewasFlow::phewas_spec(
      paste0("comparison_", type), "prs", "phenotypes_as_outcomes",
      metadata[metadata$outcome_type == type, ],
      covariates = c("age", "sex"), anchor_type = "numeric",
      anchor_transform = "none",
      eligibility = list(
        min_n = 20L, min_cases = 20L, min_controls = 20L,
        min_outcome_levels = 3L
      ),
      fdr_threshold = 0.05
    )
  })
  names(specs) <- c("binary", "continuous")
  list(data = data, metadata = metadata, specs = specs)
}

fit_phewas_comparison <- function(fixture, type, package) {
  if (package == "phewasFlow") {
    return(phewasFlow::run_phewas(
      fixture$data, fixture$specs[[type]], backend = "sequential", adjust = FALSE
    ))
  }
  arguments <- list(
    phenotypes = fixture$specs[[type]]$phenotypes$phenotype,
    genotypes = "prs", data = fixture$data, covariates = c("age", "sex"),
    cores = 1L, additive.genotypes = FALSE, min.records = 20L
  )
  if (type == "binary") {
    do.call(PheWAS::phewas_ext, c(arguments, list(method = "logistf")))
  } else {
    do.call(PheWAS::phewas, arguments)
  }
}

check_phewas_agreement <- function(flow, upstream, type) {
  flow <- as.data.frame(flow)
  if (anyDuplicated(flow$phenotype) || anyDuplicated(upstream$phenotype) ||
      !setequal(flow$phenotype, upstream$phenotype)) {
    stop("The packages returned different or duplicate phenotype sets.", call. = FALSE)
  }
  upstream <- upstream[match(flow$phenotype, upstream$phenotype), ]
  out <- data.frame(
    outcome_type = type, phenotype = flow$phenotype,
    flow_status = flow$status, upstream_note = upstream$note,
    n_flow = flow$n_complete, n_phewas = upstream$n_total,
    beta_flow = flow$estimate, beta_phewas = upstream$beta,
    se_flow = flow$std_error, se_phewas = upstream$SE,
    p_flow = flow$p_value, p_phewas = upstream$p,
    log_p_flow = flow$log_p, log_p_phewas = log(upstream$p)
  )
  out$beta_abs_difference <- abs(out$beta_flow - out$beta_phewas)
  out$se_abs_difference <- abs(out$se_flow - out$se_phewas)
  out$p_abs_difference <- abs(out$p_flow - out$p_phewas)
  out$log_p_abs_difference <- abs(out$log_p_flow - out$log_p_phewas)
  out$p_comparison <- ifelse(
    out$p_phewas == 0 & out$p_flow > 0, "upstream_reported_zero",
    ifelse(is.finite(out$log_p_abs_difference) & out$log_p_abs_difference <= 1e-4,
           "within_tolerance", "outside_tolerance")
  )
  out$counts_match <- out$n_flow == out$n_phewas
  if (type == "binary") {
    out$counts_match <- out$counts_match &
      flow$n_event == upstream$n_cases & flow$n_nonevent == upstream$n_controls
  }
  out$pass <- out$flow_status == "ok" & out$counts_match &
    is.finite(out$beta_abs_difference) & out$beta_abs_difference <= 1e-6 &
    is.finite(out$se_abs_difference) & out$se_abs_difference <= 1e-6 &
    is.finite(out$log_p_abs_difference) & out$log_p_abs_difference <= 1e-4
  out$pass[is.na(out$pass)] <- FALSE
  out
}

run_phewas_comparison <- function(directory, n = 1000L, per_type = 4L,
                                  repetitions = 3L, seed = 20260912L) {
  for (package in c("PheWAS", "phewasFlow")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install ", package, " before running the comparison.", call. = FALSE)
    }
  }
  if (length(repetitions) != 1L || !is.finite(repetitions) ||
      repetitions < 1L || repetitions != as.integer(repetitions)) {
    stop("repetitions must be a positive integer.", call. = FALSE)
  }
  if (file.exists(directory) && !dir.exists(directory)) {
    stop("Output path is an existing file.", call. = FALSE)
  }
  if (dir.exists(directory) && length(list.files(directory, all.files = TRUE,
                                                no.. = TRUE))) {
    stop("Choose a new or empty output directory.", call. = FALSE)
  }
  fixture <- make_phewas_comparison(n, per_type, seed)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  packages <- c("PheWAS", "phewasFlow")
  agreements <- timings <- list()
  for (type in names(fixture$specs)) {
    # Load paths and warm up each method once outside the measured repeats.
    fitted <- setNames(lapply(packages, function(package) {
      suppressMessages(fit_phewas_comparison(fixture, type, package))
    }), packages)
    warmup <- check_phewas_agreement(fitted$phewasFlow, fitted$PheWAS, type)
    warmup$repetition <- 0L
    agreements[[length(agreements) + 1L]] <- warmup
    for (i in seq_len(repetitions)) {
      order <- if (i %% 2L) packages else rev(packages)
      for (package in order) {
        invisible(gc())
        elapsed <- system.time({
          fitted[[package]] <- suppressMessages(
            fit_phewas_comparison(fixture, type, package)
          )
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          outcome_type = type, package = package, repetition = i,
          elapsed_seconds = elapsed, n = n, phenotypes = per_type, workers = 1L
        )
      }
      agreement <- check_phewas_agreement(fitted$phewasFlow, fitted$PheWAS, type)
      agreement$repetition <- i
      agreements[[length(agreements) + 1L]] <- agreement
    }
  }
  agreement <- do.call(rbind, agreements)
  timing <- do.call(rbind, timings)
  summary <- do.call(rbind, lapply(split(timing, list(timing$outcome_type,
                                                    timing$package)), function(x) {
    data.frame(outcome_type = x$outcome_type[[1L]], package = x$package[[1L]],
               median_seconds = median(x$elapsed_seconds),
               min_seconds = min(x$elapsed_seconds), max_seconds = max(x$elapsed_seconds))
  }))
  write_tsv <- function(x, name) {
    utils::write.table(x, file.path(directory, name), sep = "\t",
                       row.names = FALSE, quote = TRUE, na = "NA")
  }
  write_tsv(agreement, "agreement.tsv")
  write_tsv(timing, "timings.tsv")
  write_tsv(summary, "runtime-summary.tsv")
  upstream <- utils::packageDescription("PheWAS")
  versions <- vapply(c(packages, "logistf", "MASS", "dplyr", "tidyr"),
                     function(x) as.character(utils::packageVersion(x)), character(1))
  metadata <- list(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n = n, phenotypes_per_type = per_type, repetitions = repetitions, seed = seed,
    versions = as.list(versions), phewas_revision = upstream$RemoteSha,
    R = R.version.string, platform = R.version$platform,
    tolerances = list(beta_absolute = 1e-6, se_absolute = 1e-6,
                      natural_log_p_absolute = 1e-4),
    all_passed = all(agreement$pass)
  )
  jsonlite::write_json(metadata, file.path(directory, "metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines(capture.output(utils::sessionInfo()), file.path(directory, "session-info.txt"))
  report <- c(
    "# Synthetic PheWAS comparison", "",
    sprintf("%d participants; %d phenotypes per outcome type; %d measured repeats.",
            n, per_type, repetitions),
    "One worker; one unmeasured warmup per package and outcome type; alternating execution order.",
    "", paste("All numerical checks passed:", all(agreement$pass)), "",
    sprintf("Maximum absolute coefficient difference: %.6g.", max(agreement$beta_abs_difference)),
    sprintf("Maximum absolute standard-error difference: %.6g.", max(agreement$se_abs_difference)),
    sprintf("Maximum absolute p-value difference: %.6g.", max(agreement$p_abs_difference)),
    sprintf("Maximum absolute natural-log p-value difference: %.6g.", max(agreement$log_p_abs_difference)),
    sprintf("Phenotypes with zero upstream p-values and positive phewasFlow p-values: %d.",
            sum(agreement$repetition == 0L & agreement$p_comparison == "upstream_reported_zero", na.rm = TRUE)),
    "", "| Outcome | Package | Median seconds | Min | Max |",
    "|---|---|---:|---:|---:|",
    sprintf("| %s | %s | %.3f | %.3f | %.3f |", summary$outcome_type,
            summary$package, summary$median_seconds, summary$min_seconds, summary$max_seconds),
    "", "Versions and upstream revision are recorded in metadata.json; session-info.txt records the R session.",
    "These timings include each scan's validation and result assembly, but exclude simulation, specification construction, plotting, disk I/O, and multiple-testing correction.",
    "Both binary methods use logistf with profile likelihood; phewasFlow also explicitly computes a predictor likelihood-ratio test and uses a larger iteration limit.",
    "Agreement on these forward scans does not establish greater accuracy, phenotype validity, or performance on large, rare, separated, count, ordinal, or reverse-direction analyses."
  )
  writeLines(report, file.path(directory, "README.md"))
  if (!all(agreement$pass)) {
    stop("Numerical agreement failed; inspect agreement.tsv and metadata.json.", call. = FALSE)
  }
  invisible(list(agreement = agreement, timing = timing, summary = summary))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L || length(args) > 4L) {
    stop("Usage: Rscript compare-packages.R OUTPUT_DIRECTORY [N=1000] [PER_TYPE=4] [REPETITIONS=3]",
         call. = FALSE)
  }
  values <- c(1000, 4, 3)
  if (length(args) > 1L) values[seq_len(length(args) - 1L)] <- as.numeric(args[-1L])
  run_phewas_comparison(args[[1L]], values[[1L]], values[[2L]], values[[3L]])
  cat("Comparison written to", args[[1L]], "\n")
}
