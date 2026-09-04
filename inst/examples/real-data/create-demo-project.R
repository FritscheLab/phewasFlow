#!/usr/bin/env Rscript

# Create a writable, self-contained project that exercises the same file
# interface used for a real PheWAS. The generated participants are synthetic.

find_template_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    stop("Cannot locate the template directory; pass `template_dir` explicitly.",
         call. = FALSE)
  }
  script <- sub("^--file=", "", file_arg[[1L]])
  if (!file.exists(script)) {
    script <- gsub("~+~", " ", script, fixed = TRUE)
  }
  dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
}

create_phewasflow_demo_project <- function(destination, template_dir = NULL) {
  if (!is.character(destination) || length(destination) != 1L ||
      is.na(destination) || !nzchar(destination)) {
    stop("`destination` must be one non-empty path.", call. = FALSE)
  }
  if (is.null(template_dir)) {
    template_dir <- find_template_directory()
  }
  template_dir <- normalizePath(template_dir, winslash = "/", mustWork = TRUE)
  destination <- path.expand(destination)

  if (dir.exists(destination)) {
    existing <- list.files(destination, all.files = TRUE, no.. = TRUE)
    if (length(existing)) {
      stop("Destination already exists and is not empty: ", destination,
           call. = FALSE)
    }
  } else if (!dir.create(destination, recursive = TRUE)) {
    stop("Could not create destination: ", destination, call. = FALSE)
  }
  destination <- normalizePath(destination, winslash = "/", mustWork = TRUE)

  project_files <- c(
    "README.md",
    "analysis-forward.yml",
    "analysis-reverse.yml",
    "phenotype-metadata-outcomes.tsv",
    "phenotype-metadata-predictors.tsv",
    "analysis-data-dictionary.tsv",
    "run-local.R",
    "run-local.sh"
  )
  source_paths <- file.path(template_dir, project_files)
  missing <- project_files[!file.exists(source_paths)]
  if (length(missing)) {
    stop("Template is incomplete; missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  copied <- file.copy(source_paths, destination, overwrite = FALSE)
  if (!all(copied)) {
    stop("Could not copy every project template.", call. = FALSE)
  }
  if (!file.copy(
    file.path(template_dir, "project.gitignore"),
    file.path(destination, ".gitignore"),
    overwrite = FALSE
  )) {
    stop("Could not create the project .gitignore.", call. = FALSE)
  }

  input_dir <- file.path(destination, "inputs")
  if (!dir.create(input_dir)) {
    stop("Could not create the input directory.", call. = FALSE)
  }

  example <- phewasFlow::phewas_example_data(n = 220L, seed = 401L)
  saveRDS(
    as.data.frame(example$data, stringsAsFactors = FALSE),
    file.path(input_dir, "analysis.rds"),
    version = 3L
  )
  writeLines(
    c(
      "SYNTHETIC DEMONSTRATION DATA ONLY",
      "Created by phewasFlow::phewas_example_data(n = 220, seed = 401).",
      "Replace analysis.rds before conducting a scientific analysis."
    ),
    file.path(input_dir, "DEMO_DATA.txt")
  )

  forward <- file.path(destination, "analysis-forward.yml")
  reverse <- file.path(destination, "analysis-reverse.yml")
  phewasFlow::read_phewasflow_config(forward, check_data = TRUE)
  phewasFlow::read_phewasflow_config(reverse, check_data = TRUE)

  message("Created a validated synthetic PheWAS project at: ", destination)
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  message(
    "Run: ", shQuote(rscript), " ",
    shQuote(file.path(destination, "run-local.R")),
    " analysis-forward.yml"
  )
  invisible(destination)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop(
      "Usage: Rscript create-demo-project.R DESTINATION_DIRECTORY",
      call. = FALSE
    )
  }
  create_phewasflow_demo_project(args[[1L]])
}
