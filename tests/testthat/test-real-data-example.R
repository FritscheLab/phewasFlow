real_data_example_directory <- function() {
  path <- system.file("examples", "real-data", package = "phewasFlow")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "examples", "real-data")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

test_that("real-data templates are complete, portable, and well formed", {
  directory <- real_data_example_directory()
  expected <- c(
    "README.md", "analysis-forward.yml", "analysis-reverse.yml",
    "phenotype-metadata-outcomes.tsv",
    "phenotype-metadata-predictors.tsv",
    "analysis-data-dictionary.tsv", "project.gitignore",
    "create-demo-project.R", "run-local.R", "run-local.sh"
  )
  expect_true(all(file.exists(file.path(directory, expected))))

  tsv_paths <- file.path(directory, c(
    "phenotype-metadata-outcomes.tsv",
    "phenotype-metadata-predictors.tsv",
    "analysis-data-dictionary.tsv"
  ))
  for (path in tsv_paths) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    expect_gt(length(lines), 1L)
    expect_false(startsWith(lines[[1L]], "\ufeff"))
    field_counts <- lengths(strsplit(lines, "\t", fixed = TRUE))
    expect_length(unique(field_counts), 1L)
    table <- data.table::fread(
      path, sep = "\t", na.strings = "NA", data.table = FALSE,
      check.names = FALSE
    )
    expect_identical(anyDuplicated(names(table)), 0L)
    expect_true(all(nzchar(names(table))))
  }

  outcomes <- data.table::fread(
    file.path(directory, "phenotype-metadata-outcomes.tsv"),
    sep = "\t", na.strings = "NA", data.table = FALSE
  )
  predictors <- data.table::fread(
    file.path(directory, "phenotype-metadata-predictors.tsv"),
    sep = "\t", na.strings = "NA", data.table = FALSE
  )
  dictionary <- data.table::fread(
    file.path(directory, "analysis-data-dictionary.tsv"),
    sep = "\t", na.strings = "NA", data.table = FALSE
  )
  expect_identical(names(outcomes), c(
    "phenotype", "description", "group", "groupnum", "color",
    "variable_type", "outcome_type", "reference", "levels", "scores",
    "offset"
  ))
  expect_identical(names(predictors), c(
    "phenotype", "description", "group", "groupnum", "color",
    "variable_type", "reference", "levels", "scores"
  ))
  expect_true("representation" %in% names(dictionary))
  expect_setequal(
    dictionary$column,
    c(
      "id", "pgs", "endpoint_continuous", "age", "sex",
      "pheno_binary", "pheno_continuous", "pheno_count",
      "pheno_ordinal", "log_followup", "exclude_flag"
    )
  )
  canonical <- c(
    "phenotype", "description", "group", "groupnum", "color",
    "variable_type"
  )
  expect_true(all(c(canonical, "outcome_type") %in% names(outcomes)))
  expect_true(all(canonical %in% names(predictors)))
  expect_identical(anyDuplicated(outcomes$phenotype), 0L)
  expect_identical(anyDuplicated(predictors$phenotype), 0L)
  expect_true(all(vapply(
    c(outcomes$color, predictors$color),
    function(color) {
      !inherits(try(grDevices::col2rgb(color), silent = TRUE), "try-error")
    },
    logical(1L)
  )))
  expect_true(all(!is.na(outcomes$reference[outcomes$variable_type == "binary"])))
  expect_true(all(!is.na(outcomes$levels[outcomes$variable_type == "ordinal"])))
  expect_true(all(!is.na(predictors$scores[predictors$variable_type == "ordinal"])))
  for (metadata in list(outcomes, predictors)) {
    expect_true(all(vapply(
      split(metadata$color, metadata$group),
      function(value) length(unique(value)) == 1L,
      logical(1L)
    )))
    expect_true(all(vapply(
      split(metadata$groupnum, metadata$group),
      function(value) length(unique(value)) == 1L,
      logical(1L)
    )))
  }

  portable_text <- paste(
    unlist(lapply(file.path(directory, expected), readLines, warn = FALSE)),
    collapse = "\n"
  )
  expect_no_match(portable_text, "/Users/")
  expect_no_match(portable_text, "/home/")
  expect_no_match(portable_text, "[A-Za-z]:[/\\\\]")
  expect_no_match(portable_text, "run_fingerprint")
})

test_that("the real-data runner writes lean forward and reverse results", {
  template_dir <- real_data_example_directory()
  parent <- withr::local_tempdir()
  project <- file.path(parent, "PheWAS project with spaces")

  creator <- new.env(parent = globalenv())
  sys.source(file.path(template_dir, "create-demo-project.R"), envir = creator)
  creator$create_phewasflow_demo_project(project, template_dir = template_dir)
  expect_true(file.exists(file.path(project, "run-local.R")))
  expect_setequal(
    readLines(file.path(project, ".gitignore"), warn = FALSE),
    c("/inputs/", "/phewasflow-output/", ".Rhistory", ".RData")
  )
  expect_error(
    creator$create_phewasflow_demo_project(project, template_dir = template_dir),
    "already exists and is not empty"
  )

  analysis <- readRDS(file.path(project, "inputs", "analysis.rds"))
  expect_true(is.factor(analysis$sex))
  expect_identical(nrow(analysis), 220L)
  expect_identical(anyDuplicated(analysis$id), 0L)

  forward_path <- file.path(project, "analysis-forward.yml")
  reverse_path <- file.path(project, "analysis-reverse.yml")
  forward <- read_phewasflow_config(forward_path)
  reverse <- read_phewasflow_config(reverse_path)
  expect_identical(forward$analysis$direction, "phenotypes_as_outcomes")
  expect_identical(reverse$analysis$direction, "phenotypes_as_predictors")
  expect_identical(forward$plot$significance, "bonferroni")
  expect_identical(reverse$plot$significance, "bonferroni")
  expect_identical(forward$plot$group_display, "legend")
  expect_identical(reverse$plot$group_display, "legend")
  expect_true(all(startsWith(
    unlist(forward$inputs),
    normalizePath(project, winslash = "/", mustWork = TRUE)
  )))

  unrelated <- file.path(parent, "unrelated working directory")
  dir.create(unrelated)
  withr::local_dir(unrelated)
  runner <- new.env(parent = globalenv())
  sys.source(file.path(project, "run-local.R"), envir = runner)
  invisible(utils::capture.output(
    forward_result <- suppressMessages(runner$run_phewasflow_local(
      "analysis-forward.yml", project_dir = project
    ))
  ))

  output <- file.path(project, "phewasflow-output", "forward")
  result <- read_phewas_bundle(output)
  expect_equal(result, forward_result, ignore_attr = TRUE)
  expect_identical(nrow(result), 4L)
  expect_true(all(result$status == "ok"))
  expect_true(all(c(
    "phenotype", "description", "group", "groupnum", "color"
  ) %in% names(result)))
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label",
    "category", "category_order"
  ) %in% names(result)))
  tsv_result <- data.table::fread(
    file.path(output, "results.tsv"), sep = "\t", data.table = FALSE,
    check.names = FALSE
  )
  expect_identical(names(tsv_result), names(as.data.frame(result)))
  expect_setequal(result$engine, c("logistf", "lm", "glm.nb", "clm"))
  expect_setequal(
    list.files(output),
    c(
      "phewas.png", "phewas.pdf", "phewas-volcano.png",
      "phewas-volcano.pdf", "results.rds", "results.tsv"
    )
  )
  expect_setequal(
    list.files(file.path(output, ".phewasflow"), recursive = TRUE),
    c("manifest.tsv", "shards/shard-0001.rds", "shards/shard-0002.rds")
  )
  generated <- list.files(output, recursive = TRUE, all.files = TRUE)
  expect_false(any(grepl("sha256|checksums|diagnostics|COMPLETE|review", generated)))

  invisible(utils::capture.output(
    reverse_result <- suppressMessages(runner$run_phewasflow_local(
      "analysis-reverse.yml", project_dir = project
    ))
  ))
  reverse_output <- file.path(
    project, "phewasflow-output", "reverse-continuous"
  )
  expect_identical(nrow(reverse_result), 3L)
  expect_true(all(reverse_result$status == "ok"))
  expect_true(all(reverse_result$direction == "phenotypes_as_predictors"))
  expect_setequal(
    list.files(reverse_output),
    c(
      "phewas.png", "phewas.pdf", "phewas-volcano.png",
      "phewas-volcano.pdf", "results.rds", "results.tsv"
    )
  )
  expect_setequal(
    list.files(file.path(reverse_output, ".phewasflow"), recursive = TRUE),
    c("manifest.tsv", "shards/shard-0001.rds", "shards/shard-0002.rds")
  )
  reverse_generated <- list.files(
    reverse_output, recursive = TRUE, all.files = TRUE
  )
  expect_false(any(grepl(
    "sha256|checksums|diagnostics|COMPLETE|review", reverse_generated
  )))
})

test_that("the real-data runner completes cleanly when nothing is plottable", {
  template_dir <- real_data_example_directory()
  parent <- withr::local_tempdir()
  project <- file.path(parent, "all skipped project")
  creator <- new.env(parent = globalenv())
  sys.source(file.path(template_dir, "create-demo-project.R"), envir = creator)
  creator$create_phewasflow_demo_project(project, template_dir = template_dir)

  config_path <- file.path(project, "analysis-forward.yml")
  config_lines <- readLines(config_path, warn = FALSE)
  config_lines <- sub(
    "min_n: 100", "min_n: 1000", config_lines, fixed = TRUE
  )
  writeLines(config_lines, config_path, useBytes = TRUE)

  runner <- new.env(parent = globalenv())
  sys.source(file.path(project, "run-local.R"), envir = runner)
  output_text <- utils::capture.output(
    result <- runner$run_phewasflow_local(
      "analysis-forward.yml", project_dir = project
    ),
    type = "message"
  )

  expect_true(all(result$status == "skipped"))
  expect_true(all(result$reason_code == "insufficient_n"))
  expect_match(paste(output_text, collapse = "\n"), "plots were not created")
  expect_match(paste(output_text, collapse = "\n"), "ok=0, skipped=4, error=0")
  output <- file.path(project, "phewasflow-output", "forward")
  expect_setequal(list.files(output), c("results.rds", "results.tsv"))
  expect_false(any(file.exists(file.path(
    output,
    c("phewas.png", "phewas.pdf", "phewas-volcano.png", "phewas-volcano.pdf")
  ))))
})
