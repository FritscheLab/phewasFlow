phewas_example_environment <- function(script) {
  path <- system.file("examples", "phewas", script, package = "phewasFlow")
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment
}

test_that("PheWAS phenotype definitions survive the bridge", {
  skip_if_not_installed("PheWAS")
  withr::local_options(lifecycle_verbosity = "quiet")
  example <- phewas_example_environment("use-phewas-phenotypes.R")
  bridge <- suppressMessages(example$run_phewas_bridge())
  ids <- bridge$metadata$phenotype
  expect_true(any(grepl(".", ids, fixed = TRUE)))
  expect_identical(ids, setdiff(names(bridge$phenotypes), "id"))
  expect_true(all(bridge$metadata$reference == "FALSE"))
  expect_true(any(bridge$result$status == "ok"))
  expect_true(all(bridge$result$status %in% c("ok", "skipped")))
  # Upstream exclusion mappings also produce columns with no observed cases.
  skipped <- bridge$result[bridge$result$status == "skipped", ]
  expect_gt(nrow(skipped), 0L)
  expect_true(all(skipped$reason_code == "insufficient_cases"))
  expect_true(all(skipped$n_event == 0L))
  location <- match(bridge$analysis$id, bridge$phenotypes$id)
  for (id in ids) {
    expected <- bridge$phenotypes[[id]][location]
    expect_identical(bridge$analysis[[id]], expected)
    expect_true(anyNA(expected))
    result <- bridge$result[bridge$result$phenotype == id, ]
    expect_equal(result$n_complete, sum(!is.na(expected)))
    expect_equal(result$n_event, sum(expected, na.rm = TRUE))
    expect_equal(result$n_nonevent, sum(!expected, na.rm = TRUE))
  }
})

test_that("shared Firth and linear scans agree with PheWAS", {
  skip_if_not_installed("PheWAS")
  withr::local_options(lifecycle_verbosity = "quiet")
  example <- phewas_example_environment("compare-packages.R")
  fixture <- example$make_phewas_comparison(n = 320L, per_type = 2L)
  for (type in names(fixture$specs)) {
    flow <- example$fit_phewas_comparison(fixture, type, "phewasFlow")
    upstream <- suppressMessages(example$fit_phewas_comparison(fixture, type, "PheWAS"))
    agreement <- example$check_phewas_agreement(flow, upstream, type)
    expect_equal(nrow(agreement), 2L)
    expect_true(all(agreement$pass), info = paste("Outcome type:", type))
    expect_true(all(agreement$n_flow < nrow(fixture$data)))
    # A changed estimate must not be reported as numerical agreement.
    upstream$beta[[1L]] <- upstream$beta[[1L]] + 0.01
    changed <- example$check_phewas_agreement(flow, upstream, type)
    expect_false(all(changed$pass))
    expect_error(example$check_phewas_agreement(flow, upstream[-1L, ], type),
                 "different or duplicate phenotype sets")
  }
})

test_that("the comparison exposes loss of precision in upstream Firth p-values", {
  skip_if_not_installed("PheWAS")
  # This regression documents the pinned upstream and underlying engine pair.
  skip_if(as.character(utils::packageVersion("logistf")) != "1.26.1")
  skip_if(!identical(utils::packageDescription("PheWAS")$RemoteSha,
                     "55dd1c24e228851922400cfba8d7db474565ccc7"))
  withr::local_options(lifecycle_verbosity = "quiet")
  example <- phewas_example_environment("compare-packages.R")
  fixture <- example$make_phewas_comparison(n = 5000L, per_type = 4L)
  flow <- example$fit_phewas_comparison(fixture, "binary", "phewasFlow")
  upstream <- suppressMessages(example$fit_phewas_comparison(fixture, "binary", "PheWAS"))
  agreement <- example$check_phewas_agreement(flow, upstream, "binary")
  expect_true(all(agreement$counts_match))
  expect_true(all(agreement$beta_abs_difference <= 1e-6))
  expect_true(all(agreement$se_abs_difference <= 1e-6))
  expect_true(all(is.finite(agreement$log_p_flow)))
  expect_true(any(agreement$p_comparison == "upstream_reported_zero"))
  expect_false(all(agreement$pass))
})
