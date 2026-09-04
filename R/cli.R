# Command-line interface -----------------------------------------------------

.pf_cli_usage <- function() {
  paste(
    "phewasflow <command> [options]",
    "",
    "Commands:",
    "  validate   Validate YAML, inputs, metadata, and the analysis specification",
    "  manifest   Write a deterministic shard manifest",
    "  run-shard  Run one manifest shard",
    "  combine    Combine all shards and write the results",
    "  plot       Create Manhattan and volcano plots in PNG and PDF",
    sep = "\n"
  )
}

.pf_cli_common_options <- function() {
  list(optparse::make_option(
    c("-c", "--config"),
    type = "character",
    help = "Path to the phewasFlow YAML configuration"
  ))
}

.pf_cli_parse <- function(command, args) {
  options <- .pf_cli_common_options()
  description <- switch(
    command,
    validate = "Validate a phewasFlow run configuration.",
    manifest = "Create the deterministic shard manifest.",
    `run-shard` = "Run one shard and atomically publish its result.",
    combine = "Combine completed shards.",
    plot = "Plot combined PheWAS results.",
    stop("Unknown command `", command, "`.\n\n", .pf_cli_usage(), call. = FALSE)
  )
  if (identical(command, "manifest")) {
    options <- c(options, list(optparse::make_option(
      c("-o", "--output"),
      type = "character",
      default = NULL,
      help = "Manifest TSV path [default: <output.directory>/.phewasflow/manifest.tsv]"
    )))
  }
  if (identical(command, "run-shard")) {
    options <- c(options, list(
      optparse::make_option(
        c("-s", "--shard-id"),
        dest = "shard_id",
        type = "integer",
        default = suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))),
        help = "Shard number [default: SLURM_ARRAY_TASK_ID]"
      ),
      optparse::make_option(
        c("-m", "--manifest"),
        type = "character",
        default = NULL,
        help = "Manifest TSV path"
      ),
      optparse::make_option(
        "--overwrite",
        action = "store_true",
        default = FALSE,
        help = "Replace an existing shard artifact"
      )
    ))
  }
  if (identical(command, "combine")) {
    options <- c(options, list(
      optparse::make_option(
        c("-m", "--manifest"),
        type = "character",
        default = NULL,
        help = "Manifest TSV path"
      ),
      optparse::make_option(
        c("-b", "--bundle"),
        type = "character",
        default = NULL,
        help = "Result directory [default: <output.directory>]"
      ),
      optparse::make_option(
        "--overwrite",
        action = "store_true",
        default = FALSE,
        help = "Replace existing results"
      )
    ))
  }
  if (identical(command, "plot")) {
    options <- c(options, list(
      optparse::make_option(
        c("-b", "--bundle"),
        type = "character",
        default = NULL,
        help = "Result directory [default: <output.directory>]"
      ),
      optparse::make_option(
        c("-o", "--output"),
        type = "character",
        default = NULL,
        help = paste(
          "PNG/PDF filename or stem; writes Manhattan and volcano plots",
          "in both formats [default: plot.output]"
        )
      ),
      optparse::make_option(
        "--significance",
        type = "character",
        default = NULL,
        help = "Threshold lines: bh, bonferroni, or none [default: plot.significance]"
      ),
      optparse::make_option(
        "--group-display",
        dest = "group_display",
        type = "character",
        default = NULL,
        help = "Group labels: auto, x_axis, or legend [default: plot.group_display]"
      )
    ))
  }
  parser <- optparse::OptionParser(description = description, option_list = options)
  parsed <- optparse::parse_args(parser, args = args, positional_arguments = FALSE)
  if (is.null(parsed$config) || !nzchar(parsed$config)) {
    stop("`--config` is required.", call. = FALSE)
  }
  parsed
}

.pf_cli_print_json <- function(x) {
  cat(jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  ), "\n", sep = "")
}

#' Run the phewasFlow command-line interface
#'
#' This function powers the installed `phewasflow` executable. It is also
#' useful for testing or embedding the same commands in another R process.
#'
#' @param args Character vector of command-line arguments.
#'
#' @return Zero, invisibly, on success. Errors are signalled as R conditions.
#' @export
phewasflow_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[[1L]] %in% c("help", "--help", "-h")) {
    cat(.pf_cli_usage(), "\n", sep = "")
    return(invisible(0L))
  }
  command <- args[[1L]]
  options <- .pf_cli_parse(command, args[-1L])

  if (identical(command, "validate")) {
    config <- read_phewasflow_config(options$config, check_data = TRUE)
    metadata <- .pf_read_phenotype_metadata(
      config$inputs$phenotype_metadata
    )
    .pf_cli_print_json(list(
      status = "valid",
      analysis_id = config$analysis$analysis_id,
      direction = config$analysis$direction,
      phenotype_count = nrow(metadata),
      shards = config$execution$shards
    ))
    return(invisible(0L))
  }

  config <- read_phewasflow_config(options$config, check_data = FALSE)
  if (identical(command, "manifest")) {
    manifest <- create_phewas_manifest(config, path = options$output, write = TRUE)
    .pf_cli_print_json(list(
      status = "created",
      path = attr(manifest, "manifest_path", exact = TRUE),
      phenotype_count = nrow(manifest),
      shard_count = length(unique(manifest$shard_id)),
      run_fingerprint = unique(manifest$run_fingerprint)
    ))
    return(invisible(0L))
  }
  if (identical(command, "run-shard")) {
    if (length(options$shard_id) != 1L || is.na(options$shard_id)) {
      stop("`--shard-id` is required when SLURM_ARRAY_TASK_ID is not set.", call. = FALSE)
    }
    path <- run_phewas_shard(
      config,
      shard_id = options$shard_id,
      manifest_path = options$manifest,
      overwrite = options$overwrite
    )
    .pf_cli_print_json(list(status = "complete", shard_path = path))
    return(invisible(0L))
  }
  if (identical(command, "combine")) {
    result <- combine_phewas_shards(
      config,
      manifest_path = options$manifest,
      bundle_directory = options$bundle,
      overwrite = options$overwrite
    )
    .pf_cli_print_json(list(
      status = "complete",
      results = attr(result, "bundle_directory", exact = TRUE),
      association_count = nrow(result)
    ))
    return(invisible(0L))
  }
  if (identical(command, "plot")) {
    bundle <- .pf_default(options$bundle, config$output$directory)
    output <- .pf_default(options$output, config$plot$output)
    output_paths <- .pf_plot_output_paths(output, name = "--output")
    configured_output_paths <- .pf_plot_output_paths(config$plot$output)
    allowed_paths <- if (setequal(
      .pf_path_key(output_paths), .pf_path_key(configured_output_paths)
    )) {
      configured_output_paths
    } else {
      character()
    }
    .pf_assert_safe_write_target(
      config,
      output_paths,
      purpose = "Plot",
      allow = allowed_paths
    )
    state_directory <- file.path(config$output$directory, ".phewasflow")
    if (any(vapply(
      output_paths,
      .pf_path_is_within,
      logical(1L),
      directory = state_directory
    ))) {
      stop("Plot output cannot be inside the internal `.phewasflow` directory.",
           call. = FALSE)
    }
    result <- read_phewas_bundle(bundle)
    significance <- .pf_default(options$significance, config$plot$significance)
    if (!significance %in% c("bh", "bonferroni", "none")) {
      stop("`--significance` must be bh, bonferroni, or none.", call. = FALSE)
    }
    group_display <- .pf_default(
      options$group_display, config$plot$group_display
    )
    if (!group_display %in% c("auto", "x_axis", "legend")) {
      stop("`--group-display` must be auto, x_axis, or legend.", call. = FALSE)
    }
    labels <- if (length(config$plot$label)) config$plot$label else NULL
    highlights <- if (length(config$plot$highlight)) {
      config$plot$highlight
    } else {
      NULL
    }
    manhattan <- plot_phewas_manhattan(
      result,
      significance = significance,
      group_display = group_display,
      label = labels,
      highlight = highlights
    )
    volcano <- plot_phewas_volcano(
      result,
      significance = significance,
      label = labels,
      highlight = highlights
    )
    graphs <- list(
      manhattan_png = manhattan,
      manhattan_pdf = manhattan,
      volcano_png = volcano,
      volcano_pdf = volcano
    )
    for (name in names(output_paths)) {
      save_phewas_plot(
        graphs[[name]],
        filename = output_paths[[name]],
        width = config$plot$width,
        height = config$plot$height,
        dpi = config$plot$dpi
      )
    }
    .pf_cli_print_json(list(
      status = "complete",
      plot = unname(output_paths[["manhattan_png"]]),
      plots = as.list(output_paths)
    ))
    return(invisible(0L))
  }
  stop("Unknown command `", command, "`.", call. = FALSE)
}
