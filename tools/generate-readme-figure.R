#!/usr/bin/env Rscript

if (!file.exists("DESCRIPTION") || !dir.exists("R")) {
  stop("Run this script from the phewasFlow repository root.", call. = FALSE)
}

required_packages <- c("ggplot2", "ggrepel")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Install the missing package(s): ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  sys.source(source_file, envir = .GlobalEnv)
}

set.seed(20260819)

groups <- data.frame(
  group = c(
    "Cardiovascular", "Metabolic", "Respiratory", "Gastrointestinal",
    "Renal", "Hematologic", "Neurologic", "Dermatologic"
  ),
  color = c(
    "#0072B2", "#D55E00", "#009E73", "#E69F00",
    "#56B4E9", "#CC79A7", "#6F4C9B", "#666666"
  ),
  stringsAsFactors = FALSE
)

phenotypes_per_group <- 20L
n <- nrow(groups) * phenotypes_per_group
group_index <- rep(seq_len(nrow(groups)), each = phenotypes_per_group)
within_group <- rep(seq_len(phenotypes_per_group), times = nrow(groups))

result <- data.frame(
  analysis_id = "simulated_readme_example",
  phenotype = sprintf("SIM_%03d", seq_len(n)),
  description = sprintf(
    "%s outcome %02d", groups$group[group_index], within_group
  ),
  group = groups$group[group_index],
  groupnum = group_index,
  color = groups$color[group_index],
  status = "ok",
  direction = "phenotypes_as_outcomes",
  response = sprintf("SIM_%03d", seq_len(n)),
  predictor = "simulated_biomarker_score",
  outcome_type = "binary",
  predictor_type = "numeric",
  transformation = "zscore",
  reference = "0",
  outcome_reference = "0",
  predictor_reference = NA_character_,
  effect_measure = "odds_ratio",
  phenotype_order = within_group,
  stringsAsFactors = FALSE
)

result$std_error <- stats::runif(n, min = 0.055, max = 0.105)
z_value <- stats::rnorm(n)

signals <- data.frame(
  group = c(
    "Cardiovascular", "Metabolic", "Respiratory", "Renal",
    "Hematologic", "Dermatologic"
  ),
  description = c(
    "Atrial fibrillation", "Type 2 diabetes", "Asthma",
    "Chronic kidney disease", "Iron deficiency anemia", "Eczema"
  ),
  z_value = c(6.7, -7.1, -5.5, 6.1, -4.8, 4.5),
  stringsAsFactors = FALSE
)

signal_rows <- match(
  signals$group,
  result$group[!duplicated(result$group)]
)
signal_rows <- (signal_rows - 1L) * phenotypes_per_group + 1L
result$description[signal_rows] <- signals$description
z_value[signal_rows] <- signals$z_value

result$estimate <- z_value * result$std_error
result$p_value <- 2 * stats::pnorm(-abs(z_value))
result$q_value <- stats::p.adjust(result$p_value, method = "BH")
result$neg_log10_p <- -log10(result$p_value)
result$neg_log10_q <- -log10(result$q_value)
result$native_effect <- exp(result$estimate)
result$native_conf_low <- exp(result$estimate - 1.96 * result$std_error)
result$native_conf_high <- exp(result$estimate + 1.96 * result$std_error)
result$testing_family_size <- n

class(result) <- c("phewas_result", "data.frame")
attr(result, "spec") <- list(
  analysis_id = "simulated_readme_example",
  fdr_threshold = 0.05
)

labels <- result$phenotype[match(
  c(
    "Atrial fibrillation", "Type 2 diabetes", "Asthma",
    "Chronic kidney disease"
  ),
  result$description
)]

manhattan <- plot_phewas_manhattan(
  result,
  significance = "bonferroni",
  label = labels,
  label_mode = "description",
  label_wrap_width = 18L,
  group_legend_columns = 4L,
  point_size = 1.6,
  base_size = 9.5,
  label_size = 8.4,
  title = "Manhattan view"
)

volcano <- plot_phewas_volcano(
  result,
  significance = "bonferroni",
  label = labels,
  top_label_n = 0L,
  label_mode = "description",
  label_wrap_width = 18L,
  group_legend_columns = 4L,
  point_size = 1.6,
  base_size = 9.5,
  label_size = 8.4,
  title = "Volcano view"
)

extract_bottom_legend <- function(plot) {
  table <- ggplot2::ggplotGrob(plot)
  positions <- which(grepl("^guide-box", table$layout$name))
  positions <- positions[!vapply(
    table$grobs[positions], inherits, logical(1L), what = "zeroGrob"
  )]
  if (!length(positions)) {
    stop("The example plots did not produce a legend.", call. = FALSE)
  }
  table$grobs[[positions[[1L]]]]
}

output <- "man/figures/phewasflow-example.png"
grDevices::png(
  filename = output,
  width = 2600,
  height = 1120,
  res = 200,
  bg = "white"
)

shared_legend <- extract_bottom_legend(manhattan)
manhattan_table <- ggplot2::ggplotGrob(
  manhattan + ggplot2::theme(legend.position = "none")
)
volcano_table <- ggplot2::ggplotGrob(
  volcano + ggplot2::theme(legend.position = "none")
)
if (length(manhattan_table$heights) != length(volcano_table$heights) ||
    length(manhattan_table$widths) != length(volcano_table$widths)) {
  stop("The paired plot layouts could not be aligned.", call. = FALSE)
}
shared_heights <- do.call(
  grid::unit.pmax, list(manhattan_table$heights, volcano_table$heights)
)
shared_widths <- do.call(
  grid::unit.pmax, list(manhattan_table$widths, volcano_table$widths)
)
manhattan_table$heights <- volcano_table$heights <- shared_heights
manhattan_table$widths <- volcano_table$widths <- shared_widths

grid::grid.newpage()
layout <- grid::grid.layout(
  nrow = 2L,
  ncol = 2L,
  heights = grid::unit.c(
    grid::unit(1, "null"),
    sum(shared_legend$heights) + grid::unit(8, "pt")
  ),
  widths = grid::unit(c(1, 1), "null")
)
grid::pushViewport(grid::viewport(layout = layout))
grid::pushViewport(
  grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L)
)
grid::grid.draw(manhattan_table)
grid::popViewport()
grid::pushViewport(
  grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L)
)
grid::grid.draw(volcano_table)
grid::popViewport()
grid::pushViewport(
  grid::viewport(layout.pos.row = 2L, layout.pos.col = 1:2)
)
grid::grid.draw(shared_legend)
grid::popViewport()
grid::popViewport()
grDevices::dev.off()

message("Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE))
