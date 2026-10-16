local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 5L) {
    stop("Usage: Rscript script.R <alpha_csv> <group_csv> <abundance_csv> ", "<figure_dir> <table_dir>")
  }

  alpha_csv <- args[[1]]

  group_csv <- args[[2]]

  abundance_csv <- args[[3]]

  figure_dir <- args[[4]]

  table_dir <- args[[5]]

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  required_packages <- c("ggplot2", "dplyr", "patchwork", "ragg", "svglite")

  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0L) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
  })

  group_levels <- c("large", "small")

  group_labels <- c(large = "ConF\nn = 60", small = "USF\nn = 60")

  group_colours <- c(large = "#E67E22", small = "#2AA58D")

  fraction_colours <- c(Pathogen = "#D86648", `Non-pathogen` = "#A7B6C2")

  alpha <- read.csv(alpha_csv, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  groups <- read.csv(group_csv, check.names = FALSE, stringsAsFactors = FALSE)

  abundance <- read.csv(abundance_csv, check.names = FALSE, stringsAsFactors = FALSE)

  if (!"Richness" %in% names(alpha)) {
    stop("The alpha-diversity table does not contain a Richness column.")
  }

  if (!all(c("sample", "group") %in% names(groups))) {
    stop("The group table must contain sample and group columns.")
  }

  if (!all(c("group", "category", "mean_relative_abundance_percent") %in% names(abundance))) {
    stop("The abundance table is missing required columns.")
  }

  alpha_plot_data <- data.frame(sample = rownames(alpha), Richness = as.numeric(alpha$Richness), stringsAsFactors = FALSE) %>%
    left_join(groups[, c("sample", "group")], by = "sample") %>%
    filter(group %in% group_levels, is.finite(Richness)) %>%
    mutate(group = factor(tolower(group), levels = group_levels))

  if (anyNA(alpha_plot_data$group)) {
    stop("Some alpha-diversity samples could not be assigned to large/small groups.")
  }

  group_n <- table(alpha_plot_data$group)

  if (!identical(as.integer(group_n), c(60L, 60L))) {
    stop("Expected 60 samples per group; observed: ", paste(names(group_n), as.integer(group_n), collapse = ", "))
  }

  alpha_fit <- stats::aov(Richness ~ group, data = alpha_plot_data)

  alpha_tukey <- stats::TukeyHSD(alpha_fit, "group")$group

  alpha_p_adjusted <- unname(alpha_tukey[1, "p adj"])

  alpha_letters <- if (alpha_p_adjusted < 0.05) {
    c(large = "b", small = "a")
  } else {
    c(large = "a", small = "a")
  }

  set.seed(20260811)

  alpha_plot_data <- alpha_plot_data %>%
    mutate(x = as.numeric(group), x_jitter = x + stats::runif(n(), min = -0.072, max = 0.072))

  violin_data <- bind_rows(lapply(group_levels, function(group_name) {
    values <- alpha_plot_data$Richness[as.character(alpha_plot_data$group) == group_name]
    density_i <- stats::density(values, from = min(values), to = max(values), n = 512, adjust = 1)
    group_x <- match(group_name, group_levels)
    half_width <- 0.34 * density_i$y/max(density_i$y)
    data.frame(group = factor(group_name, levels = group_levels), x = c(group_x, group_x - half_width,
      group_x), Richness = c(density_i$x[1], density_i$x, density_i$x[length(density_i$x)]))
  }))

  richness_range <- diff(range(alpha_plot_data$Richness))

  stat_labels <- alpha_plot_data %>%
    group_by(group) %>%
    summarise(y = max(Richness) + 0.075 * richness_range, .groups = "drop") %>%
    mutate(x = as.numeric(group), label = unname(alpha_letters[as.character(group)]))

  abundance_plot_data <- abundance %>%
    transmute(group = factor(tolower(trimws(group)), levels = group_levels), category = ifelse(trimws(category) ==
      "Pathogen", "Pathogen", "Non-pathogen"), value = as.numeric(mean_relative_abundance_percent)) %>%
    filter(!is.na(group), is.finite(value)) %>%
    mutate(category = factor(category, levels = c("Pathogen", "Non-pathogen")), x = as.numeric(group))

  if (nrow(abundance_plot_data) != 4L) {
    stop("Expected four group-by-category abundance rows.")
  }

  abundance_sums <- abundance_plot_data %>%
    group_by(group) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    mutate(x = as.numeric(group))

  if (any(abs(abundance_sums$total - 100) > 1e-06)) {
    stop("Relative-abundance fractions do not sum to 100% within each group.")
  }

  theme_cns_closed <- function(base_size = 7.2) {
    theme_classic(base_size = base_size, base_family = "Arial") + theme(panel.border = element_rect(colour = "#202020",
      fill = NA, linewidth = 0.48), axis.line = element_blank(), axis.ticks = element_line(colour = "#202020",
      linewidth = 0.38), axis.ticks.length = grid::unit(1.5, "mm"), axis.text = element_text(colour = "#202020",
      size = 6.8), axis.title = element_text(colour = "#202020", size = 7.4), panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white",
        colour = NA), legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent",
        colour = NA), legend.key = element_rect(fill = "transparent", colour = NA))
  }

  p_richness <- ggplot(alpha_plot_data, aes(x = x, y = Richness, group = group, fill = group)) + geom_polygon(data = violin_data,
    aes(x = x, y = Richness, group = group, fill = group), inherit.aes = FALSE, alpha = 0.35, colour = NA) +
    geom_boxplot(width = 0.115, outlier.shape = NA, alpha = 0.86, colour = "#202020", linewidth = 0.38) +
    geom_point(aes(x = x_jitter), shape = 21, size = 1.45, alpha = 0.82, stroke = 0.2, colour = "white") +
    geom_text(data = stat_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, family = "Arial",
      fontface = "bold", size = 3.1, colour = "#202020") + scale_fill_manual(values = group_colours,
    guide = "none") + scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), expand = expansion(mult = c(0,
    0))) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.11))) + labs(x = NULL, y = "Richness") +
    theme_cns_closed() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(t = 4.5,
    r = 5, b = 0, l = 4, unit = "pt"))

  p_abundance <- ggplot(abundance_plot_data, aes(x = x, y = value, fill = category)) + geom_col(width = 0.38,
    colour = "white", linewidth = 0.45) + geom_text(aes(label = sprintf("%.1f%%", value), colour = category),
    position = position_stack(vjust = 0.5), family = "Arial", fontface = "bold", size = 2.55, show.legend = FALSE) +
    scale_fill_manual(values = fraction_colours, breaks = c("Pathogen", "Non-pathogen"), labels = c("Pathogen candidate",
      "Other ASVs"), name = NULL) + scale_colour_manual(values = c(Pathogen = "white", `Non-pathogen` = "#1F2933"),
    guide = "none") + scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = unname(group_labels[group_levels]),
    expand = expansion(mult = c(0, 0))) + scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = paste0(c(0,
    25, 50, 75, 100), "%"), expand = expansion(mult = c(0, 0))) + coord_cartesian(ylim = c(0, 104), expand = FALSE) +
    labs(x = NULL, y = "Mean relative abundance") + guides(fill = guide_legend(direction = "horizontal",
    override.aes = list(colour = NA))) + theme_cns_closed() + theme(axis.text.x = element_text(colour = "#202020",
    size = 7, lineheight = 0.92, margin = margin(t = 3)), legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 6.4, colour = "#202020"), legend.spacing.x = grid::unit(1.8, "mm"),
    legend.key.width = grid::unit(4.5, "mm"), legend.key.height = grid::unit(2.6, "mm"), plot.margin = margin(t = 0,
      r = 5, b = 1, l = 4, unit = "pt"))

  combined_plot <- p_richness/p_abundance + plot_layout(heights = c(1.52, 0.92), guides = "collect") &
    theme(legend.position = "bottom")

  width_mm <- 78

  height_mm <- 132

  width_in <- width_mm/25.4

  height_in <- height_mm/25.4

  prefix <- file.path(figure_dir, "Richness_pathogen_abundance_vertical_composite_CNS")

  ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in, units = "in", res = 600,
    background = "white")

  print(combined_plot)

  grDevices::dev.off()

  ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = 600,
    background = "white", compression = "lzw")

  print(combined_plot)

  grDevices::dev.off()

  svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")

  print(combined_plot)

  grDevices::dev.off()

  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in, family = "Arial",
    bg = "white")

  print(combined_plot)

  grDevices::dev.off()

  write.csv(alpha_plot_data, file.path(table_dir, "Richness_large_small_source_data.csv"), row.names = FALSE,
    fileEncoding = "UTF-8")

  write.csv(abundance_plot_data, file.path(table_dir, "pathogen_abundance_large_small_source_data.csv"),
    row.names = FALSE, fileEncoding = "UTF-8")

  write.csv(stat_labels, file.path(table_dir, "Richness_Tukey_letters.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  expected_outputs <- paste0(prefix, c(".png", ".tiff", ".svg", ".pdf"))

  stopifnot(all(file.exists(expected_outputs)), all(file.info(expected_outputs)$size > 0), identical(levels(alpha_plot_data$group),
    group_levels), identical(levels(abundance_plot_data$group), group_levels), all(abs(abundance_sums$total -
    100) < 1e-06))

  qa_file <- file.path(table_dir, "Richness_pathogen_abundance_composite_QA.txt")

  capture.output({
    cat("Vertical CNS-style composite QA\n")
    cat("================================\n")
    cat("Final size: 78 x 132 mm; raster resolution: 600 dpi\n")
    cat("Shared x order: Large, Small\n")
    cat("Sample sizes: Large n=60; Small n=60\n")
    cat("Both panels use closed borders, aligned x scales and zero inter-panel gap.\n")
    cat("Upper x labels removed; lower panel supplies the shared x axis.\n")
    cat("Redundant 100% labels above the stacked bars were removed.\n")
    cat("The bar tops end at 100%, leaving 4% headroom below the panel border.\n")
    cat("Richness Tukey letters: ", paste(names(alpha_letters), alpha_letters, collapse = "; "), "\n",
      sep = "")
    cat("Abundance totals: ", paste(abundance_sums$group, abundance_sums$total, collapse = "; "), "\n",
      sep = "")
    cat("Outputs: PNG, TIFF, SVG and PDF\n")
  }, file = qa_file)

  cat("Completed: ", normalizePath(prefix, winslash = "/", mustWork = FALSE), "\n", sep = "")
})
