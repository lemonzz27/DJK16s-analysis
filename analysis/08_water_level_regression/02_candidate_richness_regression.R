local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")

  .libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "Rlib"), .libPaths()))

  required_packages <- c("ggplot2", "dplyr", "ragg", "svglite", "scales", "sandwich", "lmtest")

  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0L) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2L) {
    stop("Usage: Rscript script.R <pathogen_input_dir> <output_root>")
  }


  pathogen_input_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)

  output_root <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)


  pathogen_output_dir <- output_root

  dir.create(pathogen_output_dir, recursive = TRUE, showWarnings = FALSE)

  normalize_sample <- function(x) {
    gsub("\\.", "_", trimws(as.character(x)))
  }

  read_csv_safely <- function(path, row.names = NULL) {
    value <- tryCatch(read.csv(path, row.names = row.names, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE,
      check.names = FALSE), error = function(e) NULL)
    if (is.null(value)) {
      value <- read.csv(path, row.names = row.names, fileEncoding = "GBK", stringsAsFactors = FALSE,
        check.names = FALSE)
    }
    value
  }

  read_abundance <- function(path) {
    value <- as.matrix(read_csv_safely(path, row.names = 1))
    storage.mode(value) <- "double"
    if (anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
      stop("Invalid abundance values in ", basename(path))
    }
    colnames(value) <- normalize_sample(colnames(value))
    value
  }

  read_metadata <- function(input_dir) {
    group_size <- read_csv_safely(file.path(input_dir, "groupsize.csv"))
    env <- read_csv_safely(file.path(input_dir, "env.csv"))
    names(group_size)[1:2] <- c("sample", "Raw_group")
    if (!all(c("sample", "WL") %in% names(env))) {
      stop("env.csv must contain sample and WL columns: ", input_dir)
    }
    group_size %>%
      transmute(sample = normalize_sample(sample), Raw_group = tolower(as.character(Raw_group))) %>%
      left_join(env %>%
        transmute(sample = normalize_sample(sample), Water_level = as.numeric(WL)), by = "sample") %>%
      mutate(Station = sub("^([^0-9]+).*$", "\\1", sample), month = sub("^[^0-9]*([0-9]{4}).*$", "\\1",
        sample), source_date = as.Date(paste0("20", month, "01"), format = "%Y%m%d"), Fraction = recode(Raw_group,
        large = "ConF", small = "USF")) %>%
      select(sample, Station, month, source_date, Water_level, Fraction)
  }

  audit_sample_design <- function(data, label) {
    if (nrow(data) != 120L)
      stop(label, ": expected 120 samples, found ", nrow(data))
    if (n_distinct(data$month) != 15L)
      stop(label, ": expected 15 months")
    if (n_distinct(data$Station) != 4L)
      stop(label, ": expected four stations")
    if (!setequal(data$Fraction, c("ConF", "USF"))) {
      stop(label, ": both ConF and USF are required")
    }
    per_month <- data %>%
      group_by(month) %>%
      summarise(n_samples = n(), n_stations = n_distinct(Station), n_fractions = n_distinct(Fraction),
        n_water_levels = n_distinct(Water_level), .groups = "drop")
    if (any(per_month$n_samples != 8L))
      stop(label, ": expected eight samples per month")
    if (any(per_month$n_stations != 4L))
      stop(label, ": expected four stations per month")
    if (any(per_month$n_fractions != 2L))
      stop(label, ": expected two fractions per month")
    if (any(per_month$n_water_levels != 1L)) {
      stop(label, ": water level differs within a month")
    }
    invisible(per_month)
  }

  make_analysis_data <- function(raw_data, metric_columns) {
    size_data <- raw_data %>%
      mutate(Group = Fraction)
    all_data <- raw_data %>%
      mutate(Group = "All")
    combined <- bind_rows(all_data, size_data) %>%
      mutate(Group = factor(Group, levels = c("All", "ConF", "USF")))
    bind_rows(lapply(names(metric_columns), function(metric_label) {
      value_column <- metric_columns[[metric_label]]
      combined %>%
        transmute(sample, Station, month, source_date, Water_level, Fraction, Group, Metric = metric_label,
          Value = as.numeric(.data[[value_column]]))
    })) %>%
      filter(is.finite(Value), is.finite(Water_level))
  }

  format_p <- function(p) {
    if (!is.finite(p))
      return("NA")
    if (p < 0.001)
      return("< 0.001")
    paste0("= ", sprintf("%.3f", p))
  }

  fit_one_group <- function(data) {
    data <- data %>%
      arrange(Station, Fraction, source_date)
    fit <- lm(Value ~ Water_level, data = data)
    n_months <- n_distinct(data$month)
    robust_vcov <- sandwich::vcovCL(fit, cluster = data$month, type = "HC1", cadjust = TRUE)
    robust_test <- lmtest::coeftest(fit, vcov. = robust_vcov, df = n_months - 1L)
    t_critical <- qt(0.975, df = n_months - 1L)
    slope <- unname(coef(fit)[["Water_level"]])
    slope_se <- robust_test["Water_level", "Std. Error"]
    x_grid <- seq(min(data$Water_level), max(data$Water_level), length.out = 180L)
    design <- model.matrix(~Water_level, data = data.frame(Water_level = x_grid))
    predicted <- as.vector(design %*% coef(fit))
    predicted_se <- sqrt(rowSums((design %*% robust_vcov) * design))
    prediction <- data.frame(Water_level = x_grid, fitted = predicted, lower = predicted - t_critical *
      predicted_se, upper = predicted + t_critical * predicted_se)
    shapiro <- shapiro.test(residuals(fit))
    bp <- lmtest::bptest(fit)
    standardized <- rstandard(fit)
    diagnostics <- data.frame(residual_Shapiro_W = unname(shapiro$statistic), residual_Shapiro_p = shapiro$p.value,
      Breusch_Pagan_statistic = unname(bp$statistic), Breusch_Pagan_p = bp$p.value, absolute_standardized_residual_gt_3 = sum(abs(standardized) >
        3), stringsAsFactors = FALSE)
    stats <- data.frame(n_observations = nrow(data), n_months = n_months, n_stations = n_distinct(data$Station),
      slope = slope, slope_cluster_SE = slope_se, slope_cluster_CI95_lower = slope - t_critical * slope_se,
      slope_cluster_CI95_upper = slope + t_critical * slope_se, intercept = unname(coef(fit)[["(Intercept)"]]),
      R2 = summary(fit)$r.squared, adjusted_R2 = summary(fit)$adj.r.squared, p_month_clustered = robust_test["Water_level",
        "Pr(>|t|)"], Pearson_r = cor(data$Water_level, data$Value), stringsAsFactors = FALSE)
    list(stats = stats, prediction = prediction, diagnostics = diagnostics)
  }

  fit_all_models <- function(plot_data) {
    keys <- plot_data %>%
      distinct(Metric, Group) %>%
      arrange(Metric, Group)
    stats_list <- vector("list", nrow(keys))
    prediction_list <- vector("list", nrow(keys))
    diagnostics_list <- vector("list", nrow(keys))
    for (index in seq_len(nrow(keys))) {
      metric_value <- keys$Metric[[index]]
      group_value <- as.character(keys$Group[[index]])
      current <- plot_data %>%
        filter(Metric == metric_value, as.character(Group) == group_value)
      fitted <- fit_one_group(current)
      stats_list[[index]] <- bind_cols(data.frame(Metric = metric_value, Group = group_value), fitted$stats)
      prediction_list[[index]] <- fitted$prediction %>%
        mutate(Metric = metric_value, Group = group_value)
      diagnostics_list[[index]] <- bind_cols(data.frame(Metric = metric_value, Group = group_value),
        fitted$diagnostics)
    }
    list(stats = bind_rows(stats_list) %>%
      mutate(Group = factor(Group, levels = c("All", "ConF", "USF"))), prediction = bind_rows(prediction_list) %>%
      mutate(Group = factor(Group, levels = c("All", "ConF", "USF"))), diagnostics = bind_rows(diagnostics_list) %>%
      mutate(Group = factor(Group, levels = c("All", "ConF", "USF"))))
  }

  colors <- c(All = "#2C7FB8", ConF = "#E67E22", USF = "#2A9D8F")

  fills <- colors

  shapes <- c(All = 23, ConF = 24, USF = 21)

  make_plot <- function(plot_data, fits, metric_name, y_label, output_dir, file_stem) {
    current_data <- plot_data %>%
      filter(Metric == metric_name)
    displayed_points <- current_data
    current_stats <- fits$stats %>%
      filter(Metric == metric_name) %>%
      arrange(Group) %>%
      mutate(sign = ifelse(intercept < 0, "-", "+"), label = paste0(as.character(Group), ": Y = ",
        sprintf("%.2f", slope), " X ", sign, " ", sprintf("%.2f", abs(intercept)), "; R2 = ", sprintf("%.3f",
          R2), "; P ", vapply(p_month_clustered, format_p, character(1))))
    current_prediction <- fits$prediction %>%
      filter(Metric == metric_name)
    x_range <- range(current_data$Water_level)
    y_range <- range(current_data$Value)
    x_span <- max(diff(x_range), 1)
    y_span <- max(diff(y_range), max(1, abs(y_range[[1]]) * 0.1))
    x_limits <- x_range + c(-0.08, 0.08) * x_span
    y_limits <- c(max(0, y_range[[1]] - 0.1 * y_span), y_range[[2]] + 0.42 * y_span)
    current_stats <- current_stats %>%
      mutate(x = x_limits[[1]] + 0.025 * diff(x_limits), y = y_limits[[2]] - c(0.05, 0.118, 0.186) *
        diff(y_limits))
    figure <- ggplot() + geom_ribbon(data = current_prediction, aes(x = Water_level, ymin = lower, ymax = upper,
      fill = Group), alpha = 0.12, color = NA, show.legend = FALSE) + geom_line(data = current_prediction,
      aes(x = Water_level, y = fitted, color = Group), linewidth = 1.05) + geom_point(data = displayed_points,
      aes(x = Water_level, y = Value, color = Group, fill = Group, shape = Group), size = 2.35, alpha = 0.68,
      stroke = 0.62) + geom_text(data = current_stats, aes(x = x, y = y, label = label, color = Group),
      hjust = 0, vjust = 1, size = 3.35, fontface = "bold", show.legend = FALSE) + scale_color_manual(values = colors,
      breaks = names(colors), drop = FALSE) + scale_fill_manual(values = fills, breaks = names(colors),
      drop = FALSE) + scale_shape_manual(values = shapes, breaks = names(colors), drop = FALSE) + scale_x_continuous(breaks = seq(150,
      170, 5), expand = expansion(mult = c(0.01, 0.01))) + scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0.01, 0.01))) + coord_cartesian(xlim = x_limits, ylim = y_limits,
      clip = "off") + labs(x = "Water level (m)", y = y_label, color = "Model", fill = "Model", shape = "Model") +
      theme_classic(base_size = 11, base_family = "sans") + theme(axis.line = element_blank(), panel.border = element_rect(fill = NA,
      color = "black", linewidth = 0.7), panel.grid = element_blank(), axis.ticks = element_line(color = "black",
      linewidth = 0.45), axis.ticks.length = grid::unit(0.13, "cm"), axis.text = element_text(color = "black",
      size = 10), axis.title = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom",
      legend.direction = "horizontal", legend.title = element_text(face = "bold", size = 9.5), legend.text = element_text(size = 9.5),
      legend.key = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "white",
        color = NA), plot.margin = margin(7, 8, 5, 7)) + guides(color = guide_legend(override.aes = list(shape = unname(shapes),
      fill = unname(fills), alpha = 1, linewidth = 1.05, size = 3), nrow = 1), fill = "none", shape = "none")
    prefix <- file.path(output_dir, file_stem)
    ggsave(paste0(prefix, ".png"), figure, width = 7.4, height = 5.9, dpi = 600, bg = "white", device = ragg::agg_png)
    ggsave(paste0(prefix, ".svg"), figure, width = 7.4, height = 5.9, bg = "white", device = svglite::svglite)
    ggsave(paste0(prefix, ".pdf"), figure, width = 7.4, height = 5.9, bg = "white", device = cairo_pdf)
    ggsave(paste0(prefix, ".tiff"), figure, width = 7.4, height = 5.9, dpi = 600, bg = "white", compression = "lzw",
      device = ragg::agg_tiff)
    figure
  }

  write_analysis_outputs <- function(raw_data, plot_data, fits, output_dir, prefix, interpretation_lines) {
    write.csv(raw_data, file.path(output_dir, paste0(prefix, "_raw_sample_data.csv")), row.names = FALSE,
      fileEncoding = "UTF-8")
    write.csv(plot_data, file.path(output_dir, paste0(prefix, "_plot_data.csv")), row.names = FALSE,
      fileEncoding = "UTF-8")
    write.csv(fits$stats, file.path(output_dir, paste0(prefix, "_regression_results.csv")), row.names = FALSE,
      fileEncoding = "UTF-8")
    write.csv(fits$diagnostics, file.path(output_dir, paste0(prefix, "_regression_diagnostics.csv")),
      row.names = FALSE, fileEncoding = "UTF-8")
    writeLines(c(interpretation_lines, "Markers retain the original convention: All = blue diamonds, ConF = orange triangles, USF = teal circles.",
      "The All model displays and pools 120 observations; ConF and USF each display and use 60.", "OLS slopes and R2 are descriptive.",
      "P values and 95% confidence bands use HC1 standard errors clustered by month (15 clusters).",
      "Raw large/small codes are retained only for input compatibility; figures display ConF/USF."),
      file.path(output_dir, paste0(prefix, "_analysis_notes.txt")), useBytes = TRUE)
  }

  pathogen_metadata <- read_metadata(pathogen_input_dir)

  community_path <- file.path(dirname(dirname(pathogen_input_dir)), "tpm_rarefied.csv")

  community <- read_abundance(community_path)

  candidate_input <- read_abundance(file.path(pathogen_input_dir, "04_pathogen_candidate_abundance.csv"))

  candidate_ids <- rownames(candidate_input)

  if (length(candidate_ids) != 2755L) {
    stop("Expected exactly 2,755 fixed candidate ASVs; found ", length(candidate_ids))
  }

  matched_candidate_ids <- intersect(candidate_ids, rownames(community))

  if (length(matched_candidate_ids) != length(candidate_ids)) {
    stop("Candidate ASV matching failed: ", length(matched_candidate_ids), " of ", length(candidate_ids),
      " IDs matched")
  }

  community_library <- colSums(community)

  if (length(community_library) != 120L || any(community_library != 21000)) {
    stop("Expected 120 whole-community samples each rarefied to 21,000 reads")
  }

  candidate_counts <- community[matched_candidate_ids, , drop = FALSE]

  pathogen_reads <- colSums(candidate_counts)

  pathogen_richness <- colSums(candidate_counts > 0)

  pathogen_relative_abundance <- 100 * pathogen_reads/community_library

  pathogen_raw <- data.frame(sample = colnames(community), Pathogen_richness = unname(pathogen_richness[colnames(community)]),
    Pathogen_reads = unname(pathogen_reads[colnames(community)]), Community_reads = unname(community_library[colnames(community)]),
    Pathogen_relative_abundance = unname(pathogen_relative_abundance[colnames(community)]), stringsAsFactors = FALSE) %>%
    left_join(pathogen_metadata, by = "sample") %>%
    arrange(source_date, Station, Fraction)

  if (anyNA(pathogen_raw)) stop("Pathogen metadata matching produced missing values")

  audit_sample_design(pathogen_raw, "Candidate-pathogen metrics")

  pathogen_plot_data <- make_analysis_data(pathogen_raw, c(`Pathogen richness` = "Pathogen_richness"))

  pathogen_fits <- fit_all_models(pathogen_plot_data)

  make_plot(pathogen_plot_data, pathogen_fits, "Pathogen richness", "Pathogen richness", pathogen_output_dir,
    "Pathogen_richness_vs_water_level_sample_level_All_ConF_USF")

  write_analysis_outputs(pathogen_raw, pathogen_plot_data, pathogen_fits, pathogen_output_dir, "Candidate_pathogen_sample_level",
    c("Figure S15: candidate-pathogen richness versus water level.", paste0("Candidate ASVs: ", length(candidate_ids),
      "."), "Richness = number of fixed candidate ASVs detected (>0 reads) in each whole-community sample rarefied once to 21,000 reads; no second candidate-subset rarefaction.",
      "Relative abundance = candidate-ASV reads / complete-community reads * 100%."))

  message("Completed Figure S15 single-rarefaction sample-level water regressions: ", output_root)
})
