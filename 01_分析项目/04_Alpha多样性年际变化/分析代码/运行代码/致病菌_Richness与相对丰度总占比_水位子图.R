rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop("Usage: Rscript script.R <input_dir> <output_dir>")
}

input_dir <- args[[1]]

output_dir <- args[[2]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "dplyr", "patchwork", "mgcv", "nlme", "ragg", "svglite", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(20260811)

required_files <- c("community_rarefied.csv", "pathogen_rarefied.csv", "groupmonth.csv", "groupsize.csv",
  "grouplocation.csv", "WaterLevel_Monthly_Data.csv", "GAMM_WaterLevel_Predictions.csv", "GAMM_WaterLevel_Model_Statistics.csv")

required_paths <- file.path(input_dir, required_files)

if (any(!file.exists(required_paths))) {
  stop("Missing input files: ", paste(required_files[!file.exists(required_paths)], collapse = ", "))
}

normalize_sample <- function(x) {
  gsub("\\.", "_", trimws(as.character(x)))
}

parse_month <- function(x) {
  value <- gsub("[^0-9]", "", as.character(x))
  result <- rep(as.Date(NA), length(value))
  is_four <- nchar(value) == 4L
  is_six <- nchar(value) == 6L
  result[is_four] <- as.Date(paste0("20", value[is_four], "01"), format = "%Y%m%d")
  result[is_six] <- as.Date(paste0(value[is_six], "01"), format = "%Y%m%d")
  result
}

format_p <- function(p_value) {
  if (length(p_value) == 0L || !is.finite(p_value))
    return("NA")
  if (p_value < 0.001)
    return("<0.001")
  sprintf("=%.3f", p_value)
}

read_abundance <- function(path) {
  value <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
  storage.mode(value) <- "double"
  if (anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
    stop("Invalid abundance values in ", basename(path))
  }
  colnames(value) <- normalize_sample(colnames(value))
  value
}

community <- read_abundance(file.path(input_dir, "community_rarefied.csv"))

pathogen <- read_abundance(file.path(input_dir, "pathogen_rarefied.csv"))

if (!setequal(colnames(community), colnames(pathogen))) {
  stop("Community and pathogen tables do not contain the same samples.")
}

pathogen <- pathogen[, colnames(community), drop = FALSE]

community_library <- colSums(community)

if (any(community_library <= 0)) {
  stop("At least one complete-community sample has zero total abundance.")
}

group_month <- read.csv(file.path(input_dir, "groupmonth.csv"), check.names = FALSE, stringsAsFactors = FALSE)

group_size <- read.csv(file.path(input_dir, "groupsize.csv"), check.names = FALSE, stringsAsFactors = FALSE)

group_location <- read.csv(file.path(input_dir, "grouplocation.csv"), check.names = FALSE, stringsAsFactors = FALSE)

names(group_month)[1:2] <- c("sample", "month")

names(group_size)[1:2] <- c("sample", "size")

names(group_location)[1:2] <- c("sample", "site")

group_month$sample <- normalize_sample(group_month$sample)

group_size$sample <- normalize_sample(group_size$sample)

group_location$sample <- normalize_sample(group_location$sample)

sample_metrics <- data.frame(sample = colnames(community), Pathogen_richness = as.numeric(colSums(pathogen >
  0)), Pathogen_relative_abundance = as.numeric(100 * colSums(pathogen)/community_library), Community_library_size = as.numeric(community_library),
  Pathogen_total_reads = as.numeric(colSums(pathogen)), stringsAsFactors = FALSE) %>%
  left_join(group_month[, c("sample", "month")], by = "sample") %>%
  left_join(group_size[, c("sample", "size")], by = "sample") %>%
  left_join(group_location[, c("sample", "site")], by = "sample") %>%
  mutate(Date = parse_month(month), cluster = interaction(site, size, drop = TRUE)) %>%
  filter(!is.na(Date), is.finite(Pathogen_richness), is.finite(Pathogen_relative_abundance), !is.na(cluster))

if (nrow(sample_metrics) != ncol(community)) {
  stop("Some samples are missing month, site, or size metadata.")
}

if (anyDuplicated(paste(sample_metrics$cluster, sample_metrics$Date))) {
  stop("Duplicate dates occur within a site-by-size repeated-measure cluster.")
}

fit_overall_gamm <- function(data, response, transformation = "identity", k = 6L) {
  model_data <- data.frame(response_raw = as.numeric(data[[response]]), Date = as.Date(data$Date),
    cluster = factor(data$cluster), stringsAsFactors = FALSE)
  model_data$time_num <- as.numeric(model_data$Date - min(model_data$Date))
  if (transformation == "logit_percent") {
    proportion <- pmin(pmax(model_data$response_raw/100, 1e-06), 1 - 1e-06)
    model_data$y <- qlogis(proportion)
    inverse_function <- function(x) 100 * plogis(x)
  }
  else {
    model_data$y <- model_data$response_raw
    inverse_function <- identity
  }
  warning_messages <- character()
  fit <- withCallingHandlers(mgcv::gamm(y ~ s(time_num, k = k, bs = "cr"), data = model_data, random = list(cluster = ~1),
    correlation = nlme::corCAR1(form = ~time_num | cluster), method = "REML"), warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  time_grid <- seq(min(model_data$time_num), max(model_data$time_num), length.out = 240L)
  prediction <- predict(fit$gam, newdata = data.frame(time_num = time_grid), se.fit = TRUE)
  prediction_data <- data.frame(Date = min(model_data$Date) + time_grid, fit = inverse_function(as.numeric(prediction$fit)),
    lower = inverse_function(as.numeric(prediction$fit - 1.96 * prediction$se.fit)), upper = inverse_function(as.numeric(prediction$fit +
      1.96 * prediction$se.fit)), stringsAsFactors = FALSE)
  if (transformation == "identity") {
    prediction_data$lower <- pmax(prediction_data$lower, 0)
  }
  smooth_table <- summary(fit$gam)$s.table
  p_column <- grep("^p", colnames(smooth_table), ignore.case = TRUE, value = TRUE)[1]
  statistic_column <- if ("F" %in% colnames(smooth_table)) {
    "F"
  }
  else {
    colnames(smooth_table)[ncol(smooth_table) - 1L]
  }
  fitted_transformed <- as.numeric(predict(fit$gam))
  pseudo_r2 <- 1 - sum((model_data$y - fitted_transformed)^2)/sum((model_data$y - mean(model_data$y))^2)
  car1_phi <- tryCatch(as.numeric(coef(fit$lme$modelStruct$corStruct, unconstrained = FALSE)[1]), error = function(e) NA_real_)
  statistics <- data.frame(response = response, n_samples = nrow(model_data), n_months = length(unique(model_data$Date)),
    n_clusters = nlevels(model_data$cluster), transformation = transformation, method = "GAMM_REML_cubic_regression_spline_CAR1",
    k = k, edf = unname(smooth_table[1, "edf"]), statistic = unname(smooth_table[1, statistic_column]),
    p_smooth = unname(smooth_table[1, p_column]), pseudo_r2 = pseudo_r2, car1_phi = car1_phi, model_warnings = paste(unique(warning_messages),
      collapse = " | "), stringsAsFactors = FALSE)
  list(predictions = prediction_data, statistics = statistics, model = fit)
}

richness_fit <- fit_overall_gamm(sample_metrics, response = "Pathogen_richness", transformation = "identity")

relative_fit <- fit_overall_gamm(sample_metrics, response = "Pathogen_relative_abundance", transformation = "logit_percent")

richness_predictions <- richness_fit$predictions

relative_predictions <- relative_fit$predictions

model_statistics <- bind_rows(richness_fit$statistics, relative_fit$statistics)

richness_axis_max <- ceiling(max(sample_metrics$Pathogen_richness, richness_predictions$upper)/25) *
  25

relative_axis_max <- ceiling(max(sample_metrics$Pathogen_relative_abundance, relative_predictions$upper)/5) *
  5

axis_multiplier <- richness_axis_max/relative_axis_max

sample_metrics <- sample_metrics %>%
  mutate(Relative_abundance_scaled = Pathogen_relative_abundance * axis_multiplier)

relative_predictions <- relative_predictions %>%
  mutate(fit_scaled = fit * axis_multiplier, lower_scaled = lower * axis_multiplier, upper_scaled = upper *
    axis_multiplier)

monthly_summary <- sample_metrics %>%
  group_by(Date, month) %>%
  summarise(n_samples = n(), mean_richness = mean(Pathogen_richness), sd_richness = sd(Pathogen_richness),
    mean_relative_abundance = mean(Pathogen_relative_abundance), sd_relative_abundance = sd(Pathogen_relative_abundance),
    .groups = "drop")

write.csv(sample_metrics, file.path(output_dir, "Pathogen_richness_relative_abundance_sample_data.csv"),
  row.names = FALSE)

write.csv(monthly_summary, file.path(output_dir, "Pathogen_richness_relative_abundance_monthly_summary.csv"),
  row.names = FALSE)

write.csv(model_statistics, file.path(output_dir, "Pathogen_richness_relative_abundance_GAMM_statistics.csv"),
  row.names = FALSE)

write.csv(transform(richness_predictions, response = "Pathogen_richness", scale_for_plot = 1), file.path(output_dir,
  "Pathogen_richness_GAMM_predictions.csv"), row.names = FALSE)

write.csv(transform(relative_predictions, response = "Pathogen_relative_abundance", scale_for_plot = axis_multiplier),
  file.path(output_dir, "Pathogen_relative_abundance_GAMM_predictions.csv"), row.names = FALSE)

water_observed <- read.csv(file.path(input_dir, "WaterLevel_Monthly_Data.csv"), check.names = FALSE,
  stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date))

water_predictions <- read.csv(file.path(input_dir, "GAMM_WaterLevel_Predictions.csv"), check.names = FALSE,
  stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date))

water_statistics <- read.csv(file.path(input_dir, "GAMM_WaterLevel_Model_Statistics.csv"), check.names = FALSE,
  stringsAsFactors = FALSE)

if (nrow(water_statistics) != 1L) {
  stop("Water-level statistics must contain exactly one row.")
}

month_breaks <- sort(unique(c(sample_metrics$Date, water_observed$Date)))

main_month_breaks <- month_breaks[seq(1L, length(month_breaks), by = 2L)]

inset_month_breaks <- month_breaks[unique(round(seq(1, length(month_breaks), length.out = 5)))]

date_limits <- range(month_breaks, na.rm = TRUE)

date_limits_plot <- date_limits + c(-7, 7)

theme_cns_closed <- function() {
  theme_classic(base_family = "sans", base_size = 12.5) + theme(panel.background = element_rect(fill = "white",
    colour = NA), panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8), axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.45), axis.ticks.length = grid::unit(1.6,
      "mm"), axis.text = element_text(colour = "black", size = 8.2), axis.title = element_text(colour = "black",
      size = 9.5, face = "bold"), plot.title = element_text(colour = "black", size = 10.5, face = "bold",
      hjust = 0), plot.subtitle = element_text(colour = "#333333", size = 7.7, lineheight = 1.08,
      margin = margin(b = 5)), legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 8.2),
    legend.key.width = grid::unit(10, "mm"), legend.key.height = grid::unit(3.5, "mm"), plot.margin = margin(7,
      10, 5, 7))
}

metric_colors <- c(`Pathogen richness` = "#6A3D9A", `Pathogen relative abundance` = "#B8860B")

subtitle_text <- paste0("Overall GAMM: richness P", format_p(richness_fit$statistics$p_smooth), "; relative abundance P",
  format_p(relative_fit$statistics$p_smooth), " (logit scale); n=", nrow(sample_metrics))

p_main <- ggplot() + geom_ribbon(data = richness_predictions, aes(x = Date, ymin = lower, ymax = upper),
  fill = unname(metric_colors["Pathogen richness"]), alpha = 0.14, colour = NA) + geom_ribbon(data = relative_predictions,
  aes(x = Date, ymin = lower_scaled, ymax = upper_scaled), fill = unname(metric_colors["Pathogen relative abundance"]),
  alpha = 0.14, colour = NA) + geom_line(data = richness_predictions, aes(x = Date, y = fit, colour = "Pathogen richness",
  linetype = "Pathogen richness"), linewidth = 1.02, lineend = "round") + geom_line(data = relative_predictions,
  aes(x = Date, y = fit_scaled, colour = "Pathogen relative abundance", linetype = "Pathogen relative abundance"),
  linewidth = 1.02, lineend = "round") + geom_point(data = sample_metrics, aes(x = Date, y = Pathogen_richness,
  colour = "Pathogen richness"), shape = 16, size = 1.55, alpha = 0.46, position = position_jitter(width = 4.2,
  height = 0, seed = 20260811)) + geom_point(data = sample_metrics, aes(x = Date, y = Relative_abundance_scaled,
  colour = "Pathogen relative abundance"), shape = 18, size = 1.65, alpha = 0.52, position = position_jitter(width = 4.2,
  height = 0, seed = 20260812)) + scale_colour_manual(values = metric_colors, breaks = names(metric_colors),
  name = NULL) + scale_linetype_manual(values = c(`Pathogen richness` = "solid", `Pathogen relative abundance` = "22"),
  breaks = names(metric_colors), name = NULL) + scale_x_date(limits = date_limits_plot, breaks = main_month_breaks,
  date_labels = "%Y-%m", expand = expansion(mult = c(0.025, 0.025))) + scale_y_continuous(limits = c(0,
  richness_axis_max), expand = expansion(mult = c(0, 0.025)), name = "Pathogen richness", sec.axis = sec_axis(~./axis_multiplier,
  name = "Total pathogen relative abundance (%)", breaks = scales::pretty_breaks(n = 6))) + labs(title = "Temporal dynamics of potential pathogens",
  subtitle = subtitle_text, x = "Month") + guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = c(1.05,
  1.05), shape = c(16, 18), alpha = 1)), linetype = "none") + theme_cns_closed() + theme(axis.text.x = element_text(angle = 45,
  hjust = 1, vjust = 1, size = 11.2), axis.title.y.left = element_text(colour = "black"), axis.text.y.left = element_text(colour = "black"),
  axis.ticks.y.left = element_line(colour = "black"), axis.title.y.right = element_text(colour = "black",
    margin = margin(l = 6)), axis.text.y.right = element_text(colour = "black"), axis.ticks.y.right = element_line(colour = "black"),
  legend.position = "top", legend.direction = "horizontal", legend.justification = "left", legend.margin = margin(0,
    0, 1, 0))

p_water_inset <- ggplot(water_observed, aes(x = Date, y = WL)) + geom_ribbon(data = water_predictions,
  aes(x = Date, ymin = lower, ymax = upper), inherit.aes = FALSE, fill = "#BDBDBD", alpha = 0.45, colour = NA) +
  geom_line(data = water_predictions, aes(x = Date, y = fit), inherit.aes = FALSE, colour = "#4D4D4D",
    linewidth = 0.95) + geom_point(colour = "#4D4D4D", fill = "white", shape = 21, stroke = 0.45,
  size = 2.15, alpha = 0.95) + scale_x_date(limits = date_limits_plot, breaks = inset_month_breaks,
  date_labels = "%Y-%m", expand = expansion(mult = c(0.025, 0.025))) + labs(title = "Water level",
  subtitle = paste0("GAMM P", format_p(water_statistics$smooth_p), "; 95% CI"), x = NULL, y = "Water level (m)") +
  theme_cns_closed() + theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7),
  axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 45, hjust = 1,
    vjust = 1, size = 7), axis.title = element_text(colour = "black", size = 6.2, face = "bold"),
  axis.ticks.length = grid::unit(0.8, "mm"), plot.title = element_text(colour = "black", size = 7.2,
    face = "bold", hjust = 0, margin = margin(l = 12)), plot.subtitle = element_text(colour = "#333333",
    size = 5.4, lineheight = 1, margin = margin(b = 2)), plot.margin = margin(4, 4, 4, 4))

combined_plot <- p_main + inset_element(p_water_inset, left = 0.045, bottom = 0.52, right = 0.49, top = 0.965,
  align_to = "panel", on_top = TRUE, clip = TRUE) + plot_annotation(caption = paste0("Lines show overall GAMM fits; shaded ribbons indicate 95% confidence intervals.\n",
  "Relative abundance is the summed pathogen fraction of the complete rarefied community."), theme = theme(plot.caption = element_text(family = "sans",
  size = 7.2, colour = "#333333", hjust = 0, lineheight = 1.05, margin = margin(t = 4)), plot.background = element_rect(fill = "white",
  colour = NA)))

prefix <- file.path(output_dir, "Pathogen_richness_relative_abundance_with_WL_inset")

figure_width <- 7.6

figure_height <- 5.9

ggsave(paste0(prefix, ".png"), combined_plot, width = figure_width, height = figure_height, dpi = 600,
  bg = "white", device = ragg::agg_png)

ggsave(paste0(prefix, ".svg"), combined_plot, width = figure_width, height = figure_height, bg = "white",
  device = svglite::svglite)

ggsave(paste0(prefix, ".pdf"), combined_plot, width = figure_width, height = figure_height, bg = "white",
  device = cairo_pdf)

ggsave(paste0(prefix, ".tiff"), combined_plot, width = figure_width, height = figure_height, dpi = 600,
  bg = "white", compression = "lzw", device = ragg::agg_tiff)

capture.output({
  cat("Pathogen dual-axis temporal trend audit\n")
  cat("=======================================\n")
  cat("Samples: ", nrow(sample_metrics), "\n", sep = "")
  cat("Months: ", length(unique(sample_metrics$Date)), "\n", sep = "")
  cat("Repeated-measure clusters: ", nlevels(sample_metrics$cluster), "\n", sep = "")
  cat("Complete-community library sizes: ", min(sample_metrics$Community_library_size), " to ", max(sample_metrics$Community_library_size),
    "\n", sep = "")
  cat("Pathogen richness range: ", min(sample_metrics$Pathogen_richness), " to ", max(sample_metrics$Pathogen_richness),
    "\n", sep = "")
  cat("Pathogen relative abundance range (%): ", sprintf("%.3f", min(sample_metrics$Pathogen_relative_abundance)),
    " to ", sprintf("%.3f", max(sample_metrics$Pathogen_relative_abundance)), "\n", sep = "")
  cat("Secondary-axis multiplier: ", sprintf("%.8f", axis_multiplier), "\n", sep = "")
  cat("Richness smooth P: ", richness_fit$statistics$p_smooth, "\n", sep = "")
  cat("Relative-abundance smooth P: ", relative_fit$statistics$p_smooth, "\n", sep = "")
  cat("Relative-abundance model uses a logit transformation.\n")
}, file = file.path(output_dir, "Pathogen_richness_relative_abundance_analysis_log.txt"))

message("Completed: ", normalizePath(output_dir, winslash = "/"))
