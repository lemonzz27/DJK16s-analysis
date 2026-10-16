options(stringsAsFactors = FALSE)

required_packages <- c("mgcv", "nlme")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

s <- mgcv::s

args <- commandArgs(trailingOnly = TRUE)

all_args <- commandArgs(trailingOnly = FALSE)

script_arg <- grep("^--file=", all_args, value = TRUE)

script_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/icamp_assembly_plots")

input_name <- "DJK16s_bin48_pairwise_process_importance.csv"

input_candidates <- c(file.path(dirname(script_dir), input_name), file.path(dirname(dirname(script_dir)),
  input_name), file.path(getwd(), input_name))

existing_input <- input_candidates[file.exists(input_candidates)]

default_input <- if (length(existing_input) > 0L) existing_input[1L] else input_candidates[1L]

input_file <- if (length(args) >= 1L) args[1L] else default_input

output_dir <- if (length(args) >= 2L) args[2L] else script_dir

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

input_file <- normalizePath(input_file, winslash = "/")

output_dir <- normalizePath(output_dir, winslash = "/")

process_cols <- c("Heterogeneous.Selection", "Homogeneous.Selection", "Dispersal.Limitation", "Homogenizing.Dispersal",
  "Drift.and.Others")

group_order <- c("ConF", "USF", "All")

group_colors <- c(ConF = "#D55E00", USF = "#009E73", All = "#0072B2")

parse_sample_metadata <- function(samples) {
  cleaned <- gsub("_", "", trimws(samples), fixed = TRUE)
  pattern <- "^([A-Za-z]+)([0-9]{2})([0-9]{2})(0[12])$"
  matches <- regexec(pattern, cleaned)
  pieces <- regmatches(cleaned, matches)
  valid <- lengths(pieces) == 5L
  if (!all(valid)) {
    bad <- unique(samples[!valid])
    stop("Could not parse sample name(s): ", paste(utils::head(bad, 8L), collapse = ", "), ". Expected <site><YY><MM><01|02>.")
  }
  values <- do.call(rbind, lapply(pieces, function(x) x[2:5]))
  month_number <- as.integer(values[, 3L])
  if (any(month_number < 1L | month_number > 12L)) {
    stop("Invalid month code found in sample names.")
  }
  size_code <- values[, 4L]
  size_group <- ifelse(size_code == "01", "USF", "ConF")
  data.frame(Sample = samples, Site = values[, 1L], Date = as.Date(paste0("20", values[, 2L], "-",
    values[, 3L], "-01")), Size = size_group, Node = paste(values[, 1L], size_group, sep = "_"),
    check.names = FALSE)
}

canonical_pair <- function(first, second) {
  ifelse(first <= second, paste(first, second, sep = "--"), paste(second, first, sep = "--"))
}

format_p <- function(value) {
  if (length(value) == 0L || is.na(value) || !is.finite(value))
    return("NA")
  if (value < 0.001)
    return("<0.001")
  sprintf("=%.3f", value)
}

extract_smooth <- function(model, row_index = 1L) {
  smooth_table <- summary(model$gam)$s.table
  if (is.null(smooth_table) || nrow(smooth_table) < row_index) {
    return(c(edf = NA_real_, statistic = NA_real_, p_value = NA_real_))
  }
  statistic_column <- if ("F" %in% colnames(smooth_table)) {
    "F"
  }
  else if ("Chi.sq" %in% colnames(smooth_table)) {
    "Chi.sq"
  }
  else {
    colnames(smooth_table)[ncol(smooth_table) - 1L]
  }
  p_column <- grep("^p", colnames(smooth_table), ignore.case = TRUE, value = TRUE)[1L]
  c(edf = unname(smooth_table[row_index, "edf"]), statistic = unname(smooth_table[row_index, statistic_column]),
    p_value = unname(smooth_table[row_index, p_column]))
}

extract_car1_phi <- function(model) {
  tryCatch(as.numeric(coef(model$lme$modelStruct$corStruct, unconstrained = FALSE)[1L]), error = function(e) NA_real_)
}

run_gamm_captured <- function(formula, data) {
  warnings <- character()
  model <- withCallingHandlers(mgcv::gamm(formula = formula, data = data, random = list(cluster = ~1),
    correlation = nlme::corCAR1(form = ~time_num | cluster), method = "REML"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(model = model, warnings = paste(unique(warnings), collapse = " | "))
}

predict_group_curve <- function(model, time_grid, group_name, ordered_levels) {
  new_data <- data.frame(time_num = time_grid, group_ord = ordered(rep(group_name, length(time_grid)),
    levels = ordered_levels))
  prediction <- predict(model$gam, newdata = new_data, se.fit = TRUE)
  data.frame(time_num = time_grid, Group = group_name, fit = exp(as.numeric(prediction$fit)), lower = exp(as.numeric(prediction$fit -
    1.96 * prediction$se.fit)), upper = exp(as.numeric(prediction$fit + 1.96 * prediction$se.fit)),
    stringsAsFactors = FALSE)
}

monthly_summary <- function(data) {
  split_data <- split(data, interaction(data$Group, data$Date, drop = TRUE))
  rows <- lapply(split_data, function(part) {
    data.frame(Date = part$Date[1L], Month = format(part$Date[1L], "%Y-%m"), Group = part$Group[1L],
      n_pairs = nrow(part), mean_ratio = mean(part$Ratio), sd_ratio = stats::sd(part$Ratio), sem_ratio = stats::sd(part$Ratio)/sqrt(nrow(part)),
      median_ratio = stats::median(part$Ratio), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$GroupOrder <- match(out$Group, group_order)
  out <- out[order(out$GroupOrder, out$Date), ]
  out$GroupOrder <- NULL
  rownames(out) <- NULL
  out
}

message("Reading: ", input_file)

pairwise <- utils::read.csv(input_file, check.names = FALSE)

required_cols <- c("sample1", "sample2", process_cols)

missing_cols <- setdiff(required_cols, names(pairwise))

if (length(missing_cols) > 0L) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

if (anyNA(pairwise[required_cols])) {
  stop("Missing values were found in required columns.")
}

for (column in process_cols) {
  pairwise[[column]] <- as.numeric(pairwise[[column]])
}

process_matrix <- as.matrix(pairwise[process_cols])

if (any(!is.finite(process_matrix))) {
  stop("Non-finite process importance values were found.")
}

row_total_error <- max(abs(rowSums(process_matrix) - 1))

if (row_total_error > 1e-08) {
  stop("Process importance rows do not sum to 1; max error = ", row_total_error)
}

deterministic <- pairwise[["Heterogeneous.Selection"]] + pairwise[["Homogeneous.Selection"]]

stochastic <- pairwise[["Dispersal.Limitation"]] + pairwise[["Homogenizing.Dispersal"]] + pairwise[["Drift.and.Others"]]

if (any(deterministic <= 0)) {
  stop("A non-positive deterministic-process value prevents ratio calculation.")
}

ratio <- stochastic/deterministic

meta1 <- parse_sample_metadata(pairwise$sample1)

meta2 <- parse_sample_metadata(pairwise$sample2)

same_month <- meta1$Date == meta2$Date

pair_cluster <- canonical_pair(meta1$Node, meta2$Node)

make_group_data <- function(group_name, keep) {
  data.frame(Date = meta1$Date[keep], Group = group_name, Ratio = ratio[keep], Deterministic = deterministic[keep],
    Stochastic = stochastic[keep], cluster = factor(pair_cluster[keep]), sample1 = pairwise$sample1[keep],
    sample2 = pairwise$sample2[keep], stringsAsFactors = FALSE)
}

all_data <- make_group_data("All", same_month)

conf_data <- make_group_data("ConF", same_month & meta1$Size == "ConF" & meta2$Size == "ConF")

usf_data <- make_group_data("USF", same_month & meta1$Size == "USF" & meta2$Size == "USF")

group_data <- rbind(conf_data, usf_data)

plot_data <- rbind(conf_data, usf_data, all_data)

plot_data$Group <- factor(plot_data$Group, levels = group_order)

duplicate_group_key <- paste(group_data$cluster, group_data$Date, sep = "|")

duplicate_all_key <- paste(all_data$cluster, all_data$Date, sep = "|")

if (anyDuplicated(duplicate_group_key) || anyDuplicated(duplicate_all_key)) {
  stop("Duplicate observations occur within a repeated pair-by-month sequence.")
}

minimum_date <- min(plot_data$Date)

maximum_date <- max(plot_data$Date)

group_data$time_num <- as.numeric(group_data$Date - minimum_date)

group_data$log_ratio <- log(group_data$Ratio)

group_data$cluster <- factor(group_data$cluster)

group_data$group_ord <- ordered(group_data$Group, levels = c("ConF", "USF"))

group_formula <- log_ratio ~ group_ord + s(time_num, k = 6L, bs = "cr") + s(time_num, by = group_ord,
  k = 6L, bs = "cr")

primary_run <- run_gamm_captured(group_formula, group_data)

reverse_group_data <- group_data

reverse_group_data$group_ord <- ordered(reverse_group_data$Group, levels = c("USF", "ConF"))

reverse_run <- run_gamm_captured(group_formula, reverse_group_data)

all_model_data <- all_data

all_model_data$time_num <- as.numeric(all_model_data$Date - minimum_date)

all_model_data$log_ratio <- log(all_model_data$Ratio)

all_model_data$cluster <- factor(all_model_data$cluster)

all_run <- run_gamm_captured(log_ratio ~ s(time_num, k = 6L, bs = "cr"), all_model_data)

time_grid <- seq(as.numeric(minimum_date - minimum_date), as.numeric(maximum_date - minimum_date), length.out = 300L)

predictions <- rbind(predict_group_curve(primary_run$model, time_grid, "ConF", c("ConF", "USF")), predict_group_curve(primary_run$model,
  time_grid, "USF", c("ConF", "USF")))

all_prediction <- predict(all_run$model$gam, newdata = data.frame(time_num = time_grid), se.fit = TRUE)

predictions <- rbind(predictions, data.frame(time_num = time_grid, Group = "All", fit = exp(as.numeric(all_prediction$fit)),
  lower = exp(as.numeric(all_prediction$fit - 1.96 * all_prediction$se.fit)), upper = exp(as.numeric(all_prediction$fit +
    1.96 * all_prediction$se.fit)), stringsAsFactors = FALSE))

predictions$Date <- minimum_date + predictions$time_num

predictions$Group <- factor(predictions$Group, levels = group_order)

conf_smooth <- extract_smooth(primary_run$model, 1L)

difference_smooth <- extract_smooth(primary_run$model, 2L)

usf_smooth <- extract_smooth(reverse_run$model, 1L)

all_smooth <- extract_smooth(all_run$model, 1L)

group_fitted <- as.numeric(predict(primary_run$model$gam))

all_fitted <- as.numeric(predict(all_run$model$gam))

group_pseudo_r2 <- 1 - sum((group_data$log_ratio - group_fitted)^2)/sum((group_data$log_ratio - mean(group_data$log_ratio))^2)

all_pseudo_r2 <- 1 - sum((all_model_data$log_ratio - all_fitted)^2)/sum((all_model_data$log_ratio - mean(all_model_data$log_ratio))^2)

statistics <- data.frame(response = "log_stochastic_to_deterministic_ratio", n_ConF = nrow(conf_data),
  n_USF = nrow(usf_data), n_All = nrow(all_data), observed_months = length(unique(plot_data$Date)),
  k = 6L, method = "GAMM_REML_cubic_regression_spline_CAR1_log_ratio", p_ConF_smooth = conf_smooth[["p_value"]],
  edf_ConF_smooth = conf_smooth[["edf"]], statistic_ConF_smooth = conf_smooth[["statistic"]], p_USF_smooth = usf_smooth[["p_value"]],
  edf_USF_smooth = usf_smooth[["edf"]], statistic_USF_smooth = usf_smooth[["statistic"]], p_difference_smooth = difference_smooth[["p_value"]],
  edf_difference_smooth = difference_smooth[["edf"]], statistic_difference_smooth = difference_smooth[["statistic"]],
  p_All_smooth = all_smooth[["p_value"]], edf_All_smooth = all_smooth[["edf"]], statistic_All_smooth = all_smooth[["statistic"]],
  car1_phi_group_model = extract_car1_phi(primary_run$model), car1_phi_All_model = extract_car1_phi(all_run$model),
  pseudo_r2_group_model_log_scale = group_pseudo_r2, pseudo_r2_All_model_log_scale = all_pseudo_r2,
  model_warnings = paste(Filter(nzchar, c(primary_run$warnings, reverse_run$warnings, all_run$warnings)),
    collapse = " | "), stringsAsFactors = FALSE)

summary_data <- monthly_summary(plot_data)

display_y_limits <- c(0, 2)

statistics$display_y_min <- display_y_limits[1L]

statistics$display_y_max <- display_y_limits[2L]

statistics$n_plot_rows_above_display_limit <- sum(plot_data$Ratio > display_y_limits[2L])

utils::write.csv(plot_data, file.path(output_dir, "stochastic_deterministic_ratio_pairwise_data.csv"),
  row.names = FALSE)

utils::write.csv(summary_data, file.path(output_dir, "stochastic_deterministic_ratio_monthly_summary.csv"),
  row.names = FALSE)

utils::write.csv(predictions, file.path(output_dir, "stochastic_deterministic_ratio_GAMM_predictions.csv"),
  row.names = FALSE)

utils::write.csv(statistics, file.path(output_dir, "stochastic_deterministic_ratio_GAMM_statistics.csv"),
  row.names = FALSE)

subtitle_line1 <- paste0("GAMM smooth (log ratio): ConF P", format_p(statistics$p_ConF_smooth), "; USF P",
  format_p(statistics$p_USF_smooth), "; All P", format_p(statistics$p_All_smooth))

subtitle_line2 <- paste0("Difference smooth (USF-ConF): P", format_p(statistics$p_difference_smooth))

draw_ratio_trend <- function() {
  par(mar = c(5.9, 5.2, 1, 1.2), family = "sans", mgp = c(2.7, 0.75, 0), xaxs = "i", yaxs = "i")
  x_limits <- c(minimum_date - 7, maximum_date + 7)
  y_limits <- display_y_limits
  plot(x = x_limits, y = y_limits, type = "n", axes = FALSE, ann = FALSE, xlim = x_limits, ylim = y_limits)
  for (group_name in group_order) {
    curve <- predictions[predictions$Group == group_name, ]
    curve <- curve[order(curve$Date), ]
    polygon(x = c(curve$Date, rev(curve$Date)), y = c(curve$lower, rev(curve$upper)), col = grDevices::adjustcolor(group_colors[group_name],
      alpha.f = 0.18), border = NA)
  }
  abline(h = 1, col = "#666666", lty = 2, lwd = 1)
  set.seed(20260824)
  for (group_name in group_order) {
    observations <- plot_data[plot_data$Group == group_name, ]
    observations <- observations[observations$Ratio >= y_limits[1L] & observations$Ratio <= y_limits[2L],
      ]
    jittered_date <- observations$Date + stats::runif(nrow(observations), -5, 5)
    points(jittered_date, observations$Ratio, pch = 16, cex = 0.58, col = grDevices::adjustcolor(group_colors[group_name],
      alpha.f = 0.36))
  }
  for (group_name in group_order) {
    curve <- predictions[predictions$Group == group_name, ]
    curve <- curve[order(curve$Date), ]
    lines(curve$Date, curve$fit, col = group_colors[group_name], lwd = 2.6)
  }
  month_ticks_all <- sort(unique(plot_data$Date))
  month_ticks <- month_ticks_all[seq(1L, length(month_ticks_all), by = 2L)]
  axis(1, at = month_ticks, labels = FALSE, tcl = -0.22)
  user_coords <- par("usr")
  text(x = month_ticks, y = user_coords[3L] - 0.025 * (user_coords[4L] - user_coords[3L]), labels = format(month_ticks,
    "%Y-%m"), srt = 45, adj = 1, xpd = NA, cex = 1.05)
  y_ticks <- pretty(y_limits, n = 6)
  y_ticks <- y_ticks[y_ticks >= y_limits[1L] & y_ticks <= y_limits[2L]]
  axis(2, at = y_ticks, labels = format(y_ticks, trim = TRUE), las = 1, tcl = -0.22, cex.axis = 1.05)
  box(bty = "o", lwd = 0.85)
  mtext("Month", side = 1, line = 4.6, cex = 1.15, font = 2)
  mtext("Stochastic / deterministic ratio", side = 2, line = 3.3, las = 0, font = 2, cex = 1.15)
  text(x = user_coords[2L] - 0.01 * (user_coords[2L] - user_coords[1L]), y = 1 + 0.025 * (user_coords[4L] -
    user_coords[3L]), labels = "Ratio = 1", adj = 1, col = "#555555", cex = 0.9)
  legend("topright", inset = c(0.025, 0.025), yjust = 1, legend = group_order, col = group_colors[group_order],
    lwd = 2.6, pch = 16, pt.cex = 0.8, bty = "n", title = "Group", cex = 1.02, y.intersp = 1.15)
}

pdf_path <- file.path(output_dir, "iCAMP_stochastic_deterministic_ratio_trend.pdf")

png_path <- file.path(output_dir, "iCAMP_stochastic_deterministic_ratio_trend.png")

grDevices::pdf(pdf_path, width = 7.2, height = 5, family = "Helvetica", useDingbats = FALSE, bg = "white")

draw_ratio_trend()

grDevices::dev.off()

grDevices::png(png_path, width = 7.2, height = 5, units = "in", res = 600, bg = "white")

draw_ratio_trend()

grDevices::dev.off()

message("Completed successfully.")

message("  All pairwise observations: ", nrow(all_data))

message("  ConF pairwise observations: ", nrow(conf_data))

message("  USF pairwise observations: ", nrow(usf_data))

message("  Observed months: ", length(unique(plot_data$Date)), " (", format(minimum_date, "%Y-%m"), " to ",
  format(maximum_date, "%Y-%m"), ")")

message("  Ratio range: ", paste(signif(range(plot_data$Ratio), 4), collapse = " to "))

message("  Displayed y range: ", paste(display_y_limits, collapse = " to "))

message("  Plot rows above display limit (retained in GAMM): ", sum(plot_data$Ratio > display_y_limits[2L]))

message("  Output directory: ", output_dir)
