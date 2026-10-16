required_gamm_packages <- c("mgcv", "nlme")

missing_gamm_packages <- required_gamm_packages[!vapply(required_gamm_packages, requireNamespace, logical(1),
  quietly = TRUE)]

if (length(missing_gamm_packages) > 0L) {
  stop("缺少GAMM所需R包：", paste(missing_gamm_packages, collapse = ", "))
}

format_gamm_p <- function(p_value) {
  if (length(p_value) == 0L || is.na(p_value) || !is.finite(p_value)) {
    return("NA")
  }
  if (p_value < 0.001)
    return("<0.001")
  sprintf("=%.3f", p_value)
}

extract_gamm_smooth <- function(model, row_index = 1L) {
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
  p_column <- grep("^p", colnames(smooth_table), ignore.case = TRUE, value = TRUE)[1]
  c(edf = unname(smooth_table[row_index, "edf"]), statistic = unname(smooth_table[row_index, statistic_column]),
    p_value = unname(smooth_table[row_index, p_column]))
}

extract_car1_phi <- function(model) {
  tryCatch(as.numeric(coef(model$lme$modelStruct$corStruct, unconstrained = FALSE)[1]), error = function(e) NA_real_)
}

run_gamm_captured <- function(formula, data, use_random_intercept) {
  warning_messages <- character()
  arguments <- list(formula = formula, data = data, method = "REML", correlation = nlme::corCAR1(form = ~time_num |
    cluster))
  if (isTRUE(use_random_intercept)) {
    arguments$random <- list(cluster = ~1)
  }
  model <- withCallingHandlers(do.call(mgcv::gamm, arguments), warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(model = model, warnings = paste(unique(warning_messages), collapse = " | "))
}

predict_gamm_curve <- function(model, time_grid, group_name, ordered_levels) {
  new_data <- data.frame(time_num = time_grid, size_ord = ordered(rep(group_name, length(time_grid)),
    levels = ordered_levels))
  prediction <- predict(model$gam, newdata = new_data, se.fit = TRUE)
  data.frame(time_num = time_grid, size_group = group_name, fit = as.numeric(prediction$fit), lower = as.numeric(prediction$fit -
    1.96 * prediction$se.fit), upper = as.numeric(prediction$fit + 1.96 * prediction$se.fit), stringsAsFactors = FALSE)
}

fit_repeated_time_gamm <- function(data, response, date_col = "Date", group_col = "size_group", cluster_col,
  k = 6L, include_overall = TRUE, use_random_intercept = TRUE) {
  required_columns <- c(response, date_col, group_col, cluster_col)
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0L) {
    stop("GAMM输入缺少列：", paste(missing_columns, collapse = ", "))
  }
  model_data <- data.frame(y = as.numeric(data[[response]]), Date = as.Date(data[[date_col]]), size_group = as.character(data[[group_col]]),
    cluster = factor(data[[cluster_col]]), stringsAsFactors = FALSE)
  model_data <- model_data[is.finite(model_data$y) & !is.na(model_data$Date) & model_data$size_group %in%
    c("large", "small") & !is.na(model_data$cluster), , drop = FALSE]
  model_data$time_num <- as.numeric(model_data$Date - min(model_data$Date))
  model_data$size_group <- factor(model_data$size_group, levels = c("large", "small"))
  if (nrow(model_data) < 12L || length(unique(model_data$time_num)) < 6L) {
    stop(response, "的有效时间观测不足，无法拟合预设GAMM。")
  }
  duplicate_key <- paste(model_data$cluster, model_data$time_num, sep = "|")
  if (anyDuplicated(duplicate_key)) {
    stop(response, "在同一相关序列内存在重复时间点，无法使用corCAR1。")
  }
  model_data$size_ord <- ordered(as.character(model_data$size_group), levels = c("large", "small"))
  primary_formula <- y ~ size_ord + s(time_num, k = k, bs = "cr") + s(time_num, by = size_ord, k = k,
    bs = "cr")
  primary_run <- run_gamm_captured(primary_formula, model_data, use_random_intercept)
  small_reference_data <- model_data
  small_reference_data$size_ord <- ordered(as.character(small_reference_data$size_group), levels = c("small",
    "large"))
  small_reference_run <- run_gamm_captured(primary_formula, small_reference_data, use_random_intercept)
  large_smooth <- extract_gamm_smooth(primary_run$model, 1L)
  difference_smooth <- extract_gamm_smooth(primary_run$model, 2L)
  small_smooth <- extract_gamm_smooth(small_reference_run$model, 1L)
  time_grid <- seq(min(model_data$time_num), max(model_data$time_num), length.out = 200L)
  prediction_data <- rbind(predict_gamm_curve(primary_run$model, time_grid, "large", c("large", "small")),
    predict_gamm_curve(primary_run$model, time_grid, "small", c("large", "small")))
  overall_smooth <- c(edf = NA_real_, statistic = NA_real_, p_value = NA_real_)
  overall_rho <- NA_real_
  overall_pseudo_r2 <- NA_real_
  overall_warning <- ""
  if (isTRUE(include_overall)) {
    overall_formula <- y ~ s(time_num, k = k, bs = "cr")
    overall_run <- run_gamm_captured(overall_formula, model_data, use_random_intercept)
    overall_smooth <- extract_gamm_smooth(overall_run$model, 1L)
    overall_prediction <- predict(overall_run$model$gam, newdata = data.frame(time_num = time_grid),
      se.fit = TRUE)
    prediction_data <- rbind(prediction_data, data.frame(time_num = time_grid, size_group = "All",
      fit = as.numeric(overall_prediction$fit), lower = as.numeric(overall_prediction$fit - 1.96 *
        overall_prediction$se.fit), upper = as.numeric(overall_prediction$fit + 1.96 * overall_prediction$se.fit),
      stringsAsFactors = FALSE))
    fitted_overall <- as.numeric(predict(overall_run$model$gam))
    overall_pseudo_r2 <- 1 - sum((model_data$y - fitted_overall)^2)/sum((model_data$y - mean(model_data$y))^2)
    overall_rho <- extract_car1_phi(overall_run$model)
    overall_warning <- overall_run$warnings
  }
  prediction_data$Date <- min(model_data$Date) + prediction_data$time_num
  prediction_data$metric <- response
  prediction_data$size_group <- factor(prediction_data$size_group, levels = c("large", "small", "All"))
  fitted_primary <- as.numeric(predict(primary_run$model$gam))
  primary_pseudo_r2 <- 1 - sum((model_data$y - fitted_primary)^2)/sum((model_data$y - mean(model_data$y))^2)
  statistics <- data.frame(metric = response, n = nrow(model_data), observed_months = length(unique(model_data$time_num)),
    k = as.integer(k), method = "GAMM_REML_cubic_regression_spline_CAR1", p_large_smooth = large_smooth[["p_value"]],
    edf_large_smooth = large_smooth[["edf"]], statistic_large_smooth = large_smooth[["statistic"]],
    p_small_smooth = small_smooth[["p_value"]], edf_small_smooth = small_smooth[["edf"]], statistic_small_smooth = small_smooth[["statistic"]],
    p_difference_smooth = difference_smooth[["p_value"]], edf_difference_smooth = difference_smooth[["edf"]],
    statistic_difference_smooth = difference_smooth[["statistic"]], p_overall_smooth = overall_smooth[["p_value"]],
    edf_overall_smooth = overall_smooth[["edf"]], statistic_overall_smooth = overall_smooth[["statistic"]],
    car1_phi_group_model = extract_car1_phi(primary_run$model), car1_phi_overall_model = overall_rho,
    pseudo_r2_group_model = primary_pseudo_r2, pseudo_r2_overall_model = overall_pseudo_r2, model_warnings = paste(Filter(nzchar,
      c(primary_run$warnings, small_reference_run$warnings, overall_warning)), collapse = " | "),
    stringsAsFactors = FALSE)
  list(predictions = prediction_data, statistics = statistics, model_data = model_data)
}

add_gamm_bh_adjustment <- function(statistics, include_overall = TRUE) {
  p_columns <- c("p_large_smooth", "p_small_smooth", "p_difference_smooth")
  if (isTRUE(include_overall)) {
    p_columns <- c(p_columns, "p_overall_smooth")
  }
  for (column in p_columns) {
    statistics[[sub("^p_", "q_BH_", column)]] <- p.adjust(statistics[[column]], method = "BH")
  }
  statistics
}

gamm_subtitle <- function(statistics_row, include_overall = TRUE, note = NULL) {
  first_line <- paste0("GAMM smooth: large P", format_gamm_p(statistics_row$p_large_smooth), "; small P",
    format_gamm_p(statistics_row$p_small_smooth))
  if (isTRUE(include_overall)) {
    first_line <- paste0(first_line, "; All P", format_gamm_p(statistics_row$p_overall_smooth))
  }
  second_line <- paste0("Difference smooth (small-large): P", format_gamm_p(statistics_row$p_difference_smooth))
  if (!is.null(note) && nzchar(note)) {
    second_line <- paste0(second_line, "; ", note)
  }
  paste(first_line, second_line, sep = "\n")
}
