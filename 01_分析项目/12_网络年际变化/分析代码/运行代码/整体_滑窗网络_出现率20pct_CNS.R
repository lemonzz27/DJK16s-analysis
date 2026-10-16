rm(list = ls())

set.seed(123)

cat("\n>>> [INFO] 启动分析: ", as.character(Sys.time()), "\n")

packages <- c("psych", "igraph", "dplyr", "tidyr", "Hmisc", "ggplot2", "lubridate", "mgcv", "nlme")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    stop("缺少 R 包：", p, "。请先完成中央依赖准备。")
  }
}

source(normalizePath(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/12_网络年际变化/分析代码/运行代码/GAMM曲线工具.R"),
  winslash = "/", mustWork = TRUE), encoding = "UTF-8")

COR_R_CUTOFF <- 0.6

COR_P_CUTOFF <- 0.05

WINDOW_SIZE <- 3

STEP_SIZE <- 1

MIN_SAMPLES <- 5

MIN_OTU_SUM <- 0.1

MIN_PREVALENCE_PROP <- 0.2

OUTPUT_DIR <- Sys.getenv("SLIDING_NETWORK_OUTPUT_DIR", unset = "Results_prevalence20_CNS")

PLOT_DIR <- file.path(OUTPUT_DIR, "Plots")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

analysis_parameters <- data.frame(Parameter = c("WINDOW_SIZE_MONTHS", "STEP_SIZE_MONTHS", "MIN_PREVALENCE_PROP",
  "MIN_SAMPLES", "MIN_OTU_SUM", "COR_R_CUTOFF", "COR_P_CUTOFF", "P_ADJUST_METHOD", "RANDOM_SEED", "GAMM_K",
  "GAMM_METHOD", "GAMM_CORRELATION", "GAMM_CI_LEVEL"), Value = c(WINDOW_SIZE, STEP_SIZE, MIN_PREVALENCE_PROP,
  MIN_SAMPLES, MIN_OTU_SUM, COR_R_CUTOFF, COR_P_CUTOFF, "BH", 123, 5, "REML cubic regression spline",
  "continuous-time CAR(1) within size group", 0.95))

write.csv(analysis_parameters, file.path(OUTPUT_DIR, "Network_Analysis_Parameters.csv"), row.names = FALSE)

cat("\n[1/5] 正在读取数据并合并分组...\n")

tpm_data <- read.csv("tpm.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

group_month <- read.csv("groupmonth.csv", stringsAsFactors = FALSE, check.names = FALSE)

group_size <- read.csv("groupsize.csv", stringsAsFactors = FALSE, check.names = FALSE)

colnames(group_month) <- c("sample", "month_group")

colnames(group_size) <- c("sample", "size_group")

group_info <- group_month %>%
  inner_join(group_size, by = "sample")

common_samps <- intersect(group_info$sample, colnames(tpm_data))

tpm_data <- tpm_data[, common_samps, drop = FALSE]

group_info <- group_info[group_info$sample %in% common_samps, , drop = FALSE]

group_info$month_group <- as.numeric(group_info$month_group)

time_groups <- sort(unique(group_info$month_group))

cat(paste0("    识别到 ", length(time_groups), " 个时间分组: ", paste(time_groups, collapse = ", "),
  "\n"))

cat(paste0("    参数: 3 个月窗口，步长 1 个月，节点出现率阈值 ", MIN_PREVALENCE_PROP *
  100, "%\n"))

get_natural_connectivity <- function(g, weighted = FALSE) {
  if (vcount(g) == 0)
    return(NA)
  if (weighted && "weight" %in% edge_attr_names(g)) {
    adj <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
    adj <- abs(adj)
  }
  else {
    adj <- as_adjacency_matrix(g, sparse = FALSE)
    adj[adj != 0] <- 1
  }
  evals <- tryCatch(eigen(adj, only.values = TRUE)$values, error = function(e) NA)
  if (all(is.na(evals)))
    return(NA)
  nc <- log(mean(exp(evals)))
  return(nc)
}

get_global_efficiency <- function(g) {
  n <- vcount(g)
  if (n < 2)
    return(0)
  d_mat <- distances(g, weights = NA)
  inv_d <- 1/d_mat
  inv_d[is.infinite(inv_d)] <- 0
  diag(inv_d) <- NA
  eff <- mean(inv_d, na.rm = TRUE)
  if (is.nan(eff))
    eff <- 0
  return(eff)
}

get_vulnerability <- function(g) {
  nodes_num <- vcount(g)
  if (nodes_num < 3)
    return(0)
  E_global <- get_global_efficiency(g)
  if (E_global == 0)
    return(0)
  top_nodes <- names(sort(degree(g), decreasing = TRUE)[1:min(nodes_num, 3)])
  max_vul <- 0
  for (node in top_nodes) {
    g_rem <- delete_vertices(g, node)
    E_rem <- get_global_efficiency(g_rem)
    v <- (E_global - E_rem)/E_global
    if (!is.na(v) && v > max_vul) {
      max_vul <- v
    }
  }
  return(max_vul)
}

build_and_extract_metrics <- function(sub_tpm, window_name, size_group_name, date_val) {
  valid_otus <- rowSums(sub_tpm, na.rm = TRUE) > MIN_OTU_SUM
  sub_tpm <- sub_tpm[valid_otus, , drop = FALSE]
  min_prevalence <- ceiling(ncol(sub_tpm) * MIN_PREVALENCE_PROP)
  sub_tpm <- sub_tpm[rowSums(sub_tpm > 0, na.rm = TRUE) >= min_prevalence, , drop = FALSE]
  if (nrow(sub_tpm) < 10)
    return(NULL)
  res <- rcorr(as.matrix(t(sub_tpm)), type = "spearman")
  r_mat <- res$r
  p_mat <- res$P
  r_mat[is.na(r_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  p_vec <- p_mat[upper.tri(p_mat)]
  p_adj <- p.adjust(p_vec, method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat))
  p_adj_mat[upper.tri(p_adj_mat)] <- p_adj
  p_adj_mat <- p_adj_mat + t(p_adj_mat) - 1
  diag(p_adj_mat) <- 1
  mask <- (abs(r_mat) >= COR_R_CUTOFF) & (p_adj_mat < COR_P_CUTOFF)
  r_mat[!mask] <- 0
  diag(r_mat) <- 0
  if (sum(abs(r_mat) > 0) == 0)
    return(NULL)
  g <- graph_from_adjacency_matrix(r_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  num_nodes <- vcount(g)
  num_edges <- ecount(g)
  if (num_nodes < 3 || num_edges == 0)
    return(NULL)
  weights <- E(g)$weight
  density_val <- edge_density(g)
  avg_path_len <- mean_distance(g, directed = FALSE, weights = NA)
  diam_val <- diameter(g, directed = FALSE, weights = NA)
  trans_global <- transitivity(g, type = "global")
  trans_local <- mean(transitivity(g, type = "local"), na.rm = TRUE)
  avg_degree <- mean(degree(g), na.rm = TRUE)
  avg_betweenness <- mean(betweenness(g, normalized = TRUE, weights = NA), na.rm = TRUE)
  avg_edge_betweenness <- mean(edge_betweenness(g, directed = FALSE, weights = NA), na.rm = TRUE)
  avg_closeness <- mean(closeness(g, normalized = TRUE, weights = NA), na.rm = TRUE)
  eig <- tryCatch(mean(eigen_centrality(g, weights = NA)$vector, na.rm = TRUE), error = function(e) NA)
  info_centrality_sum <- sum(closeness(g, normalized = TRUE, weights = NA), na.rm = TRUE)
  mod_obj <- tryCatch(cluster_fast_greedy(as_undirected(g), weights = abs(weights)), error = function(e) NULL)
  modularity_val <- if (!is.null(mod_obj)) {
    modularity(mod_obj)
  }
  else {
    NA
  }
  cores <- coreness(g)
  k_core_mean <- mean(cores)
  efficiency_val <- get_global_efficiency(g)
  cohesion_pos <- sum(weights[weights > 0], na.rm = TRUE)
  cohesion_neg <- abs(sum(weights[weights < 0], na.rm = TRUE))
  prop_negative <- ifelse(num_edges > 0, sum(weights < 0)/num_edges, 0)
  robust_unweight <- get_natural_connectivity(g, weighted = FALSE)
  robust_weight <- get_natural_connectivity(g, weighted = TRUE)
  vulnerability_val <- get_vulnerability(g)
  stability_val <- 1 - vulnerability_val
  return(data.frame(Window = window_name, Date = date_val, size_group = size_group_name, Node = num_nodes,
    Edge = num_edges, Density = density_val, Distance = avg_path_len, Diameter = diam_val, Transitivity_global = trans_global,
    Transitivity_local = trans_local, Degree = avg_degree, Betweenness = avg_betweenness, Betweenness_edge = avg_edge_betweenness,
    Closeness = avg_closeness, Eigen_centrality = eig, Network_info_centrality = info_centrality_sum,
    Modularity = modularity_val, K_core_mean = k_core_mean, Network_efficiency = efficiency_val,
    Cohesion_Positive = cohesion_pos, Cohesion_Negative = cohesion_neg, Proportion_Negative = prop_negative,
    Robustness_unweight = robust_unweight, Robustness_weight = robust_weight, Vulnerability = vulnerability_val,
    Stability = stability_val))
}

cat(paste0("\n[2/5] 开始执行滑动窗口计算，窗口 = ", WINDOW_SIZE, " 个月...\n"))

total_windows <- length(time_groups) - WINDOW_SIZE + 1

if (total_windows < 1) {
  stop("时间分组数量少于窗口大小，无法进行滑动窗口分析。")
}

results_list <- list()

for (i in seq(1, total_windows, by = STEP_SIZE)) {
  cur_groups <- time_groups[i:(i + WINDOW_SIZE - 1)]
  win_label <- paste(cur_groups[1], cur_groups[length(cur_groups)], sep = "-")
  mid_month <- cur_groups[ceiling(WINDOW_SIZE/2)]
  win_date <- as.Date(paste0("20", mid_month, "01"), format = "%Y%m%d")
  cat(paste0("    -> 正在处理窗口 [", win_label, "]\n"))
  win_info <- group_info %>%
    filter(month_group %in% cur_groups)
  samps_large <- win_info %>%
    filter(size_group == "large") %>%
    pull(sample)
  if (length(samps_large) >= MIN_SAMPLES) {
    res_large <- build_and_extract_metrics(sub_tpm = tpm_data[, samps_large, drop = FALSE], window_name = win_label,
      size_group_name = "large", date_val = win_date)
    if (!is.null(res_large)) {
      results_list[[paste0(win_label, "_large")]] <- res_large
    }
  }
  samps_small <- win_info %>%
    filter(size_group == "small") %>%
    pull(sample)
  if (length(samps_small) >= MIN_SAMPLES) {
    res_small <- build_and_extract_metrics(sub_tpm = tpm_data[, samps_small, drop = FALSE], window_name = win_label,
      size_group_name = "small", date_val = win_date)
    if (!is.null(res_small)) {
      results_list[[paste0(win_label, "_small")]] <- res_small
    }
  }
}

cat("\n[3/5] 正在整合网络指标表格...\n")

final_results <- bind_rows(results_list)

if (nrow(final_results) == 0) {
  stop(paste0("没有成功构建任何网络。\n", "建议检查：\n", "1. MIN_SAMPLES 是否太高；\n",
    "2. COR_R_CUTOFF 是否太高；\n", "3. COR_P_CUTOFF 是否太严格；\n", "4. 每个窗口内 ConF/USF 样本数是否足够；\n",
    "5. tpm.csv、groupmonth.csv、groupsize.csv 的样本名是否一致。"))
}

out_file <- file.path(OUTPUT_DIR, "Sliding_Window_Topology_Metrics_ConF_USF.csv")

final_results_export <- final_results

final_results_export$size_group <- dplyr::recode(as.character(final_results_export$size_group), large = "ConF",
  small = "USF")

write.csv(final_results_export, out_file, row.names = FALSE)

cat(paste0("    指标表已保存至: ", out_file, "\n"))

cat("\n[4/5] 正在生成 ConF 和 USF 两组拓扑参数趋势图...\n")

final_results$size_group <- factor(final_results$size_group, levels = c("large", "small"))

metrics_to_plot <- colnames(final_results)[4:ncol(final_results)]

gamm_fit_list <- lapply(metrics_to_plot, function(var) {
  plot_data <- final_results %>%
    filter(!is.na(.data[[var]]), !is.infinite(.data[[var]]))
  fit_repeated_time_gamm(data = plot_data, response = var, cluster_col = "size_group", k = 5L, include_overall = FALSE,
    use_random_intercept = FALSE)
})

names(gamm_fit_list) <- metrics_to_plot

gamm_statistics <- bind_rows(lapply(gamm_fit_list, `[[`, "statistics"))

gamm_statistics <- add_gamm_bh_adjustment(gamm_statistics, include_overall = FALSE)

gamm_predictions <- bind_rows(lapply(gamm_fit_list, `[[`, "predictions"))

gamm_statistics_export <- gamm_statistics

names(gamm_statistics_export) <- stringr::str_replace_all(names(gamm_statistics_export), c(large = "ConF",
  small = "USF", difference = "difference_USF_minus_ConF"))

gamm_predictions_export <- gamm_predictions

gamm_predictions_export$size_group <- dplyr::recode(as.character(gamm_predictions_export$size_group),
  large = "ConF", small = "USF")

write.csv(gamm_statistics_export, file.path(OUTPUT_DIR, "GAMM_Network_Model_Statistics.csv"), row.names = FALSE)

write.csv(gamm_predictions_export, file.path(OUTPUT_DIR, "GAMM_Network_Predictions.csv"), row.names = FALSE)

theme_cns_closed <- function() {
  theme_classic(base_family = "sans", base_size = 12.5) + theme(panel.background = element_rect(fill = "white",
    colour = NA), panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8), axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.45), axis.ticks.length = grid::unit(1.6,
      "mm"), axis.text = element_text(colour = "black", size = 8.2), axis.title = element_text(colour = "black",
      size = 9.5, face = "bold"), plot.title = element_text(colour = "black", size = 10.5, face = "bold",
      hjust = 0), plot.subtitle = element_text(colour = "#333333", size = 7.7, lineheight = 1.08,
      margin = margin(b = 5)), plot.caption = element_text(colour = "#333333", size = 7.5, hjust = 0,
      margin = margin(t = 4)), legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 8.2),
    legend.key.width = grid::unit(10, "mm"), legend.key.height = grid::unit(3.5, "mm"), plot.margin = margin(7,
      8, 5, 7))
}

metric_labels <- c(Node = "Nodes", Edge = "Edges", Density = "Density", Distance = "Average path length",
  Diameter = "Diameter", Transitivity_global = "Global transitivity", Transitivity_local = "Mean local transitivity",
  Degree = "Mean degree", Betweenness = "Mean node betweenness", Betweenness_edge = "Mean edge betweenness",
  Closeness = "Mean closeness", Eigen_centrality = "Mean eigenvector centrality", Network_info_centrality = "Network information centrality",
  Modularity = "Modularity", K_core_mean = "Mean K-core", Network_efficiency = "Global efficiency",
  Cohesion_Positive = "Positive cohesion", Cohesion_Negative = "Negative cohesion", Proportion_Negative = "Proportion of negative edges",
  Robustness_unweight = "Natural connectivity (unweighted)", Robustness_weight = "Natural connectivity (weighted)",
  Vulnerability = "Vulnerability", Stability = "Stability")

group_colors <- c(large = "#D55E00", small = "#009E73")

group_labels <- c(large = "ConF", small = "USF")

month_breaks <- sort(unique(final_results$Date))

month_breaks_display <- month_breaks[seq(1L, length(month_breaks), by = 2L)]

date_limits_plot <- range(month_breaks, na.rm = TRUE) + c(-7, 7)

network_gamm_subtitle <- function(statistics_row) {
  paste0("GAMM smooth: ConF P", format_gamm_p(statistics_row$p_large_smooth), "; USF P", format_gamm_p(statistics_row$p_small_smooth),
    "\nDifference smooth (USF-ConF): P", format_gamm_p(statistics_row$p_difference_smooth))
}

for (var in metrics_to_plot) {
  plot_data <- final_results %>%
    filter(!is.na(.data[[var]]), !is.infinite(.data[[var]]))
  curve_data <- gamm_predictions %>%
    filter(metric == var)
  statistics_row <- gamm_statistics %>%
    filter(metric == var)
  p <- ggplot(plot_data, aes(x = Date, y = .data[[var]], color = size_group, fill = size_group)) +
    geom_ribbon(data = curve_data, aes(x = Date, ymin = lower, ymax = upper, fill = size_group),
      inherit.aes = FALSE, alpha = 0.18, color = NA) + geom_line(data = curve_data, aes(x = Date,
    y = fit, color = size_group), inherit.aes = FALSE, linewidth = 0.95) + geom_point(size = 1.65,
    alpha = 0.48, position = position_jitter(width = 5, height = 0, seed = 20260831)) + scale_color_manual(values = group_colors,
    breaks = names(group_colors), labels = group_labels, name = "Group") + scale_fill_manual(values = group_colors,
    breaks = names(group_colors), labels = group_labels, name = "Group") + scale_x_date(limits = date_limits_plot,
    breaks = month_breaks_display, date_labels = "%Y-%m", expand = expansion(mult = c(0.025, 0.025))) +
    labs(title = paste0("Network topology (", unname(metric_labels[[var]]), ")"), subtitle = network_gamm_subtitle(statistics_row),
      x = "Month", y = unname(metric_labels[[var]]), caption = paste0("Lines show GAMM fits; shaded ribbons indicate 95% confidence intervals.\n",
        "Three-month overlapping windows; step = 1 month; ASV prevalence threshold = 20%.")) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1, alpha = 1,
      size = 2)), fill = "none") + theme_cns_closed() + theme(axis.text.x = element_text(angle = 45,
    hjust = 1, vjust = 1, size = if (var %in% c("Degree", "Robustness_weight", "Density"))
      16.5
    else 11.2, margin = margin(t = 4)), axis.title.y = element_text(size = if (var %in% c("Degree",
    "Robustness_weight", "Density"))
    17.5
  else 9.5, face = "bold", margin = margin(r = 10)), legend.position = "top", legend.direction = "horizontal",
    legend.justification = "left")
  save_path <- paste0(PLOT_DIR, "/", var, "_Trend_ConF_USF.png")
  ggsave(save_path, plot = p, width = 7.2, height = 5.9, dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", save_path), plot = p, width = 7.2, height = 5.9, device = cairo_pdf,
    bg = "white")
}

cat("\n[5/5] 任务全部完成!\n")

cat(">>> 已删除整体网络 All，仅保留 ConF 和 USF。\n")

cat(">>> 参数: 3 个月窗口，步长 1 个月，节点出现率阈值 20%。\n")

cat(">>> 拓扑指标表已保存至: ", out_file, "\n", sep = "")

cat(">>> 参数记录已保存至: ", file.path(OUTPUT_DIR, "Network_Analysis_Parameters.csv"), "\n",
  sep = "")

cat(">>> GAMM统计表已保存至: ", file.path(OUTPUT_DIR, "GAMM_Network_Model_Statistics.csv"), "\n",
  sep = "")

cat(">>> 趋势图（600 dpi PNG/PDF）已保存至: ", PLOT_DIR, "\n\n", sep = "")
