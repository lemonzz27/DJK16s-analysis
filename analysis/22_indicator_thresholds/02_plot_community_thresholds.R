options(stringsAsFactors = FALSE, scipen = 999)

try(Sys.setlocale("LC_ALL", ""), silent = TRUE)

ROOT_DIR <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT_DIR <- file.path(ROOT_DIR, "01_分析项目", "36_TITAN2致病菌ASV水平_水位")

TABLE_DIR <- file.path(PROJECT_DIR, "结果表格", "WL", "ASV")

FIGURE_ROOT <- file.path(PROJECT_DIR, "结果图片", "WL", "ASV")

PLOTDATA_DIR <- file.path(TABLE_DIR, "plotdata_figure1_community_sumz")

FIGURE_DIR <- file.path(FIGURE_ROOT, "01_community_sumz_ASV")

dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(unique(c(file.path(ROOT_DIR, "99_运行环境", "R_library"), .libPaths())))

required_packages <- c("ggplot2", "dplyr", "tidyr", "ragg", "svglite")

available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!available)) {
  stop("缺少R包：", paste(required_packages[!available], collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

read_csv <- function(name) {
  path <- file.path(TABLE_DIR, name)
  if (!file.exists(path))
    stop("缺少正式TITAN2结果表：", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

save_four_formats <- function(plot, prefix, width, height) {
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(paste0(prefix, ".tiff"), plot, width = width, height = height, dpi = 600, bg = "white", compression = "lzw",
    device = ragg::agg_tiff)
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height, bg = "white", device = svglite::svglite)
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, bg = "white", device = cairo_pdf)
}

config <- read_csv("04_analysis_config.csv")

config_value <- function(parameter) {
  x <- config$value[match(parameter, config$parameter)]
  if (length(x) != 1L || is.na(x))
    stop("配置表缺少参数：", parameter)
  as.character(x)
}

if (config_value("run_mode") != "formal" || as.integer(config_value("num_permutations")) != 500L || as.integer(config_value("n_bootstraps")) !=
  500L || config_value("taxon_level") != "ASV" || toupper(config_value("taxonomy_annotation_used")) !=
  "FALSE") {
  stop("绘图只允许使用未用分类注释的ASV级500/500正式结果。")
}

screening <- read_csv("05_all_candidate_ASV_input_screening.csv")

analyzed <- read_csv("10_all_analyzed_ASV_TITAN2_results.csv")

robust <- read_csv("11_robust_indicator_ASVs.csv")

profile <- read_csv("15_community_sumz_profile_by_split.csv")

thresholds <- read_csv("16_community_thresholds.csv")

n_screened_pool <- nrow(screening)

n_analyzed <- nrow(analyzed)

n_robust <- nrow(robust)

if (n_screened_pool != 2755L || n_analyzed != 476L || sum(as.logical(screening$included_in_TITAN2)) !=
  476L) {
  stop("ASV筛选池或正式入模ASV数不符合审计基线。")
}

if (n_robust != sum(as.logical(analyzed$robust_indicator))) {
  stop("稳健ASV表与全部入模ASV结果不一致。")
}

if (!setequal(as.character(robust$taxon_id), as.character(analyzed$taxon_id[as.logical(analyzed$robust_indicator)]))) {
  stop("稳健ASV表与全部结果中的稳健ASV ID集合不一致。")
}

required_profile <- c("change_point", "sumz_negative_all_taxa", "sumz_positive_all_taxa", "sumz_negative_robust_filtered",
  "sumz_positive_robust_filtered")

required_threshold <- c("threshold_type", "taxa_scope", "direction", "change_point_observed", "change_point_median",
  "robust_taxa_n")

if (!all(required_profile %in% names(profile))) {
  stop("community sum(z)原始表缺少必要列。")
}

if (!all(required_threshold %in% names(thresholds))) {
  stop("community阈值表缺少必要列。")
}

profile_long <- dplyr::filter(dplyr::mutate(tidyr::pivot_longer(dplyr::select(profile, dplyr::all_of(required_profile)),
  cols = dplyr::starts_with("sumz_"), names_to = "series", values_to = "sumz"), direction = ifelse(grepl("negative",
  series), "negative_z-", "positive_z+"), taxa_scope = ifelse(grepl("robust_filtered", series), "robust_filtered",
  "all_taxa"), change_point = as.numeric(change_point), sumz = as.numeric(sumz)), is.finite(change_point),
  is.finite(sumz))

profile_unique <- dplyr::arrange(dplyr::summarise(dplyr::group_by(profile_long, change_point, direction,
  taxa_scope), sumz = max(sumz), n_tied_splits = dplyr::n(), .groups = "drop"), direction, taxa_scope,
  change_point)

make_shape_preserving_curve <- function(data, n_grid = 700L) {
  data <- data[order(data$change_point), , drop = FALSE]
  if (nrow(data) < 3L || anyDuplicated(data$change_point)) {
    stop("每条曲线至少需要3个唯一候选切点。")
  }
  interpolation <- splinefun(x = data$change_point, y = data$sumz, method = "monoH.FC")
  x_grid <- sort(unique(c(seq(min(data$change_point), max(data$change_point), length.out = n_grid),
    data$change_point)))
  curve_y <- as.numeric(interpolation(x_grid))
  interval_index <- findInterval(x_grid, data$change_point, all.inside = TRUE)
  interval_index[interval_index >= nrow(data)] <- nrow(data) - 1L
  endpoint_low <- pmin(data$sumz[interval_index], data$sumz[interval_index + 1L])
  endpoint_high <- pmax(data$sumz[interval_index], data$sumz[interval_index + 1L])
  curve_y <- pmin(pmax(curve_y, endpoint_low), endpoint_high)
  data.frame(change_point = x_grid, sumz = curve_y, direction = data$direction[1L], taxa_scope = data$taxa_scope[1L],
    stringsAsFactors = FALSE)
}

profile_curves <- dplyr::mutate(dplyr::bind_rows(lapply(dplyr::group_split(dplyr::group_by(profile_unique,
  direction, taxa_scope)), make_shape_preserving_curve)), profile_group = interaction(direction, taxa_scope,
  drop = TRUE), ribbon_ymin = pmin(sumz, 0), ribbon_ymax = pmax(sumz, 0))

robust_thresholds <- dplyr::transmute(dplyr::filter(thresholds, threshold_type %in% c("fsumz-", "fsumz+"),
  taxa_scope == "robust_filtered", direction %in% c("negative_z-", "positive_z+"), robust_taxa_n >
    0, is.finite(as.numeric(change_point_observed))), threshold_type, direction, taxa_scope, threshold = as.numeric(change_point_observed),
  threshold_median = as.numeric(change_point_median), robust_taxa_n = as.integer(robust_taxa_n))

if (nrow(robust_thresholds) == 0L) {
  stop("正式分析没有可绘制的稳健ASV群落阈值。")
}

if (any(table(robust_thresholds$direction) > 1L)) {
  stop("每个响应方向只能有一个观察稳健群落阈值。")
}

for (i in seq_len(nrow(robust_thresholds))) {
  direction_i <- robust_thresholds$direction[i]
  column_i <- if (direction_i == "negative_z-") {
    "sumz_negative_robust_filtered"
  }
  else {
    "sumz_positive_robust_filtered"
  }
  raw_y <- as.numeric(profile[[column_i]])
  if (!any(is.finite(raw_y))) {
    stop("原始filtered sum(z)曲线没有有限值。")
  }
  peak_x <- as.numeric(profile$change_point)[abs(raw_y - max(raw_y, na.rm = TRUE)) < 1e-10]
  if (!any(abs(peak_x - robust_thresholds$threshold[i]) < 1e-08)) {
    stop("观察到的稳健群落阈值不是原始filtered sum(z)峰。")
  }
}

threshold_heights <- dplyr::bind_rows(lapply(seq_len(nrow(robust_thresholds)), function(i) {
  row_i <- robust_thresholds[i, ]
  curve_i <- profile_curves[profile_curves$direction == row_i$direction & profile_curves$taxa_scope ==
    "robust_filtered", , drop = FALSE]
  row_i$sumz_at_threshold <- approx(x = curve_i$change_point, y = curve_i$sumz, xout = row_i$threshold,
    rule = 2)$y
  row_i
}))

plot_y_range <- diff(range(c(0, profile_curves$sumz), finite = TRUE))

label_offset <- max(0.035 * plot_y_range, 0.5)

threshold_heights <- dplyr::mutate(threshold_heights, segment_ymin = pmin(0, sumz_at_threshold), segment_ymax = pmax(0,
  sumz_at_threshold), threshold_label = sprintf("%.2f m", threshold + 1e-08), label_y = sumz_at_threshold +
  label_offset)

write_csv(profile, file.path(PLOTDATA_DIR, "01_raw_community_sumz_profile_by_split.csv"))

write_csv(profile_unique, file.path(PLOTDATA_DIR, "02_unique_candidate_change_point_profile.csv"))

write_csv(profile_curves, file.path(PLOTDATA_DIR, "03_shape_preserving_curve_plotdata.csv"))

write_csv(thresholds, file.path(PLOTDATA_DIR, "04_raw_community_thresholds.csv"))

write_csv(threshold_heights, file.path(PLOTDATA_DIR, "05_plotted_robust_thresholds.csv"))

write_csv(data.frame(parameter = c("candidate_ASV_screening_pool", "analyzed_ASV_after_screening", "robust_ASV",
  "num_permutations", "n_bootstraps", "smoothing_method", "threshold_source", "label_offset_sumz_units"),
  value = c(n_screened_pool, n_analyzed, n_robust, 500, 500, "shape_preserving_monoH.FC_display_only",
    "observed_TITAN2_filtered_sumz_raw_maximum", label_offset), stringsAsFactors = FALSE), file.path(PLOTDATA_DIR,
  "06_figure_parameters_and_QA.csv"))

write_csv(data.frame(file = c("01_raw_community_sumz_profile_by_split.csv", "02_unique_candidate_change_point_profile.csv",
  "03_shape_preserving_curve_plotdata.csv", "04_raw_community_thresholds.csv", "05_plotted_robust_thresholds.csv",
  "06_figure_parameters_and_QA.csv"), purpose_cn = c("TITAN2输出的全部候选切点及四条原始sum(z)序列",
  "同一候选切点按每条序列取最大值后的包络数据", "图中四条形状保持平滑曲线的逐点数据，可直接重画",
  "四类总体阈值、bootstrap中位数及5%至95%分位数", "图中黑色竖虚线和标签的准确位置、方向、峰高及稳健ASV数",
  "正式运行、ASV数量、平滑和阈值来源的绘图QA"), stringsAsFactors = FALSE), file.path(PLOTDATA_DIR,
  "00_plotdata_guide.csv"))

direction_colors <- c(`negative_z-` = "#2C7FB8", `positive_z+` = "#C62F43")

direction_labels <- c(`negative_z-` = "Negative (z-)", `positive_z+` = "Positive (z+)")

scope_linetypes <- c(all_taxa = "22", robust_filtered = "solid")

scope_alphas <- c(all_taxa = 0.27, robust_filtered = 1)

scope_labels <- c(all_taxa = sprintf("All analyzed ASVs after screening (n = %d)", n_analyzed), robust_filtered = sprintf("Robust ASVs (purity & reliability >= 0.95; n = %d)",
  n_robust))

robust_ribbons <- dplyr::filter(profile_curves, taxa_scope == "robust_filtered")

sumz_plot <- ggplot() + geom_hline(yintercept = 0, color = "#9A9A9A", linewidth = 0.38) + geom_ribbon(data = robust_ribbons,
  aes(x = change_point, ymin = ribbon_ymin, ymax = ribbon_ymax, group = direction, fill = direction),
  alpha = 0.18, colour = NA, show.legend = FALSE) + geom_line(data = profile_curves, aes(x = change_point,
  y = sumz, group = profile_group, colour = direction, linetype = taxa_scope, alpha = taxa_scope),
  linewidth = 0.86, lineend = "round") + geom_segment(data = threshold_heights, aes(x = threshold,
  xend = threshold, y = segment_ymin, yend = segment_ymax), colour = "black", linewidth = 0.4, linetype = "22",
  lineend = "round") + geom_text(data = threshold_heights, aes(x = threshold, y = label_y, label = threshold_label,
  colour = direction), fontface = "bold", size = 3, show.legend = FALSE) + scale_colour_manual(values = direction_colors,
  breaks = names(direction_labels), labels = direction_labels, name = "Direction") + scale_fill_manual(values = direction_colors) +
  scale_linetype_manual(values = scope_linetypes, breaks = names(scope_labels), labels = scope_labels,
    name = "ASV scope") + scale_alpha_manual(values = scope_alphas, guide = "none") + scale_x_continuous(breaks = seq(150,
  170, 5), expand = expansion(mult = c(0.015, 0.025))) + scale_y_continuous(expand = expansion(mult = c(0.06,
  0.12))) + labs(title = "Community-level TITAN2 sum(z) profiles", subtitle = paste0("Candidate-pathogen ASVs along the water-level gradient; ",
  "black dashed lines mark observed robust thresholds"), x = "Water level (m)", y = "Community sum(z)") +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", alpha = 1, linewidth = 1.2)),
    linetype = guide_legend(order = 2, override.aes = list(colour = "#4A4A4A", alpha = c(0.32, 1),
      linewidth = c(0.9, 1.15)))) + theme_classic(base_family = "sans", base_size = 12.5) + theme(panel.border = element_rect(fill = NA,
  colour = "black", linewidth = 0.4), axis.line = element_blank(), axis.ticks = element_line(colour = "black",
  linewidth = 0.42), axis.ticks.length = grid::unit(1.5, "mm"), axis.text = element_text(colour = "black",
  size = 10.8), axis.title = element_text(colour = "black", size = 12, face = "bold"), panel.grid.major.x = element_line(colour = "#E6E6E6",
  linewidth = 0.3), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold",
  size = 14.2, hjust = 0), plot.subtitle = element_text(colour = "#4D4D4D", size = 10.3, hjust = 0,
  margin = margin(b = 5)), legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
  legend.background = element_rect(fill = "transparent", colour = NA), legend.key = element_rect(fill = "transparent",
    colour = NA), legend.title = element_text(face = "bold", size = 10.5), legend.text = element_text(size = 9.8),
  legend.key.width = grid::unit(10, "mm"), legend.key.height = grid::unit(3.6, "mm"), plot.background = element_rect(fill = "white",
    colour = NA), panel.background = element_rect(fill = "white", colour = NA), plot.margin = margin(8,
    10, 7, 8))

output_prefix <- file.path(FIGURE_DIR, "TITAN2_WL_ASV_community_sumz_smooth")

save_four_formats(sumz_plot, output_prefix, width = 8.6, height = 6)

expected_figure_files <- paste0(output_prefix, c(".png", ".tiff", ".svg", ".pdf"))

if (any(!file.exists(expected_figure_files)) || any(file.info(expected_figure_files)$size <= 0)) {
  stop("第一张图的四种输出格式未全部成功保存。")
}

capture.output(sessionInfo(), file = file.path(PLOTDATA_DIR, "sessionInfo_plot.txt"))

message("完成ASV级第一张图：入模ASV ", n_analyzed, "；稳健ASV ", n_robust, "；观察稳健阈值 ",
  paste(sprintf("%s %.3f m", robust_thresholds$direction, robust_thresholds$threshold), collapse = "; "),
  "。")
