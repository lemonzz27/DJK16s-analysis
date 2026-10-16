rm(list = ls())

if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"), silent = TRUE)
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

script_dir <- if (length(script_path) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_path[[1]]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

args <- commandArgs(trailingOnly = TRUE)

input_xlsx <- if (length(args) >= 1L) {
  args[[1]]
} else {
  file.path(script_dir, "环境因子_按月份平均_17个月.xlsx")
}

output_dir <- if (length(args) >= 2L) {
  args[[2]]
} else {
  file.path(script_dir, "output_monthly_environment")
}

if (!file.exists(input_xlsx)) {
  stop("找不到输入文件：", input_xlsx, "\n请把 环境因子_按月份平均_17个月.xlsx 放到脚本同目录，",
    "或在命令行把它作为第 1 个参数传入。")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(output_dir, "individual"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(output_dir, "pages"), recursive = TRUE, showWarnings = FALSE)

input_is_csv <- tolower(tools::file_ext(input_xlsx)) == "csv"

required_packages <- c("ggplot2", "ragg", "cowplot")

if (!input_is_csv) required_packages <- c(required_packages, "readxl")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))

if (input_is_csv) {
  input_sheet <- "CSV"
  raw <- read.csv(input_xlsx, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  sheets <- readxl::excel_sheets(input_xlsx)
  input_sheet <- if ("按月平均" %in% sheets)
    "按月平均"
  else sheets[[1]]
  raw <- readxl::read_excel(input_xlsx, sheet = input_sheet)
  raw <- as.data.frame(raw, check.names = FALSE, stringsAsFactors = FALSE)
}

if (!"Month" %in% names(raw)) {
  stop("输入表必须包含 Month 列，格式例如 2025-01。")
}

raw <- raw[, !grepl("^Point($|\\.)", names(raw)), drop = FALSE]

raw$Month <- trimws(as.character(raw$Month))

raw$Month_start <- as.Date(paste0(raw$Month, "-01"), format = "%Y-%m-%d")

if (anyNA(raw$Month_start)) {
  bad <- raw$Month[is.na(raw$Month_start)]
  stop("以下 Month 无法识别为 YYYY-MM：", paste(bad, collapse = ", "))
}

if (anyDuplicated(raw$Month)) {
  stop("Month 列存在重复月份，请确保每个月只有一行月均值。")
}

raw <- raw[order(raw$Month_start), , drop = FALSE]

row.names(raw) <- NULL

axis_months <- data.frame(Month_start = seq(min(raw$Month_start), max(raw$Month_start), by = "month"))

axis_months$Month <- format(axis_months$Month_start, "%Y-%m")

axis_months$Point <- seq_len(nrow(axis_months))

axis_months_display <- axis_months[seq(1L, nrow(axis_months), by = 2L), , drop = FALSE]

plot_data <- merge(axis_months, raw[, setdiff(names(raw), "Month_start"), drop = FALSE], by = "Month",
  all.x = TRUE, sort = FALSE)

plot_data <- plot_data[match(axis_months$Month, plot_data$Month), , drop = FALSE]

plot_data$Month_start <- axis_months$Month_start

plot_data$Point <- axis_months$Point

interval_boundaries <- seq(0.5, nrow(axis_months) + 0.5, by = 1)

write.csv(plot_data, file.path(output_dir, "monthly_environment_full_axis.csv"), row.names = FALSE, na = "")

exclude_columns <- c("Month", "Month_start", "Point", "N", "Longitude", "Latitude")

metrics_to_plot <- setdiff(names(plot_data), exclude_columns)

metrics_to_plot <- metrics_to_plot[metrics_to_plot %in% names(plot_data)]

if (length(metrics_to_plot) == 0L) {
  stop("在输入表中没有找到可绘制的指标列。")
}

for (metric in metrics_to_plot) {
  plot_data[[metric]] <- suppressWarnings(as.numeric(plot_data[[metric]]))
}

y_labels <- c(WT = "Water temperature (°C)", pH = "pH", DO = "DO (mg L⁻¹)", CODmn = "CODMn (mg L⁻¹)",
  BOD5 = "BOD₅ (mg L⁻¹)", NH3N = "NH₃-N (mg L⁻¹)", TN = "TN (mg L⁻¹)", TP = "TP (mg L⁻¹)",
  Cu = "Cu (mg L⁻¹)", Zn = "Zn (mg L⁻¹)", F = "F (mg L⁻¹)", Se = "Se (mg L⁻¹)", As = "As (mg L⁻¹)",
  Pb = "Pb (mg L⁻¹)", Chla = "Chl-a (mg L⁻¹)", Q = "Q (m³ s⁻¹)", WL = "Water level (m)",
  AIRT = "Air temperature (°C)", MWL = "MWL", RSV = "RSV", Qin = "Qin (m³ s⁻¹)", Qout = "Qout (m³ s⁻¹)",
  Qnet = "Qnet (m³ s⁻¹)", NFR = "NFR", HRT = "HRT (d)")

plot_titles <- c(WT = "Water temperature", pH = "pH", DO = "Dissolved oxygen", CODmn = "CODMn", BOD5 = "BOD₅",
  NH3N = "NH₃-N", TN = "Total nitrogen", TP = "Total phosphorus", Cu = "Copper", Zn = "Zinc", F = "Fluoride",
  Se = "Selenium", As = "Arsenic", Pb = "Lead", Chla = "Chlorophyll-a", Q = "Flow", WL = "Water level",
  AIRT = "Air temperature", MWL = "MWL", RSV = "RSV", Qin = "Qin", Qout = "Qout", Qnet = "Qnet", NFR = "NFR",
  HRT = "Hydraulic retention time")

metric_colors <- c(TN = "#4C78A8", Chla = "#59A14F", HRT = "#F28E2B", Qnet = "#E15759")

default_line_color <- "#4C78A8"

base_publication_theme <- theme_bw(base_size = 14, base_family = "sans") + theme(panel.grid.major = element_line(color = "#D9D9D9",
  linewidth = 0.35), panel.grid.minor = element_blank(), panel.border = element_rect(color = "#333333",
  fill = NA, linewidth = 0.55), axis.line = element_blank(), axis.ticks = element_line(color = "#666666",
  linewidth = 0.35), axis.ticks.x = element_blank(), axis.text = element_text(color = "#4D4D4D"), axis.text.x = element_text(angle = 45,
  hjust = 1, vjust = 1, size = 12.2), axis.text.y = element_text(size = 11.5), axis.title.x = element_text(color = "#333333"),
  axis.title.y = element_text(color = "#333333", margin = margin(r = 8)), plot.title = element_text(size = 17.5,
    face = "bold", hjust = 0.5, color = "#1A1A1A", margin = margin(b = 7)), plot.margin = margin(t = 16,
    r = 10, b = 8, l = 8))

make_metric_plot <- function(metric, show_x_text = TRUE) {
  df <- plot_data[, c("Month", "Point", metric), drop = FALSE]
  names(df)[3] <- "Value"
  line_color <- if (metric %in% names(metric_colors)) {
    unname(metric_colors[metric])
  }
  else {
    default_line_color
  }
  point_fill <- line_color
  observed <- df[is.finite(df$Value), c("Month", "Point", "Value"), drop = FALSE]
  if (nrow(observed) == 0L) {
    stop("指标 ", metric, " 没有可绘制的数值。")
  }
  if (nrow(observed) >= 2L) {
    smooth_curve <- as.data.frame(stats::spline(x = observed$Point, y = observed$Value, xout = seq(min(observed$Point),
      max(observed$Point), length.out = nrow(axis_months) * 40L), method = "natural"))
  }
  else {
    smooth_curve <- data.frame(x = observed$Point, y = observed$Value)
  }
  finite_values <- observed$Value
  curve_values <- smooth_curve$y[is.finite(smooth_curve$y)]
  y_range_values <- c(finite_values, curve_values)
  y_breaks <- pretty(y_range_values, n = 6)
  y_lower <- min(y_breaks)
  y_upper <- max(y_breaks)
  if (isTRUE(all.equal(y_lower, y_upper))) {
    y_lower <- y_lower - 0.5
    y_upper <- y_upper + 0.5
  }
  p <- ggplot() + geom_vline(xintercept = axis_months$Point, color = "#E6E6E6", linewidth = 0.35) +
    geom_line(data = smooth_curve, aes(x = x, y = y), linewidth = 1.05, color = line_color, lineend = "round") +
    geom_point(data = observed, aes(x = Point, y = Value), size = 2.7, shape = 21, fill = point_fill,
      color = "#2B2B2B", stroke = 0.4) + scale_x_continuous(breaks = axis_months_display$Point,
    labels = axis_months_display$Month, limits = range(interval_boundaries), expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, limits = c(y_lower, y_upper), expand = expansion(mult = c(0.03,
      0.06))) + labs(title = if (metric %in% names(plot_titles))
    unname(plot_titles[metric])
  else metric, x = NULL, y = if (metric %in% names(y_labels))
    unname(y_labels[metric])
  else metric) + base_publication_theme
  if (metric == "Qnet" && min(finite_values) < 0 && max(finite_values) > 0) {
    p <- p + geom_hline(yintercept = 0, color = "#333333", linewidth = 0.55, linetype = "dashed")
  }
  if (!show_x_text) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  p
}

plot_list <- setNames(vector("list", length(metrics_to_plot)), metrics_to_plot)

for (metric in metrics_to_plot) {
  p <- make_metric_plot(metric)
  plot_list[[metric]] <- p
  png_file <- file.path(output_dir, "individual", paste0(metric, "_monthly_curve.png"))
  pdf_file <- file.path(output_dir, "individual", paste0(metric, "_monthly_curve.pdf"))
  ggsave(png_file, p, width = 9.6, height = 4.8, dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(pdf_file, p, width = 9.6, height = 4.8, bg = "white", device = grDevices::cairo_pdf)
}

focus_metrics <- c("TN", "Chla", "HRT", "Qnet")

focus_metrics <- focus_metrics[focus_metrics %in% names(plot_list)]

if (length(focus_metrics) > 0L) {
  focus_plot <- cowplot::plot_grid(plotlist = plot_list[focus_metrics], ncol = 2, align = "hv", axis = "tblr")
  ggsave(file.path(output_dir, "TN_Chla_HRT_Qnet_combined.png"), focus_plot, width = 12, height = 8.2,
    dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, "TN_Chla_HRT_Qnet_combined.pdf"), focus_plot, width = 12, height = 8.2,
    bg = "white", device = grDevices::cairo_pdf)
}

page_size <- 4L

page_groups <- split(metrics_to_plot, ceiling(seq_along(metrics_to_plot)/page_size))

for (i in seq_along(page_groups)) {
  these_metrics <- page_groups[[i]]
  page_plot <- cowplot::plot_grid(plotlist = plot_list[these_metrics], ncol = 2, align = "hv", axis = "tblr",
    labels = LETTERS[seq_along(these_metrics)], label_size = 12)
  page_name <- sprintf("monthly_metrics_page_%02d", i)
  ggsave(file.path(output_dir, "pages", paste0(page_name, ".png")), page_plot, width = 12, height = 8.2,
    dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, "pages", paste0(page_name, ".pdf")), page_plot, width = 12, height = 8.2,
    bg = "white", device = grDevices::cairo_pdf)
}

message("读取工作表：", input_sheet)

message("输入文件：", normalizePath(input_xlsx, winslash = "/"))

message("输出目录：", normalizePath(output_dir, winslash = "/"))

message("已绘制指标：", paste(metrics_to_plot, collapse = ", "))

message("重点四联图：TN, Chla, HRT, Qnet（无 A/B/C/D 标记；4 个不同配色）")

message("缺测月份保留为空；绘图使用 natural spline 跨缺测月份连续连接，不补虚假观测点。")
