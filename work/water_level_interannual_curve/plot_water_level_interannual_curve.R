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

project_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

project_library <- file.path(project_dir, "99_运行环境", "R_library")

if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)

input_xlsx <- if (length(args) >= 1L) args[[1]] else file.path(project_dir, "00_输入文件", "03_环境数据",
  "2024至2026年丹江口水位线.xlsx")

output_dir <- if (length(args) >= 2L) args[[2]] else file.path(script_dir, "output")

qnet_csv <- if (length(args) >= 3L) args[[3]] else file.path(project_dir, "outputs", "20260825_env_hydro_air_update",
  "tables", "hydrology_monthly_source_24.csv")

if (!file.exists(input_xlsx)) {
  stop("找不到输入文件：", input_xlsx)
}

if (!file.exists(qnet_csv)) {
  stop("找不到 Qnet 月表：", qnet_csv)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "ragg", "xml2", "cowplot")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))

read_water_level_xlsx <- function(path) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    raw <- readxl::read_excel(path, sheet = 1, col_names = c("Date", "Water_level_m"), skip = 1)
    return(data.frame(Date_raw = raw[[1]], Water_level_m = suppressWarnings(as.numeric(raw[[2]]))))
  }
  extract_dir <- tempfile("water_level_xlsx_")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
  archive_files <- utils::unzip(path, list = TRUE)$Name
  required_file <- "xl/worksheets/sheet1.xml"
  if (!required_file %in% archive_files) {
    stop("Excel 中未找到第一个工作表。")
  }
  files_to_extract <- required_file
  if ("xl/sharedStrings.xml" %in% archive_files) {
    files_to_extract <- c(files_to_extract, "xl/sharedStrings.xml")
  }
  utils::unzip(path, files = files_to_extract, exdir = extract_dir)
  shared_strings <- character(0)
  shared_path <- file.path(extract_dir, "xl", "sharedStrings.xml")
  if (file.exists(shared_path)) {
    shared_doc <- xml2::read_xml(shared_path)
    xml2::xml_ns_strip(shared_doc)
    shared_strings <- vapply(xml2::xml_find_all(shared_doc, ".//si"), xml2::xml_text, character(1))
  }
  sheet_path <- file.path(extract_dir, "xl", "worksheets", "sheet1.xml")
  sheet_doc <- xml2::read_xml(sheet_path)
  xml2::xml_ns_strip(sheet_doc)
  cells <- xml2::xml_find_all(sheet_doc, ".//sheetData/row/c")
  cell_value <- function(cell) {
    cell_type <- xml2::xml_attr(cell, "t")
    if (!is.na(cell_type) && cell_type == "inlineStr") {
      return(xml2::xml_text(xml2::xml_find_first(cell, "./is")))
    }
    value_node <- xml2::xml_find_first(cell, "./v")
    if (inherits(value_node, "xml_missing"))
      return(NA_character_)
    value <- xml2::xml_text(value_node)
    if (!is.na(cell_type) && cell_type == "s") {
      index <- suppressWarnings(as.integer(value)) + 1L
      if (!is.na(index) && index >= 1L && index <= length(shared_strings)) {
        return(shared_strings[[index]])
      }
    }
    value
  }
  refs <- xml2::xml_attr(cells, "r")
  rows <- suppressWarnings(as.integer(sub("^[A-Z]+", "", refs)))
  cols <- sub("[0-9]+$", "", refs)
  values <- vapply(cells, cell_value, character(1))
  date_rows <- rows[cols == "A" & rows >= 2L]
  level_rows <- rows[cols == "B" & rows >= 2L]
  complete_rows <- sort(intersect(date_rows, level_rows))
  date_values <- values[match(paste0("A", complete_rows), refs)]
  level_values <- values[match(paste0("B", complete_rows), refs)]
  data.frame(Date_raw = suppressWarnings(as.numeric(date_values)), Water_level_m = suppressWarnings(as.numeric(level_values)))
}

water_level <- read_water_level_xlsx(input_xlsx)

convert_to_date <- function(x) {
  if (inherits(x, "Date"))
    return(x)
  if (inherits(x, "POSIXt"))
    return(as.Date(x))
  if (is.numeric(x))
    return(as.Date(x, origin = "1899-12-30"))
  text <- trimws(as.character(x))
  text <- gsub("/", "-", text, fixed = TRUE)
  text[nchar(text) == 7L] <- paste0(text[nchar(text) == 7L], "-01")
  as.Date(text)
}

water_level$Date <- convert_to_date(water_level$Date_raw)

water_level <- water_level[!is.na(water_level$Date) & !is.na(water_level$Water_level_m), c("Date", "Water_level_m")]

water_level <- water_level[order(water_level$Date), , drop = FALSE]

row.names(water_level) <- NULL

if (nrow(water_level) < 2L) {
  stop("至少需要两个完整的日期—水位观测值。")
}

if (anyDuplicated(water_level$Date)) {
  stop("发现重复日期，请先确认每个日期只对应一个水位点。")
}

water_level$Month_start <- as.Date(format(water_level$Date, "%Y-%m-01"))

water_level$Month_label <- format(water_level$Month_start, "%Y-%m")

if (anyDuplicated(water_level$Month_label)) {
  stop("同一个月份存在多个水位观测值，无法建立一月一区间的横轴。")
}

axis_months <- data.frame(Month_start = seq(min(water_level$Month_start), max(water_level$Month_start),
  by = "month"))

axis_months$Month_label <- format(axis_months$Month_start, "%Y-%m")

axis_months$Point <- seq_len(nrow(axis_months))

axis_months_display <- axis_months

water_level$Point <- axis_months$Point[match(water_level$Month_label, axis_months$Month_label)]

interval_boundaries <- seq(0.5, nrow(axis_months) + 0.5, by = 1)

full_axis_data <- axis_months

full_axis_data$Water_level_m <- water_level$Water_level_m[match(full_axis_data$Month_label, water_level$Month_label)]

full_axis_data$Observed <- !is.na(full_axis_data$Water_level_m)

write.csv(full_axis_data, file.path(output_dir, "water_level_full_month_axis.csv"), row.names = FALSE,
  na = "")

qnet_raw <- read.csv(qnet_csv, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("HydroMonth", "Qnet") %in% names(qnet_raw))) {
  stop("Qnet 月表必须包含 HydroMonth 和 Qnet 两列。")
}

qnet_raw <- data.frame(Month_label = trimws(as.character(qnet_raw$HydroMonth)), Net_inflow_m3_s = -suppressWarnings(as.numeric(qnet_raw$Qnet)))

qnet_raw <- qnet_raw[!is.na(qnet_raw$Net_inflow_m3_s), , drop = FALSE]

if (anyDuplicated(qnet_raw$Month_label)) {
  stop("Qnet 月表中存在重复月份。")
}

change_data <- merge(axis_months, qnet_raw, by = "Month_label", all.x = TRUE, sort = FALSE)

change_data <- change_data[match(axis_months$Month_label, change_data$Month_label), ]

if (anyNA(change_data$Net_inflow_m3_s)) {
  stop("以下水位月份缺少 Qnet：", paste(change_data$Month_label[is.na(change_data$Net_inflow_m3_s)],
    collapse = ", "))
}

change_data$Flow_state <- ifelse(change_data$Net_inflow_m3_s > 0, "Net inflow", ifelse(change_data$Net_inflow_m3_s <
  0, "Net outflow", "Balanced"))

change_data$Flow_state <- factor(change_data$Flow_state, levels = c("Net inflow", "Net outflow", "Balanced"))

write.csv(change_data, file.path(output_dir, "net_inflow_source_data_aligned.csv"), row.names = FALSE)

smooth_curve <- as.data.frame(stats::spline(x = water_level$Point, y = water_level$Water_level_m, xout = seq(min(water_level$Point),
  max(water_level$Point), length.out = nrow(axis_months) * 40L), method = "natural"))

y_step <- 5

y_lower <- floor(min(water_level$Water_level_m)/y_step) * y_step

y_upper <- ceiling(max(water_level$Water_level_m)/y_step) * y_step

y_breaks <- seq(y_lower, y_upper, by = y_step)

line_color <- "#86B7E3"

water_level_plot <- ggplot() + geom_vline(xintercept = interval_boundaries, color = "#E6E6E6", linewidth = 0.35) +
  geom_line(data = smooth_curve, aes(x = x, y = y, color = "Water level"), linewidth = 0.85, lineend = "round") +
  geom_point(data = water_level, aes(x = Point, y = Water_level_m, fill = "Water level"), size = 2.5,
    shape = 21, color = "#1A1A1A", stroke = 0.35) + scale_color_manual(values = c(`Water level` = line_color),
  name = NULL) + scale_fill_manual(values = c(`Water level` = line_color), name = NULL) + scale_x_continuous(breaks = axis_months_display$Point,
  labels = axis_months_display$Month_label, limits = range(interval_boundaries), expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, labels = function(x) sprintf("%.2f", x), limits = c(y_lower,
    y_upper + 0.5), expand = c(0, 0)) + labs(x = NULL, y = "Water level (m)") + theme_classic(base_size = 13.5,
  base_family = "sans") + theme(axis.line = element_blank(), panel.border = element_rect(color = "#BFBFBF",
  fill = NA, linewidth = 0.45), axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "#BFBFBF",
  linewidth = 0.4), axis.text = element_text(color = "#595959"), axis.text.x = element_text(angle = 60,
  hjust = 1, vjust = 1, size = 11.5), axis.text.y = element_text(size = 11.5), axis.title.y = element_text(color = "#333333",
  margin = margin(r = 8)), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "#E6E6E6",
  linewidth = 0.4), panel.grid.minor = element_blank(), legend.position = "right", legend.text = element_text(color = "#595959",
  size = 10), legend.key = element_blank(), plot.margin = margin(t = 10, r = 18, b = 8, l = 10))

qnet_breaks <- pretty(c(change_data$Net_inflow_m3_s, 0), n = 7)

qnet_lower <- min(qnet_breaks)

qnet_upper <- max(qnet_breaks)

flow_colors <- c(`Net inflow` = "#4393C3", `Net outflow` = "#D98C45", Balanced = "#8C8C8C")

qnet_plot <- ggplot(change_data) + geom_vline(xintercept = interval_boundaries, color = "#E6E6E6", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "#333333", linewidth = 0.6) + geom_col(aes(x = Point, y = Net_inflow_m3_s,
  fill = Flow_state), width = 0.72, color = "#1A1A1A", linewidth = 0.35) + scale_fill_manual(values = flow_colors,
  breaks = c("Net inflow", "Net outflow"), drop = FALSE, name = NULL) + scale_x_continuous(breaks = axis_months_display$Point,
  labels = axis_months_display$Month_label, limits = range(interval_boundaries), expand = c(0, 0)) +
  scale_y_continuous(breaks = qnet_breaks, labels = function(x) format(x, scientific = FALSE, big.mark = ",",
    trim = TRUE), limits = c(qnet_lower, qnet_upper), expand = c(0, 0)) + labs(x = NULL, y = expression("Net inflow" ~
  (m^3 ~ s^{
    -1
  }))) + theme_classic(base_size = 13.5, base_family = "sans") + theme(axis.line = element_blank(),
  panel.border = element_rect(color = "#BFBFBF", fill = NA, linewidth = 0.45), axis.ticks.x = element_blank(),
  axis.ticks.y = element_line(color = "#BFBFBF", linewidth = 0.4), axis.text = element_text(color = "#595959"),
  axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 11.5), axis.text.y = element_text(size = 11.5),
  axis.title.y = element_text(color = "#333333", margin = margin(r = 8)), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.4, linetype = "dashed"), panel.grid.minor = element_blank(),
  legend.position = "right", legend.text = element_text(color = "#595959", size = 10), legend.key = element_blank(),
  plot.margin = margin(t = 10, r = 18, b = 8, l = 10))

aligned_standalone <- cowplot::align_plots(water_level_plot, qnet_plot, align = "v", axis = "lr")

water_level_aligned <- cowplot::ggdraw(aligned_standalone[[1]])

qnet_aligned <- cowplot::ggdraw(aligned_standalone[[2]])

qnet_for_combining <- qnet_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  plot.margin = margin(t = 10, r = 18, b = 6, l = 10))

water_level_plot_for_combining <- water_level_plot + theme(plot.margin = margin(t = 6, r = 18, b = 8,
  l = 10))

aligned_combined <- cowplot::align_plots(qnet_for_combining, water_level_plot_for_combining, align = "v",
  axis = "lr")

combined_plot <- cowplot::plot_grid(aligned_combined[[1]], aligned_combined[[2]], ncol = 1, rel_heights = c(0.72,
  1.38))

line_png <- file.path(output_dir, "water_level_interannual_scatter_curve.png")

line_pdf <- file.path(output_dir, "water_level_interannual_scatter_curve.pdf")

bar_png <- file.path(output_dir, "net_inflow_bar_aligned.png")

bar_pdf <- file.path(output_dir, "net_inflow_bar_aligned.pdf")

combined_png <- file.path(output_dir, "water_level_and_net_inflow_combined.png")

combined_pdf <- file.path(output_dir, "water_level_and_net_inflow_combined.pdf")

ggsave(line_png, water_level_aligned, width = 9.6, height = 4.8, dpi = 600, bg = "white", device = ragg::agg_png)

ggsave(line_pdf, water_level_aligned, width = 9.6, height = 4.8, bg = "white", device = grDevices::cairo_pdf)

ggsave(bar_png, qnet_aligned, width = 9.6, height = 3.8, dpi = 600, bg = "white", device = ragg::agg_png)

ggsave(bar_pdf, qnet_aligned, width = 9.6, height = 3.8, bg = "white", device = grDevices::cairo_pdf)

ggsave(combined_png, combined_plot, width = 9.6, height = 7.6, dpi = 600, bg = "white", device = ragg::agg_png)

ggsave(combined_pdf, combined_plot, width = 9.6, height = 7.6, bg = "white", device = grDevices::cairo_pdf)

message("水位曲线图：", normalizePath(line_png, winslash = "/"))

message("净流入柱状图：", normalizePath(bar_png, winslash = "/"))

message("上下拼接图：", normalizePath(combined_png, winslash = "/"))
