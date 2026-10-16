options(stringsAsFactors = FALSE, scipen = 999)

ROOT_DIR <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT_DIR <- file.path(ROOT_DIR, "01_分析项目", "34_蓝藻门结构方程模型")

INPUT_FILE <- file.path(PROJECT_DIR, "原始数据", "蓝藻门SEM_每样本原始数据_120样本.csv")

OUTPUT_DIR <- file.path(PROJECT_DIR, "结果图片", "05_蓝藻门月均相对丰度纵向曲线")

TABLE_DIR <- file.path(PROJECT_DIR, "结果表格", "05_蓝藻门月均相对丰度纵向曲线")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(unique(c(file.path(ROOT_DIR, "99_运行环境", "R_library"), .libPaths())))

required_packages <- c("ggplot2", "ragg", "svglite")

available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!available)) {
  stop("Missing R packages: ", paste(required_packages[!available], collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))

if (!file.exists(INPUT_FILE)) stop("Input file not found: ", INPUT_FILE)

raw <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- c("sample", "site", "month_code", "date", "size_fraction", "cyanobacteria_percent",
  "taxonomy_source")

if (!all(required_columns %in% names(raw))) {
  stop("Input data are missing required columns: ", paste(setdiff(required_columns, names(raw)), collapse = ", "))
}

if (nrow(raw) != 120L || anyDuplicated(raw$sample)) {
  stop("Expected 120 unique sample rows.")
}

raw$date <- as.Date(raw$date)

raw$cyanobacteria_percent <- as.numeric(raw$cyanobacteria_percent)

if (anyNA(raw$date) || any(!is.finite(raw$cyanobacteria_percent))) {
  stop("Date or cyanobacterial relative-abundance values are invalid.")
}

monthly_split <- split(raw, format(raw$date, "%Y-%m"))

monthly <- do.call(rbind, lapply(monthly_split, function(z) {
  data.frame(Month = format(z$date[[1L]], "%Y-%m"), Date = as.Date(format(z$date[[1L]], "%Y-%m-01")),
    N_samples = nrow(z), N_sites = length(unique(z$site)), N_ConF = sum(tolower(z$size_fraction) ==
      "large"), N_USF = sum(tolower(z$size_fraction) == "small"), Mean_relative_abundance_percent = mean(z$cyanobacteria_percent),
    SD_relative_abundance_percent = stats::sd(z$cyanobacteria_percent), SE_relative_abundance_percent = stats::sd(z$cyanobacteria_percent)/sqrt(nrow(z)),
    stringsAsFactors = FALSE)
}))

monthly <- monthly[order(monthly$Date), , drop = FALSE]

row.names(monthly) <- NULL

if (any(monthly$N_samples != 8L) || any(monthly$N_sites != 4L) || any(monthly$N_ConF != 4L) || any(monthly$N_USF !=
  4L)) {
  stop("Each observed month must contain 4 sites, 4 ConF samples, and 4 USF samples.")
}

axis_dates <- seq(min(monthly$Date), max(monthly$Date), by = "month")

axis_data <- data.frame(Date = axis_dates, Month = format(axis_dates, "%Y-%m"), Point = seq_along(axis_dates),
  stringsAsFactors = FALSE)

plot_data <- merge(axis_data, monthly, by = c("Date", "Month"), all.x = TRUE, sort = FALSE)

plot_data <- plot_data[match(axis_data$Month, plot_data$Month), , drop = FALSE]

plot_data$Point <- axis_data$Point

observed <- plot_data[is.finite(plot_data$Mean_relative_abundance_percent), , drop = FALSE]

smooth <- as.data.frame(stats::spline(x = observed$Point, y = observed$Mean_relative_abundance_percent,
  xout = seq(min(observed$Point), max(observed$Point), length.out = 800L), method = "natural"))

names(smooth) <- c("Point", "Mean_relative_abundance_percent")

write.csv(monthly, file.path(TABLE_DIR, "Cyanobacteria_monthly_mean_relative_abundance.csv"), row.names = FALSE,
  na = "", fileEncoding = "UTF-8")

write.csv(plot_data, file.path(TABLE_DIR, "Cyanobacteria_monthly_complete_axis_plot_data.csv"), row.names = FALSE,
  na = "", fileEncoding = "UTF-8")

line_color <- "#2C7FB8"

x_values <- c(observed$Mean_relative_abundance_percent, smooth$Mean_relative_abundance_percent)

x_breaks <- pretty(c(0, x_values), n = 5)

x_breaks <- x_breaks[x_breaks >= 0]

x_upper <- max(x_breaks)

if (x_upper < max(x_values)) {
  x_upper <- max(pretty(c(0, max(x_values)), n = 5))
}

p <- ggplot() + geom_hline(yintercept = axis_data$Point, color = "#E6E6E6", linewidth = 0.38) + geom_vline(xintercept = x_breaks,
  color = "#E6E6E6", linewidth = 0.38) + geom_path(data = smooth, aes(x = Mean_relative_abundance_percent,
  y = Point), color = line_color, linewidth = 1.2, lineend = "round") + geom_point(data = observed,
  aes(x = Mean_relative_abundance_percent, y = Point), size = 3.4, shape = 21, fill = line_color, color = "#2B2B2B",
  stroke = 0.48) + scale_x_continuous(breaks = x_breaks, limits = c(0, x_upper), expand = expansion(mult = c(0,
  0.04))) + scale_y_continuous(breaks = axis_data$Point, labels = axis_data$Month, limits = c(0.5,
  nrow(axis_data) + 0.5), expand = c(0, 0)) + labs(title = "Cyanobacterial relative\nabundance", x = "Cyanobacterial relative\nabundance (%)",
  y = NULL) + theme_bw(base_size = 14, base_family = "sans") + theme(panel.grid = element_blank(),
  panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.65), axis.line = element_blank(),
  axis.ticks = element_line(color = "#666666", linewidth = 0.42), axis.text.x = element_text(color = "#4D4D4D",
    size = 12.5), axis.text.y = element_text(color = "#4D4D4D", size = 12.8), axis.title.x = element_text(color = "#333333",
    size = 13.5, face = "bold", lineheight = 1, margin = margin(t = 9)), plot.title = element_text(size = 16,
    face = "bold", hjust = 0.5, color = "#1A1A1A", lineheight = 1, margin = margin(b = 10)), plot.margin = margin(t = 14,
    r = 12, b = 10, l = 10), plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white",
    color = NA))

prefix <- file.path(OUTPUT_DIR, "Cyanobacteria_relative_abundance_vertical_curve")

figure_width <- 3.6

figure_height <- 9.4

ggsave(paste0(prefix, ".png"), p, width = figure_width, height = figure_height, dpi = 600, bg = "white",
  device = ragg::agg_png)

ggsave(paste0(prefix, ".tiff"), p, width = figure_width, height = figure_height, dpi = 600, bg = "white",
  compression = "lzw", device = ragg::agg_tiff)

ggsave(paste0(prefix, ".svg"), p, width = figure_width, height = figure_height, bg = "white", device = svglite::svglite)

ggsave(paste0(prefix, ".pdf"), p, width = figure_width, height = figure_height, bg = "white", device = cairo_pdf)

qa <- data.frame(check = c("input_has_120_unique_samples", "taxonomy_is_SILVA_Cyanobacteriota_or_Cyanobacteria",
  "each_observed_month_has_8_samples", "each_observed_month_has_4_sites", "each_observed_month_has_4_ConF_and_4_USF",
  "observed_months_are_not_imputed", "figure_is_narrow_and_tall", "all_four_output_formats_exist"),
  passed = c(nrow(raw) == 120L && !anyDuplicated(raw$sample), all(grepl("Cyanobacteriota/Cyanobacteria",
    raw$taxonomy_source, fixed = TRUE)), all(monthly$N_samples == 8L), all(monthly$N_sites == 4L),
    all(monthly$N_ConF == 4L & monthly$N_USF == 4L), sum(is.finite(plot_data$Mean_relative_abundance_percent)) ==
      nrow(monthly), figure_width < figure_height && figure_width/figure_height < 0.4, all(file.exists(paste0(prefix,
      c(".png", ".tiff", ".svg", ".pdf"))))), stringsAsFactors = FALSE)

write.csv(qa, file.path(TABLE_DIR, "Cyanobacteria_vertical_curve_QA.csv"), row.names = FALSE, fileEncoding = "UTF-8")

if (!all(qa$passed)) stop("Figure QA failed; inspect the QA table.")

message("Completed cyanobacterial monthly relative-abundance curve: ", nrow(monthly), " observed months; ",
  nrow(raw), " samples; ", "4 sites x 2 size fractions per observed month.")
