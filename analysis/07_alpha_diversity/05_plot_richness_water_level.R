required_packages <- c("ggplot2", "dplyr", "patchwork")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(20260809)

required_files <- c("alpha1.csv", "groupmonth.csv", "groupsize.csv", "GAMM_Alpha_Model_Statistics.csv",
  "GAMM_Alpha_Predictions.csv", "WaterLevel_Monthly_Data.csv", "GAMM_WaterLevel_Model_Statistics.csv",
  "GAMM_WaterLevel_Predictions.csv")

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("缺少以下文件，请先运行整体 Alpha 和水位 GAMM 脚本：\n", paste(missing_files,
    collapse = "\n"))
}

alpha_data <- read.csv("alpha1.csv", stringsAsFactors = FALSE, check.names = FALSE)

names(alpha_data)[1] <- "sample"

group_month <- read.csv("groupmonth.csv", stringsAsFactors = FALSE, check.names = FALSE)

group_size <- read.csv("groupsize.csv", stringsAsFactors = FALSE, check.names = FALSE)

names(group_month)[1:2] <- c("sample", "month_group")

names(group_size)[1:2] <- c("sample", "size_group")

alpha_observed <- alpha_data %>%
  inner_join(group_month, by = "sample") %>%
  inner_join(group_size, by = "sample") %>%
  mutate(Date = as.Date(paste0("20", month_group, "01"), format = "%Y%m%d"), size_group = tolower(trimws(as.character(size_group)))) %>%
  filter(!is.na(Date), is.finite(Richness), size_group %in% c("large", "small"))

if (nrow(alpha_observed) == 0L) {
  stop("没有可用于绘制 Richness 趋势的有效数据。")
}

alpha_all <- alpha_observed

alpha_all$size_group <- "All"

alpha_plot_data <- bind_rows(alpha_observed, alpha_all) %>%
  mutate(size_group = factor(size_group, levels = c("large", "small", "All")))

alpha_predictions <- read.csv("GAMM_Alpha_Predictions.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(metric == "Richness") %>%
  mutate(Date = as.Date(Date), size_group = factor(size_group, levels = c("large", "small", "All")))

alpha_statistics <- read.csv("GAMM_Alpha_Model_Statistics.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(metric == "Richness")

water_observed <- read.csv("WaterLevel_Monthly_Data.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Date = as.Date(Date))

water_predictions <- read.csv("GAMM_WaterLevel_Predictions.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Date = as.Date(Date))

water_statistics <- read.csv("GAMM_WaterLevel_Model_Statistics.csv", stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(alpha_statistics) != 1L) {
  stop("Richness 的 GAMM 统计结果不是唯一一行。")
}

if (nrow(water_statistics) != 1L) {
  stop("水位 GAMM 统计结果不是唯一一行。")
}

format_p <- function(p_value) {
  if (length(p_value) == 0L || !is.finite(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  sprintf("=%.3f", p_value)
}

alpha_subtitle <- paste0("GAMM smooth: ConF P", format_p(alpha_statistics$p_large_smooth), "; USF P",
  format_p(alpha_statistics$p_small_smooth), "; All P", format_p(alpha_statistics$p_overall_smooth),
  "\nDifference smooth (USF-ConF): P", format_p(alpha_statistics$p_difference_smooth))

water_subtitle <- paste0("GAMM smooth: P", format_p(water_statistics$smooth_p), "; k=", water_statistics$k,
  "; REML; continuous-time CAR(1); n=", water_statistics$n_months, " months")

month_breaks <- sort(unique(c(alpha_observed$Date, water_observed$Date)))

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
      margin = margin(b = 5)), plot.tag = element_text(colour = "black", size = 12, face = "bold"),
    plot.tag.position = "topleft", legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 8.2),
    legend.key.width = grid::unit(10, "mm"), legend.key.height = grid::unit(3.5, "mm"), plot.margin = margin(7,
      8, 5, 7))
}

group_colors <- c(large = "#D55E00", small = "#009E73", All = "#0072B2")

group_labels <- c(large = "ConF", small = "USF", All = "All")

p_richness_cns <- ggplot(alpha_plot_data, aes(x = Date, y = Richness, colour = size_group)) + geom_ribbon(data = alpha_predictions,
  aes(x = Date, ymin = lower, ymax = upper, fill = size_group), inherit.aes = FALSE, alpha = 0.18,
  colour = NA, show.legend = FALSE) + geom_line(data = alpha_predictions, aes(x = Date, y = fit, colour = size_group),
  inherit.aes = FALSE, linewidth = 0.95) + geom_point(size = 1.65, alpha = 0.48, position = position_jitter(width = 5,
  height = 0, seed = 20260809)) + scale_colour_manual(values = group_colors, breaks = names(group_colors),
  labels = group_labels, name = "Group") + scale_fill_manual(values = group_colors, breaks = names(group_colors),
  labels = group_labels, name = "Group") + scale_x_date(limits = date_limits_plot, breaks = main_month_breaks,
  date_labels = "%Y-%m", expand = expansion(mult = c(0.025, 0.025))) + labs(title = "Alpha diversity (Richness)",
  subtitle = alpha_subtitle, x = "Month", y = "Richness") + guides(colour = guide_legend(nrow = 1,
  byrow = TRUE, override.aes = list(linewidth = 1, alpha = 1)), fill = "none") + theme_cns_closed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11.2), legend.position = "top",
    legend.direction = "horizontal", legend.justification = "left")

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
    size = 5.4, lineheight = 1, margin = margin(b = 2)), plot.tag = element_text(colour = "black",
    size = 8, face = "bold"), plot.margin = margin(4, 4, 4, 4))

p_water_transparent <- p_water_inset + theme(panel.background = element_rect(fill = NA, colour = NA),
  plot.background = element_rect(fill = NA, colour = NA))

combined_plot <- p_richness_cns + inset_element(p_water_inset, left = 0.045, bottom = 0.52, right = 0.49,
  top = 0.965, align_to = "panel", on_top = TRUE, clip = TRUE) + plot_annotation(caption = paste0("Lines show GAMM fits; shaded ribbons indicate 95% confidence intervals. ",
  "Richness observations: n=120; water-level observations: n=15 months."), theme = theme(plot.caption = element_text(family = "sans",
  size = 7.5, colour = "#333333", hjust = 0, margin = margin(t = 4)), plot.background = element_rect(fill = "white",
  colour = NA)))

print(combined_plot)

ggsave("Richness_WaterLevel_Combined_CNS.png", plot = combined_plot, width = 7.2, height = 5.9, dpi = 600,
  bg = "white")

ggsave("Richness_WaterLevel_Combined_CNS.pdf", plot = combined_plot, width = 7.2, height = 5.9, device = cairo_pdf,
  bg = "white")

ggsave("WaterLevel_Inset_Transparent.png", plot = p_water_transparent, width = 3.6, height = 2.6, dpi = 600,
  bg = "transparent")

ggsave("WaterLevel_Inset_Transparent.pdf", plot = p_water_transparent, width = 3.6, height = 2.6, device = cairo_pdf,
  bg = "transparent")

cat("已输出：\n", "1. Richness_WaterLevel_Combined_CNS.png\n", "2. Richness_WaterLevel_Combined_CNS.pdf\n",
  "3. WaterLevel_Inset_Transparent.png\n", "4. WaterLevel_Inset_Transparent.pdf\n")
