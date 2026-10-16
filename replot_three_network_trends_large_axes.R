.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

TABLE_DIR <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "replot_inputs")

PLOT_DIR <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "replot_outputs")

metrics <- read.csv(file.path(TABLE_DIR, "Sliding_Window_Topology_Metrics_ConF_USF.csv"), check.names = FALSE,
  stringsAsFactors = FALSE)

predictions <- read.csv(file.path(TABLE_DIR, "GAMM_Network_Predictions.csv"), check.names = FALSE, stringsAsFactors = FALSE)

statistics <- read.csv(file.path(TABLE_DIR, "GAMM_Network_Model_Statistics.csv"), check.names = FALSE,
  stringsAsFactors = FALSE)

metrics$Date <- as.Date(metrics$Date)

predictions$Date <- as.Date(predictions$Date)

metrics$size_group <- factor(metrics$size_group, levels = c("ConF", "USF"))

predictions$size_group <- factor(predictions$size_group, levels = c("ConF", "USF"))

selected_metrics <- c("Degree", "Robustness_weight", "Density")

metric_labels <- c(Degree = "Mean degree", Robustness_weight = "Natural connectivity (weighted)", Density = "Density")

group_colors <- c(ConF = "#D55E00", USF = "#009E73")

format_gamm_p <- function(p_value) {
  if (length(p_value) == 0L || is.na(p_value) || !is.finite(p_value)) {
    return("NA")
  }
  if (p_value < 0.001)
    return("<0.001")
  sprintf("=%.3f", p_value)
}

subtitle_for <- function(row) {
  paste0("GAMM smooth: ConF P", format_gamm_p(row$p_ConF_smooth), "; USF P", format_gamm_p(row$p_USF_smooth),
    "\nDifference smooth (USF-ConF): P", format_gamm_p(row$p_difference_USF_minus_ConF_smooth))
}

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
      8, 5, 10))
}

month_breaks <- sort(unique(metrics$Date))

month_breaks_display <- month_breaks[seq(1L, length(month_breaks), by = 2L)]

date_limits_plot <- range(month_breaks, na.rm = TRUE) + c(-7, 7)

for (var in selected_metrics) {
  plot_data <- metrics %>%
    filter(!is.na(.data[[var]]), !is.infinite(.data[[var]]))
  curve_data <- predictions %>%
    filter(metric == var)
  statistics_row <- statistics %>%
    filter(metric == var)
  p <- ggplot(plot_data, aes(x = Date, y = .data[[var]], color = size_group, fill = size_group)) +
    geom_ribbon(data = curve_data, aes(x = Date, ymin = lower, ymax = upper, fill = size_group),
      inherit.aes = FALSE, alpha = 0.18, color = NA) + geom_line(data = curve_data, aes(x = Date,
    y = fit, color = size_group), inherit.aes = FALSE, linewidth = 0.95) + geom_point(size = 1.65,
    alpha = 0.48, position = position_jitter(width = 5, height = 0, seed = 20260831)) + scale_color_manual(values = group_colors,
    breaks = names(group_colors), name = "Group") + scale_fill_manual(values = group_colors, breaks = names(group_colors),
    name = "Group") + scale_x_date(limits = date_limits_plot, breaks = month_breaks_display, date_labels = "%Y-%m",
    expand = expansion(mult = c(0.025, 0.025))) + labs(title = paste0("Network topology (", unname(metric_labels[[var]]),
    ")"), subtitle = subtitle_for(statistics_row), x = "Month", y = unname(metric_labels[[var]]),
    caption = paste0("Lines show GAMM fits; shaded ribbons indicate 95% confidence intervals.\n",
      "Three-month overlapping windows; step = 1 month; ASV prevalence threshold = 20%.")) + guides(color = guide_legend(nrow = 1,
    byrow = TRUE, override.aes = list(linewidth = 1, alpha = 1, size = 2)), fill = "none") + theme_cns_closed() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16.5, margin = margin(t = 4)),
      axis.title.y = element_text(size = 17.5, face = "bold", margin = margin(r = 10)), legend.position = "top",
      legend.direction = "horizontal", legend.justification = "left")
  save_path <- file.path(PLOT_DIR, paste0(var, "_Trend_ConF_USF.png"))
  ggsave(save_path, plot = p, width = 7.2, height = 5.9, dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", save_path), plot = p, width = 7.2, height = 5.9, device = cairo_pdf,
    bg = "white")
}

message("Re-rendered: ", paste(selected_metrics, collapse = ", "))
