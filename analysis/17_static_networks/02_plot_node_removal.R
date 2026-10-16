rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

table_dir <- file.path(project_root, "01_分析项目/11_静态共现网络/结果表格/SILVA", "整体_静态网络_出现率20pct_增强版")

output_dir <- file.path(project_root, "01_分析项目/11_静态共现网络/结果图片/SILVA", "整体_静态网络_出现率20pct_增强版")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_data <- read.csv(file.path(table_dir, "Node_Removal_Stability_Summary.csv"), check.names = FALSE)

raw_data <- read.csv(file.path(table_dir, "Node_Removal_Stability_Raw.csv"), check.names = FALSE)

robustness_50 <- read.csv(file.path(table_dir, "Node_Removal_50pct_Robustness.csv"), check.names = FALSE)

comparison_50 <- read.csv(file.path(table_dir, "Node_Removal_50pct_Group_Comparison.csv"), check.names = FALSE)

group_levels <- c("ConF", "USF")

group_colors <- c(ConF = "#D55E00", USF = "#009E73")

group_shapes <- c(ConF = 16, USF = 17)

single_trajectory <- raw_data[raw_data$Replicate == 1, , drop = FALSE]

summary_data$Group <- factor(summary_data$Group, levels = group_levels)

raw_data$Group <- factor(raw_data$Group, levels = group_levels)

robustness_50$Group <- factor(robustness_50$Group, levels = group_levels)

theme_network_panel <- function() {
  theme_classic(base_family = "sans", base_size = 11.5) + theme(panel.background = element_rect(fill = "#F1F5F8",
    colour = NA), panel.border = element_rect(fill = NA, colour = "#B8C5D0", linewidth = 0.55), axis.line = element_blank(),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.35), axis.text = element_text(colour = "#333333",
      size = 9.5), axis.title = element_text(colour = "#202124", size = 10.5), plot.title = element_text(colour = "#202124",
      face = "bold", size = 10.8, hjust = 0.5, lineheight = 0.95, margin = margin(b = 4)), legend.title = element_blank(),
    legend.text = element_text(size = 10), legend.key.width = grid::unit(5.5, "mm"), plot.margin = margin(4,
      4, 3, 4))
}

p_degree <- ggplot(single_trajectory, aes(x = RemovalPercent, y = AverageDegree, colour = Group, shape = Group,
  group = Group)) + geom_line(linewidth = 0.82) + geom_point(size = 1.05, stroke = 0.22) + scale_colour_manual(values = group_colors,
  breaks = group_levels) + scale_shape_manual(values = group_shapes, breaks = group_levels) + scale_x_continuous(limits = c(0,
  100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.02, 0.02))) + scale_y_continuous(expand = expansion(mult = c(0.02,
  0.08))) + labs(title = "Average degree\n(one fixed trajectory)", x = "Node removal (%)", y = "Average degree",
  colour = NULL, shape = NULL) + guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.9,
  size = 2.1, shape = unname(group_shapes[group_levels]))), shape = "none") + theme_network_panel()

q_value <- comparison_50$BH_q_value[comparison_50$Metric == "NaturalConnectivityRetention"]

if (length(q_value) != 1L || !is.finite(q_value)) {
  stop("NaturalConnectivityRetention BH q-value was not found.")
}

box_max <- max(robustness_50$NaturalConnectivityRetention, na.rm = TRUE)

bracket_y <- box_max * 1.07

tick_height <- box_max * 0.035

annotation_y <- box_max * 1.145

p_robustness <- ggplot(robustness_50, aes(x = Group, y = NaturalConnectivityRetention, colour = Group,
  fill = Group)) + gghalves::geom_half_violin(side = "l", trim = FALSE, width = 0.92, alpha = 0.34,
  colour = NA) + geom_boxplot(width = 0.19, alpha = 0.58, outlier.shape = NA, linewidth = 0.62, colour = "black") +
  geom_jitter(width = 0.105, height = 0, size = 1.28, alpha = 0.72, shape = 21, stroke = 0.32, colour = "white") +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.45, colour = "#333333") +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - tick_height, linewidth = 0.45,
    colour = "#333333") + annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y -
  tick_height, linewidth = 0.45, colour = "#333333") + annotate("text", x = 1.5, y = annotation_y,
  label = paste0("Simulation BH q = ", format(q_value, scientific = TRUE, digits = 2)), size = 3.5,
  family = "sans", colour = "#202124") + scale_colour_manual(values = group_colors, breaks = group_levels) +
  scale_fill_manual(values = group_colors, breaks = group_levels) + scale_y_continuous(limits = c(0,
  box_max * 1.24), expand = expansion(mult = c(0, 0.01))) + labs(title = "Robustness\n(50% removal)",
  x = NULL, y = "Natural connectivity retention", colour = NULL, fill = NULL) + guides(colour = "none",
  fill = "none") + theme_classic(base_family = "sans", base_size = 11.8) + theme(panel.border = element_rect(fill = NA,
  colour = "black", linewidth = 0.75), axis.line = element_blank(), axis.text = element_text(colour = "#202124",
  size = 10.2), axis.text.x = element_text(face = "bold", size = 11.5), axis.title.y = element_text(colour = "#202124",
  size = 10.8), plot.title = element_text(colour = "#202124", face = "bold", size = 11.2, hjust = 0.5,
  lineheight = 0.95, margin = margin(b = 4)), plot.margin = margin(4, 5, 3, 5))

combined_plot <- (p_degree | p_robustness) + plot_layout(widths = c(1.32, 0.88), guides = "collect") +
  plot_annotation(title = "Static network stability under random node removal", caption = paste0("Left: one fixed seeded random-removal trajectory, shown to retain stepwise fluctuations.\n",
    "Right: 100 nested random trajectories at 50% removal; BH q describes simulations, not biological replicates."),
    theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 11.3, colour = "#202124",
      hjust = 0), plot.caption = element_text(family = "sans", size = 8.8, colour = "#444444",
      hjust = 0, margin = margin(t = 5)), plot.margin = margin(6, 7, 5, 7))) & theme(legend.position = "bottom")

png_file <- file.path(output_dir, "Static_Network_Stability_Composite_ConF_USF.png")

pdf_file <- file.path(output_dir, "Static_Network_Stability_Composite_ConF_USF.pdf")

svg_file <- file.path(output_dir, "Static_Network_Stability_Composite_ConF_USF.svg")

ggsave(png_file, plot = combined_plot, width = 5.35, height = 4.75, dpi = 600, bg = "white")

ggsave(pdf_file, plot = combined_plot, width = 5.35, height = 4.75, device = grDevices::cairo_pdf, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(svg_file, plot = combined_plot, width = 5.35, height = 4.75, device = svglite::svglite, bg = "white")
} else {
  ggsave(svg_file, plot = combined_plot, width = 5.35, height = 4.75, device = "svg", bg = "white")
}

message("Saved: ", png_file)

message("Saved: ", pdf_file)

message("Saved: ", svg_file)
