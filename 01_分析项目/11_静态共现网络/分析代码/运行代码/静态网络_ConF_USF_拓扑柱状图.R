rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

project_root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

input_file <- file.path(project_root, "01_分析项目/11_静态共现网络/结果表格/SILVA", "整体_静态网络_出现率20pct_增强版/Network_Global_Metrics.csv")

output_dir <- file.path(project_root, "01_分析项目/11_静态共现网络/结果图片/SILVA", "整体_静态网络_出现率20pct_增强版")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

network_metrics <- read.csv(input_file, check.names = FALSE)

required_groups <- c("ConF", "USF")

if (!all(required_groups %in% network_metrics$Group)) {
  stop("Network_Global_Metrics.csv must contain ConF and USF.")
}

metric_spec <- data.frame(Column = c("Network_nodes_nonisolated", "Edges", "Average_degree", "Density_network_nodes"),
  Label = c("Number of\nnodes", "Number of\nedges", "Average\ndegree", "Network\ndensity"), Digits = c(0,
    0, 1, 3), stringsAsFactors = FALSE)

plot_data <- do.call(rbind, lapply(seq_len(nrow(metric_spec)), function(i) {
  values <- network_metrics[match(required_groups, network_metrics$Group), metric_spec$Column[i]]
  data.frame(Group = required_groups, Metric = metric_spec$Label[i], Value = as.numeric(values), Digits = metric_spec$Digits[i],
    stringsAsFactors = FALSE)
}))

plot_data$Group <- factor(plot_data$Group, levels = required_groups)

plot_data$Metric <- factor(plot_data$Metric, levels = metric_spec$Label)

plot_data$ValueLabel <- mapply(function(value, digits) {
  formatC(value, format = "f", digits = digits, big.mark = ",")
}, plot_data$Value, plot_data$Digits, USE.NAMES = FALSE)

group_colors <- c(ConF = "#D55E00", USF = "#009E73")

p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) + geom_col(width = 0.4, colour = "white",
  linewidth = 0.35) + geom_text(aes(label = ValueLabel), vjust = -0.42, size = 4, family = "sans",
  colour = "#202124") + facet_wrap(~Metric, nrow = 2, ncol = 2, scales = "free_y") + scale_fill_manual(values = group_colors,
  breaks = required_groups, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = NULL, title = "Static network topology", caption = paste0("Bars show descriptive values from one static network per group;\n",
    "no biological-replicate significance test is implied.")) + guides(fill = guide_legend(nrow = 1,
  byrow = TRUE, override.aes = list(colour = NA))) + theme_classic(base_family = "sans", base_size = 12.5) +
  theme(panel.background = element_rect(fill = "#FFF5F0", colour = NA), panel.border = element_rect(fill = NA,
    colour = "#E8CFC4", linewidth = 0.55), panel.spacing.x = grid::unit(3.2, "mm"), panel.spacing.y = grid::unit(4,
    "mm"), strip.background = element_blank(), strip.text = element_text(colour = "#202124", face = "bold",
    size = 10.8, lineheight = 0.95, margin = margin(b = 4)), axis.line = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x = element_blank(), axis.ticks.y = element_line(colour = "#444444", linewidth = 0.35),
    axis.text.y = element_text(colour = "#333333", size = 10.2), legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 10.5), legend.key.width = grid::unit(5.5, "mm"), legend.spacing.x = grid::unit(2,
      "mm"), plot.title = element_text(colour = "#202124", face = "bold", size = 13, hjust = 0),
    plot.caption = element_text(colour = "#444444", size = 9, hjust = 0, margin = margin(t = 5)),
    plot.margin = margin(7, 8, 5, 7))

png_file <- file.path(output_dir, "Static_Network_Topology_Barplot_ConF_USF.png")

pdf_file <- file.path(output_dir, "Static_Network_Topology_Barplot_ConF_USF.pdf")

svg_file <- file.path(output_dir, "Static_Network_Topology_Barplot_ConF_USF.svg")

ggsave(png_file, plot = p, width = 3.8, height = 5, dpi = 600, bg = "white")

ggsave(pdf_file, plot = p, width = 3.8, height = 5, device = grDevices::cairo_pdf, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(svg_file, plot = p, width = 3.8, height = 5, device = svglite::svglite, bg = "white")
} else {
  ggsave(svg_file, plot = p, width = 3.8, height = 5, device = "svg", bg = "white")
}

message("Saved: ", png_file)

message("Saved: ", pdf_file)

message("Saved: ", svg_file)
