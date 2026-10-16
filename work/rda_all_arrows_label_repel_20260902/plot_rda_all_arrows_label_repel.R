rm(list = ls())

.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/r_packages_extra"), file.path(Sys.getenv("DJK_PROJECT_ROOT",
  unset = "."), "work/figure2_hydro_core/Rlib"), .libPaths()))

required_packages <- c("ggplot2", "ggrepel", "patchwork")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()

project_root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

default_input_dir <- file.path(project_root, "01_分析项目", "14_RDA_CCA", "结果表格", "基础_注释无关",
  "整体_RDA_CCA_抽平去Q版", "outputs", "RDA_CCA_rarefied_noQ")

default_output_dir <- file.path(project_root, "01_分析项目", "14_RDA_CCA", "结果图片", "基础_注释无关",
  "整体_RDA_CCA_抽平去Q版", "outputs", "RDA_CCA_rarefied_noQ", "正方形三组RDA")

input_dir <- Sys.getenv("RDA_SQUARE_INPUT_DIR", unset = default_input_dir)

output_dir <- Sys.getenv("RDA_SQUARE_OUTPUT_DIR", unset = default_output_dir)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

GROUP_COLORS <- c(All = "#2C7FB8", small = "#2A9D8F", large = "#E67E22")

GROUP_LABELS <- c(All = "All", small = "USF", large = "ConF")

SET_ORDER <- c("All", "large", "small")

summary_table <- read.csv(file.path(input_dir, "RDA_CCA_summary_noQ.csv"), check.names = FALSE, stringsAsFactors = FALSE)

summary_table <- summary_table[summary_table$method == "RDA", , drop = FALSE]

make_hull <- function(df) {
  unique_xy <- unique(df[, c("Axis1", "Axis2")])
  if (nrow(unique_xy) < 3L)
    return(unique_xy)
  unique_xy[chull(unique_xy$Axis1, unique_xy$Axis2), , drop = FALSE]
}

make_rda_plot <- function(set_name, show_legend = TRUE) {
  site <- read.csv(file.path(input_dir, paste0(set_name, "_RDA_sample_scores.csv")), check.names = FALSE,
    stringsAsFactors = FALSE)
  env <- read.csv(file.path(input_dir, paste0(set_name, "_RDA_env_scores.csv")), check.names = FALSE,
    stringsAsFactors = FALSE)
  stats <- summary_table[summary_table$set == set_name, , drop = FALSE]
  if (nrow(stats) != 1L)
    stop("Missing RDA summary for set: ", set_name)
  site$plot_group <- factor(set_name, levels = names(GROUP_COLORS))
  hull <- make_hull(site)
  hull$plot_group <- factor(set_name, levels = names(GROUP_COLORS))
  x_range <- range(c(site$Axis1, env$Axis1), na.rm = TRUE)
  y_range <- range(c(site$Axis2, env$Axis2), na.rm = TRUE)
  x_pad <- diff(x_range) * 0.2
  y_pad <- diff(y_range) * 0.2
  x_limits <- x_range + c(-x_pad, x_pad)
  y_limits <- y_range + c(-y_pad, y_pad)
  env$arrow_length <- sqrt(env$Axis1^2 + env$Axis2^2)
  env$highlight <- FALSE
  eligible <- which(is.finite(env$p_adj_BH) & env$p_adj_BH < 0.05)
  if (length(eligible) > 0L) {
    highlight_order <- eligible[order(env$p_adj_BH[eligible], -env$arrow_length[eligible], env$variable[eligible])]
    env$highlight[head(highlight_order, 5L)] <- TRUE
  }
  shaft_fraction <- seq(0, 1.04, length.out = 50L)
  arrow_obstacles <- do.call(rbind, lapply(seq_len(nrow(env)), function(i) {
    data.frame(repel_x = env$Axis1[i] * shaft_fraction, repel_y = env$Axis2[i] * shaft_fraction,
      repel_label = "", stringsAsFactors = FALSE)
  }))
  label_points <- data.frame(repel_x = env$Axis1, repel_y = env$Axis2, repel_label = env$label, stringsAsFactors = FALSE)
  repel_data <- rbind(arrow_obstacles, label_points)
  stat_x <- x_limits[1] + 0.035 * diff(x_limits)
  stat_y <- y_limits[2] - 0.035 * diff(y_limits)
  stat_label <- sprintf("R^2 == %.3f*','~~P == %.3f", stats$R2, stats$p_value)
  legend_df <- data.frame(Axis1 = NA_real_, Axis2 = NA_real_, plot_group = factor(SET_ORDER, levels = names(GROUP_COLORS)))
  p <- ggplot(site, aes(Axis1, Axis2)) + geom_hline(yintercept = 0, linewidth = 0.3, color = "grey76") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey76") + geom_polygon(data = hull, aes(fill = plot_group,
    group = plot_group), alpha = 0.14, color = NA, inherit.aes = TRUE) + geom_point(aes(fill = plot_group,
    color = plot_group), shape = 21, size = 3.2, stroke = 0.55, alpha = 0.88) + geom_point(data = legend_df,
    aes(Axis1, Axis2, fill = plot_group, color = plot_group), inherit.aes = FALSE, shape = 21, size = 3.2,
    alpha = 0, na.rm = TRUE, show.legend = show_legend) + geom_segment(data = env[!env$highlight,
    , drop = FALSE], aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, linewidth = 0.55,
    color = "black", arrow = grid::arrow(length = grid::unit(0.17, "cm"))) + geom_segment(data = env[env$highlight,
    , drop = FALSE], aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, linewidth = 1.05,
    color = "#E6A700", arrow = grid::arrow(length = grid::unit(0.19, "cm"))) + ggrepel::geom_text_repel(data = repel_data,
    aes(repel_x, repel_y, label = repel_label), inherit.aes = FALSE, seed = 123, size = 3, color = "black",
    box.padding = grid::unit(0.34, "lines"), point.padding = grid::unit(0.18, "lines"), min.segment.length = 0,
    segment.color = "grey42", segment.size = 0.25, segment.alpha = 0.78, force = 2.2, force_pull = 0.58,
    max.overlaps = Inf, max.time = 18, max.iter = 2e+05, direction = "both", xlim = x_limits, ylim = y_limits) +
    annotate("text", x = stat_x, y = stat_y, label = stat_label, hjust = 0, vjust = 1, size = 4.1,
      fontface = "plain", parse = TRUE) + scale_fill_manual(values = GROUP_COLORS, breaks = SET_ORDER,
    labels = GROUP_LABELS[SET_ORDER], drop = FALSE) + scale_color_manual(values = GROUP_COLORS, breaks = SET_ORDER,
    labels = GROUP_LABELS[SET_ORDER], drop = FALSE) + scale_x_continuous(limits = x_limits, expand = c(0,
    0)) + scale_y_continuous(limits = y_limits, expand = c(0, 0)) + labs(x = sprintf("RDA1 (%.1f%%)",
    stats$axis1_percent), y = sprintf("RDA2 (%.1f%%)", stats$axis2_percent), fill = "Group", color = "Group") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1)), color = "none") +
    theme_classic(base_size = 11) + theme(aspect.ratio = 1, axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"), legend.position = if (show_legend)
      "bottom"
    else "none", legend.direction = "horizontal", legend.title = element_text(face = "bold"), panel.border = element_rect(color = "black",
      fill = NA, linewidth = 0.65), plot.margin = margin(7, 7, 7, 7))
  p
}

plots_no_legend <- lapply(SET_ORDER, make_rda_plot, show_legend = FALSE)

names(plots_no_legend) <- SET_ORDER

for (set_name in SET_ORDER) {
  ggsave(file.path(output_dir, paste0(set_name, "_RDA_square.png")), plots_no_legend[[set_name]], width = 6,
    height = 6, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, paste0(set_name, "_RDA_square.pdf")), plots_no_legend[[set_name]], width = 6,
    height = 6, bg = "white")
}

plots_combined <- list(make_rda_plot("All", show_legend = TRUE), make_rda_plot("large", show_legend = FALSE),
  make_rda_plot("small", show_legend = FALSE))

combined <- wrap_plots(plots_combined, nrow = 1, guides = "collect")

ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.png"), combined, width = 15.5, height = 5.7,
  dpi = 600, bg = "white")

ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.pdf"), combined, width = 15.5, height = 5.7,
  bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.svg"), combined, width = 15.5,
    height = 5.7, bg = "white")
}

message("Square RDA panels saved to: ", output_dir)
