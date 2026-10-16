rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

analysis_root <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/16_距离衰减")

input_dir <- file.path(analysis_root, "结果表格", "基础_注释无关", "整体_距离衰减", "outputs",
  "decay_analysis")

pair_file <- file.path(input_dir, "decay_pairwise_plot_data.csv")

mantel_file <- file.path(input_dir, "mantel_decay_results.csv")

figure_dir <- file.path(analysis_root, "结果图片", "基础_注释无关", "整体_距离衰减",
  "outputs", "decay_All_Large_Small_combined")

table_dir <- file.path(analysis_root, "结果表格", "基础_注释无关", "整体_距离衰减", "outputs",
  "decay_All_Large_Small_combined")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pair_file)) stop("Pairwise data not found: ", pair_file)

if (!file.exists(mantel_file)) stop("Mantel results not found: ", mantel_file)

group_levels <- c("All", "ConF", "USF")

group_labels <- c(all = "All", large = "ConF", small = "USF")

distance_levels <- c("Environmental distance", "Geographic distance", "Temporal distance")

distance_axis_labels <- c(`Environmental distance` = "Environmental distance\n(standardized Euclidean; Q excluded)",
  `Geographic distance` = "Geographic distance (km)", `Temporal distance` = "Temporal distance (months)")

pair_wide <- read.csv(pair_file, stringsAsFactors = FALSE, check.names = FALSE)

plot_data <- pair_wide %>%
  transmute(Group = factor(unname(group_labels[tolower(trimws(Group))]), levels = group_levels), Pair_id,
    Bray_Curtis_similarity = 1 - as.numeric(Bray_Curtis), Environmental_distance = as.numeric(Environmental_distance),
    Geographic_distance = as.numeric(Geographic_distance_km), Temporal_distance = as.numeric(Temporal_distance_month)) %>%
  pivot_longer(cols = c(Environmental_distance, Geographic_distance, Temporal_distance), names_to = "Distance_type",
    values_to = "Distance") %>%
  mutate(Distance_type = recode(Distance_type, Environmental_distance = "Environmental distance", Geographic_distance = "Geographic distance",
    Temporal_distance = "Temporal distance"), Distance_type = factor(Distance_type, levels = distance_levels)) %>%
  filter(!is.na(Group), is.finite(Distance), is.finite(Bray_Curtis_similarity))

mantel_data <- read.csv(mantel_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Group = factor(unname(group_labels[tolower(trimws(Group))]), levels = group_levels), Distance_type = factor(Distance_type,
    levels = distance_levels)) %>%
  filter(!is.na(Group), !is.na(Distance_type)) %>%
  transmute(Group, Distance_type, N_sample, N_taxa, Env_vars_used, Mantel_r_dissimilarity = Mantel_r,
    Mantel_r_similarity = -Mantel_r, Mantel_P = P_value, Mantel_P_BH = P_adj_BH)

env_description <- unique(mantel_data$Env_vars_used[mantel_data$Distance_type == "Environmental distance"])

if (any(grepl("(^|; )Q($|;)", env_description))) {
  stop("The input still contains Q in environmental distance. ", "Please rerun 整体_距离衰减.R after excluding Q.")
}

fit_statistics <- plot_data %>%
  group_by(Distance_type, Group) %>%
  group_modify(~{
    fit <- lm(Bray_Curtis_similarity ~ Distance, data = .x)
    sm <- summary(fit)
    tibble(N_pair = nrow(.x), Intercept = unname(coef(fit)[1]), Slope = unname(coef(fit)[2]), R_squared = sm$r.squared,
      LM_P = sm$coefficients["Distance", "Pr(>|t|)"])
  }) %>%
  ungroup() %>%
  left_join(mantel_data, by = c("Distance_type", "Group"))

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3))))
}

panel_x <- plot_data %>%
  group_by(Distance_type) %>%
  summarise(x = min(Distance, na.rm = TRUE) + 0.035 * diff(range(Distance, na.rm = TRUE)), .groups = "drop")

annotation_data <- fit_statistics %>%
  left_join(panel_x, by = "Distance_type") %>%
  mutate(label = paste0(as.character(Group), ": ", "R² = ", sprintf("%.3f", R_squared), "; Mantel r = ",
    sprintf("%.3f", Mantel_r_similarity), "; P ", format_p(Mantel_P)), y = case_when(Group == "All" ~
    0.985, Group == "ConF" ~ 0.94, TRUE ~ 0.895))

write.csv(plot_data, file.path(table_dir, "distance_decay_All_Large_Small_plot_data.csv"), row.names = FALSE)

write.csv(fit_statistics, file.path(table_dir, "distance_decay_All_Large_Small_statistics.csv"), row.names = FALSE)

group_colors <- c(All = "#0072B2", ConF = "#D55E00", USF = "#009E73")

plot_data <- plot_data %>%
  arrange(factor(Group, levels = c("All", "ConF", "USF")), Distance_type)

p <- ggplot(plot_data, aes(x = Distance, y = Bray_Curtis_similarity, colour = Group, fill = Group)) +
  geom_point(size = 1.25, alpha = 0.14, stroke = 0, show.legend = FALSE) + geom_smooth(method = "lm",
  formula = y ~ x, se = TRUE, linewidth = 1, alpha = 0.12, fullrange = TRUE) + geom_text(data = annotation_data,
  aes(x = x, y = y, label = label, colour = Group), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.8,
  family = "sans", show.legend = FALSE) + facet_wrap(~Distance_type, nrow = 1, scales = "free_x", labeller = as_labeller(distance_axis_labels)) +
  scale_colour_manual(values = group_colors, breaks = group_levels, drop = FALSE) + scale_fill_manual(values = group_colors,
  breaks = group_levels, drop = FALSE) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
  expand = expansion(mult = c(0.01, 0.015))) + scale_x_continuous(expand = expansion(mult = c(0.015,
  0.025))) + labs(title = "Distance-decay relationships of microbial community", x = NULL, y = "Bray-Curtis similarity",
  colour = "Group", fill = "Group") + guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.1),
  nrow = 1, byrow = TRUE), fill = "none") + theme_classic(base_size = 13.5, base_family = "sans") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65), axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.45), axis.text = element_text(colour = "black",
      size = 11.2), axis.title.y = element_text(colour = "black", size = 13), strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 12.5, face = "bold", margin = margin(b = 5)),
    plot.title = element_text(colour = "black", size = 14.5, face = "plain", hjust = 0, margin = margin(b = 5)),
    legend.position = "top", legend.justification = "left", legend.box.just = "left", legend.title = element_text(face = "bold",
      size = 12), legend.text = element_text(size = 11), legend.key.width = grid::unit(1.25, "cm"),
    legend.key.height = grid::unit(0.4, "cm"), legend.margin = margin(0, 0, 2, 0), panel.spacing.x = grid::unit(0.45,
      "cm"), plot.margin = margin(8, 10, 8, 8))

file_stub <- "distance_decay_All_Large_Small_3distances"

png_file <- file.path(figure_dir, paste0(file_stub, ".png"))

pdf_file <- file.path(figure_dir, paste0(file_stub, ".pdf"))

svg_file <- file.path(figure_dir, paste0(file_stub, ".svg"))

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(png_file, p, width = 12, height = 4.9, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
} else {
  ggsave(png_file, p, width = 12, height = 4.9, units = "in", dpi = 600, bg = "white")
}

ggsave(pdf_file, p, width = 12, height = 4.9, units = "in", device = cairo_pdf, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(svg_file, p, width = 12, height = 4.9, units = "in", device = svglite::svglite, bg = "white")
}

print(fit_statistics)

message("Figure written to: ", figure_dir)

message("Tables written to: ", table_dir)
