rm(list = ls())

set.seed(123)

args <- commandArgs(trailingOnly = TRUE)

script_args <- commandArgs(trailingOnly = FALSE)

file_arg <- grep("^--file=", script_args, value = TRUE)

script_dir <- Sys.getenv("DJK16S_ACTIVE_ROOT", unset = file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA/decay"))

input_dir <- if (length(args) >= 1L) args[[1]] else script_dir

output_dir <- if (length(args) >= 2L) {
  args[[2]]
} else {
  file.path(script_dir, "outputs", "water_level_distance_decay")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 9999L

MANTEL_METHOD <- "spearman"

MIN_TOTAL_ABUNDANCE <- 0.1

required_packages <- c("dplyr", "ggplot2", "vegan", "scales", "ragg", "svglite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(vegan)
})

read_csv_smart <- function(path, ...) {
  value <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM",
    ...), error = function(e) NULL)
  if (is.null(value)) {
    value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "GBK",
      ...)
  }
  value
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p))))
}

save_plot <- function(plot_object, prefix, width, height) {
  ggsave(file.path(output_dir, paste0(prefix, ".png")), plot_object, width = width, height = height,
    dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, paste0(prefix, ".pdf")), plot_object, width = width, height = height,
    bg = "white", device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(prefix, ".svg")), plot_object, width = width, height = height,
    bg = "white", device = svglite::svglite)
}

tpm_path <- file.path(input_dir, "tpm.csv")

env_path <- file.path(input_dir, "env.csv")

size_path <- file.path(input_dir, "groupsize.csv")

for (path in c(tpm_path, env_path, size_path)) {
  if (!file.exists(path))
    stop("Input file not found: ", path)
}

tpm_raw <- read_csv_smart(tpm_path)

env <- read_csv_smart(env_path)

group_size <- read_csv_smart(size_path)

names(group_size)[1:2] <- c("sample", "size_group")

env$sample <- trimws(as.character(env$sample))

group_size$sample <- trimws(as.character(group_size$sample))

group_size$size_group <- tolower(trimws(as.character(group_size$size_group)))

if (!"WL" %in% names(env)) stop("env.csv does not contain a WL column.")

env$WL <- suppressWarnings(as.numeric(env$WL))

taxon_column <- names(tpm_raw)[1]

taxon_ids <- make.unique(as.character(tpm_raw[[taxon_column]]))

tpm_raw[[taxon_column]] <- NULL

tpm_matrix <- as.matrix(tpm_raw)

mode(tpm_matrix) <- "numeric"

rownames(tpm_matrix) <- taxon_ids

tpm_matrix <- tpm_matrix[rowSums(tpm_matrix, na.rm = TRUE) > MIN_TOTAL_ABUNDANCE, , drop = FALSE]

abundance <- as.data.frame(t(tpm_matrix), check.names = FALSE)

abundance$sample <- rownames(abundance)

metadata <- env %>%
  inner_join(group_size, by = "sample") %>%
  filter(size_group %in% c("large", "small"), !is.na(WL))

common_samples <- intersect(metadata$sample, abundance$sample)

if (length(common_samples) < 3L) stop("Fewer than three complete samples remain.")

metadata <- metadata %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

abundance <- abundance %>%
  filter(sample %in% metadata$sample) %>%
  arrange(match(sample, metadata$sample))

rownames(abundance) <- abundance$sample

abundance$sample <- NULL

group_order <- c("All", "ConF", "USF")

analyse_group <- function(group_label) {
  if (group_label == "All") {
    sub_metadata <- metadata
  }
  else {
    raw_group <- c(ConF = "large", USF = "small")[[group_label]]
    sub_metadata <- metadata %>%
      filter(size_group == raw_group)
  }
  sub_abundance <- abundance[sub_metadata$sample, , drop = FALSE]
  sub_abundance <- sub_abundance[, colSums(sub_abundance, na.rm = TRUE) > 0, drop = FALSE]
  community_distance <- vegan::vegdist(sub_abundance, method = "bray")
  water_level_distance <- as.dist(abs(outer(sub_metadata$WL, sub_metadata$WL, "-")))
  mantel_result <- vegan::mantel(community_distance, water_level_distance, method = MANTEL_METHOD,
    permutations = N_PERM, na.rm = TRUE)
  pair_index <- which(lower.tri(matrix(FALSE, nrow = nrow(sub_metadata), ncol = nrow(sub_metadata))),
    arr.ind = TRUE)
  pair_table <- data.frame(Group = group_label, Sample_1 = sub_metadata$sample[pair_index[, 1]], Sample_2 = sub_metadata$sample[pair_index[,
    2]], Water_level_1_m = sub_metadata$WL[pair_index[, 1]], Water_level_2_m = sub_metadata$WL[pair_index[,
    2]], Water_level_distance_m = as.vector(water_level_distance), Bray_Curtis_dissimilarity = as.vector(community_distance),
    stringsAsFactors = FALSE) %>%
    mutate(Bray_Curtis_similarity = 1 - Bray_Curtis_dissimilarity)
  linear_model <- lm(Bray_Curtis_similarity ~ Water_level_distance_m, data = pair_table)
  model_summary <- summary(linear_model)
  slope_row <- coef(model_summary)["Water_level_distance_m", ]
  statistics <- data.frame(Group = group_label, N_samples = nrow(sub_metadata), N_taxa = ncol(sub_abundance),
    N_pairs = nrow(pair_table), Water_level_min_m = min(sub_metadata$WL), Water_level_max_m = max(sub_metadata$WL),
    Mantel_method = MANTEL_METHOD, Mantel_permutations = N_PERM, Mantel_r_dissimilarity = unname(mantel_result$statistic),
    Mantel_r_similarity = -unname(mantel_result$statistic), Mantel_P = mantel_result$signif, LM_slope_per_m = unname(slope_row["Estimate"]),
    LM_slope_SE = unname(slope_row["Std. Error"]), LM_slope_P_descriptive = unname(slope_row["Pr(>|t|)"]),
    LM_R2_descriptive = unname(model_summary$r.squared), stringsAsFactors = FALSE)
  list(pair = pair_table, statistics = statistics)
}

results <- lapply(group_order, analyse_group)

pair_data <- bind_rows(lapply(results, `[[`, "pair")) %>%
  mutate(Group = factor(Group, levels = group_order))

statistics <- bind_rows(lapply(results, `[[`, "statistics")) %>%
  mutate(Group = factor(Group, levels = group_order), Mantel_P_BH = p.adjust(Mantel_P, method = "BH")) %>%
  arrange(Group)

sample_audit <- metadata %>%
  transmute(Sample = sample, Size = size_group, Water_level_m = WL, Library_total = rowSums(abundance[metadata$sample,
    , drop = FALSE], na.rm = TRUE))

write.csv(pair_data, file.path(output_dir, "water_level_decay_pairwise_data.csv"), row.names = FALSE)

write.csv(statistics, file.path(output_dir, "water_level_decay_statistics.csv"), row.names = FALSE)

write.csv(sample_audit, file.path(output_dir, "water_level_decay_sample_audit.csv"), row.names = FALSE)

palette <- c(All = "#2C7FB8", ConF = "#E67E22", USF = "#2A9D8F")

annotation <- statistics %>%
  mutate(x = Inf, y = c(0.98, 0.925, 0.87), label = paste0(as.character(Group), ": Mantel r = ", sprintf("%.3f",
    Mantel_r_similarity), ", BH P ", format_p(Mantel_P_BH), "; R2 = ", sprintf("%.3f", LM_R2_descriptive)))

base_theme <- theme_classic(base_size = 11) + theme(axis.line = element_blank(), panel.border = element_rect(color = "black",
  fill = NA, linewidth = 0.65), panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
  panel.grid.minor = element_blank(), axis.text = element_text(color = "black"), axis.title = element_text(color = "black",
    face = "bold"), legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(face = "bold"),
  plot.title = element_text(face = "bold", size = 13, hjust = 0), plot.subtitle = element_text(color = "#4D4D4D",
    size = 10), plot.margin = margin(8, 10, 8, 8))

p_overlay <- ggplot(pair_data, aes(x = Water_level_distance_m, y = Bray_Curtis_similarity, color = Group)) +
  geom_point(size = 0.72, alpha = 0.1, show.legend = TRUE) + geom_smooth(aes(fill = Group), method = "lm",
  formula = y ~ x, se = TRUE, linewidth = 1.05, alpha = 0.13, show.legend = FALSE) + geom_text(data = annotation,
  aes(x = x, y = y, label = label, color = Group), inherit.aes = FALSE, hjust = 1.03, vjust = 1, size = 3.25,
  fontface = "bold", show.legend = FALSE) + scale_color_manual(values = palette, breaks = group_order,
  drop = FALSE) + scale_fill_manual(values = palette, breaks = group_order, drop = FALSE) + scale_y_continuous(limits = c(0,
  1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.01))) + scale_x_continuous(expand = expansion(mult = c(0.01,
  0.04))) + labs(x = "Water-level distance |WL1 - WL2| (m)", y = "Bray-Curtis similarity", color = "Group",
  fill = "Group", title = "Water-level distance decay", subtitle = "Lines are descriptive linear fits; Mantel tests use 9,999 permutations") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.2), nrow = 1), fill = "none") +
  base_theme

save_plot(p_overlay, "water_level_distance_decay_All_Large_Small", 6.6, 5.5)

p_faceted <- ggplot(pair_data, aes(x = Water_level_distance_m, y = Bray_Curtis_similarity, color = Group)) +
  geom_point(size = 0.68, alpha = 0.16) + geom_smooth(aes(fill = Group), method = "lm", formula = y ~
  x, se = TRUE, linewidth = 0.95, alpha = 0.15, show.legend = FALSE) + facet_wrap(~Group, nrow = 1) +
  scale_color_manual(values = palette, drop = FALSE) + scale_fill_manual(values = palette, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + labs(x = "Water-level distance |WL1 - WL2| (m)",
  y = "Bray-Curtis similarity", title = "Water-level distance decay", subtitle = "Group-specific views of the same pairwise analysis") +
  guides(color = "none", fill = "none") + base_theme + theme(strip.background = element_rect(fill = "white",
  color = "black", linewidth = 0.55), strip.text = element_text(face = "bold", color = "black"))

save_plot(p_faceted, "water_level_distance_decay_faceted", 9.6, 4.2)

capture.output({
  cat("Water-level distance-decay analysis\n")
  cat("===================================\n")
  cat("Community input: ", normalizePath(tpm_path, winslash = "/"), "\n", sep = "")
  cat("Environmental input: ", normalizePath(env_path, winslash = "/"), "\n", sep = "")
  cat("Size input: ", normalizePath(size_path, winslash = "/"), "\n", sep = "")
  cat("Taxon filter: total abundance > ", MIN_TOTAL_ABUNDANCE, "\n", sep = "")
  cat("Mantel: ", MANTEL_METHOD, ", ", N_PERM, " permutations\n\n", sep = "")
  print(statistics)
  cat("\nLibrary total range: ", min(sample_audit$Library_total), " to ", max(sample_audit$Library_total),
    "\n", sep = "")
}, file = file.path(output_dir, "water_level_decay_analysis_log.txt"))

message("Completed. Results written to: ", normalizePath(output_dir, winslash = "/"))
