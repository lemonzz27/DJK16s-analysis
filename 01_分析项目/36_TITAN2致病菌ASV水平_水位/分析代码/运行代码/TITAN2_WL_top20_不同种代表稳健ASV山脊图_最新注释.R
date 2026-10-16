options(stringsAsFactors = FALSE, scipen = 999)

try(Sys.setlocale("LC_ALL", ""), silent = TRUE)

ROOT_DIR <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT_DIR <- file.path(ROOT_DIR, "01_分析项目", "36_TITAN2致病菌ASV水平_水位")

TABLE_DIR <- file.path(PROJECT_DIR, "结果表格", "WL", "ASV")

PLOTDATA_DIR <- file.path(TABLE_DIR, "plotdata_figure2_top20_unique_species_ridgeline_latest_annotation_sunrise")

FIGURE_DIR <- file.path(PROJECT_DIR, "结果图片", "WL", "ASV", "02_Top20_species_ASV_ridge_sunrise_narrow_latest")

dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

DISPLAY_N <- 20L

ANNOTATION_FILE <- file.path(ROOT_DIR, "00_输入文件", "04_分类注释", "MBPD", "MBPD_ASV_taxonomy_LCA_2755.csv")

TOP1_ANNOTATION_FILE <- file.path(ROOT_DIR, "00_输入文件", "04_分类注释", "MBPD", "MBPD_ASV_taxonomy_top1_2755.csv")

.libPaths(unique(c(file.path(ROOT_DIR, "99_运行环境", "R_library"), .libPaths())))

required_packages <- c("ggplot2", "dplyr", "tidyr", "ragg", "svglite")

available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!available)) {
  stop("缺少R包：", paste(required_packages[!available], collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

read_result <- function(name) {
  path <- file.path(TABLE_DIR, name)
  if (!file.exists(path))
    stop("缺少正式TITAN2结果：", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

save_four_formats <- function(plot, prefix, width, height) {
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(paste0(prefix, ".tiff"), plot, width = width, height = height, dpi = 600, bg = "white", compression = "lzw",
    device = ragg::agg_tiff)
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height, bg = "white", device = svglite::svglite)
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, bg = "white", device = cairo_pdf)
}

normalize_rank <- function(x) {
  y <- trimws(as.character(x))
  bad <- is.na(x) | y == "" | tolower(y) %in% c("na", "n/a", "unknown", "unclassified", "unassigned",
    "uncultured", "unidentified", "not assigned", "not_assigned", "metagenome")
  y[bad] <- NA_character_
  y
}

species_epithet_label <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^['\"]+|['\"]+$", "", x)
  x <- sub("^Candidatus\\s+", "", x, ignore.case = TRUE)
  words <- strsplit(x, "\\s+")
  vapply(words, function(parts) {
    if (length(parts) >= 2L)
      parts[[2L]]
    else parts[[1L]]
  }, character(1L))
}

config <- read_result("04_analysis_config.csv")

cfg <- stats::setNames(as.character(config$value), config$parameter)

if (cfg[["run_mode"]] != "formal" || as.integer(cfg[["num_permutations"]]) != 500L || as.integer(cfg[["n_bootstraps"]]) !=
  500L || cfg[["taxon_level"]] != "ASV" || toupper(cfg[["taxonomy_annotation_used"]]) != "FALSE") {
  stop("山脊图只能使用未用分类注释的ASV级500/500正式结果。")
}

if (!file.exists(ANNOTATION_FILE)) stop("缺少最新MBPD注释文件。")

annotation <- read.csv(ANNOTATION_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_annotation <- c("ASV_ID", "Selection_Rule", "Best_Pident", "Best_Qcovs", "domain", "phylum",
  "class", "order", "family", "genus", "species")

if (!all(required_annotation %in% names(annotation))) {
  stop("最新MBPD注释缺少必要列。")
}

if (nrow(annotation) != 2755L || anyDuplicated(annotation$ASV_ID) || anyNA(annotation$ASV_ID) || any(trimws(annotation$ASV_ID) ==
  "")) {
  stop("最新MBPD注释必须是2755个唯一且非空的ASV ID。")
}

rank_columns <- c("domain", "phylum", "class", "order", "family", "genus", "species")

for (rank_i in rank_columns) {
  annotation[[rank_i]] <- normalize_rank(annotation[[rank_i]])
}

annotation$species_annotation_eligible <- !is.na(annotation$species) & !grepl("(^|\\s)bacterium(\\s|$)|(^|\\s)sp\\.(\\s|$)|oral\\s+taxon",
  annotation$species, ignore.case = TRUE, perl = TRUE)

annotation$species_dedup_key <- ifelse(annotation$species_annotation_eligible, tolower(gsub("\\s+", " ",
  trimws(annotation$species))), NA_character_)

annotation$short_ASV_ID <- substr(annotation$ASV_ID, 1L, 8L)

if (anyDuplicated(annotation$short_ASV_ID)) {
  stop("ASV前8位在最新注释中存在碰撞，不能作为短标签。")
}

genus_valid <- !is.na(annotation$genus)

species_valid <- !is.na(annotation$species)

species_without_candidatus <- sub("^Candidatus\\s+", "", annotation$species, ignore.case = TRUE)

annotation$genus_species_consistent <- ifelse(genus_valid & species_valid, startsWith(tolower(species_without_candidatus),
  paste0(tolower(annotation$genus), " ")) | tolower(species_without_candidatus) == tolower(annotation$genus),
  NA)

annotation$genus_species_conflict <- genus_valid & species_valid & !annotation$genus_species_consistent

annotation$deepest_annotation_level <- dplyr::case_when(!is.na(annotation$species) ~ "species", !is.na(annotation$genus) ~
  "genus", !is.na(annotation$family) ~ "family", !is.na(annotation$order) ~ "order", !is.na(annotation$class) ~
  "class", !is.na(annotation$phylum) ~ "phylum", !is.na(annotation$domain) ~ "domain", TRUE ~ "unresolved")

annotation$deepest_annotation <- dplyr::case_when(!is.na(annotation$species) ~ annotation$species, !is.na(annotation$genus) ~
  paste0(annotation$genus, " sp."), !is.na(annotation$family) ~ paste0(annotation$family, " (family)"),
  !is.na(annotation$order) ~ paste0(annotation$order, " (order)"), !is.na(annotation$class) ~ paste0(annotation$class,
    " (class)"), !is.na(annotation$phylum) ~ paste0(annotation$phylum, " (phylum)"), !is.na(annotation$domain) ~
    paste0(annotation$domain, " (domain)"), TRUE ~ "Unresolved ASV")

annotation$display_label <- paste0(annotation$deepest_annotation, ifelse(annotation$genus_species_conflict,
  "†", ""), " [", annotation$short_ASV_ID, "]")

annotation$annotation_method_interpreted <- "MBPD_single_best_hit_with_numeric_and_lexical_tie_breakers"

annotation_md5 <- unname(tools::md5sum(ANNOTATION_FILE))

top1_exists <- file.exists(TOP1_ANNOTATION_FILE)

top1_md5 <- if (top1_exists) {
  unname(tools::md5sum(TOP1_ANNOTATION_FILE))
} else {
  NA_character_
}

annotation_info <- file.info(ANNOTATION_FILE)

top1_info <- if (top1_exists) file.info(TOP1_ANNOTATION_FILE) else NULL

annotation_provenance <- data.frame(role = c("latest_annotation_used", "explicit_top1_comparison"), absolute_path = c(ANNOTATION_FILE,
  TOP1_ANNOTATION_FILE), exists = c(TRUE, top1_exists), bytes = c(annotation_info$size, if (top1_exists) top1_info$size else NA),
  last_modified = c(format(annotation_info$mtime, "%Y-%m-%d %H:%M:%S"), if (top1_exists) format(top1_info$mtime,
    "%Y-%m-%d %H:%M:%S") else NA), md5 = c(annotation_md5, top1_md5), identical_to_latest = c(TRUE,
    isTRUE(top1_md5 == annotation_md5)), interpretation = c(paste0("Latest file by modification time; Selection_Rule is a single-best-hit ",
    "rule, so labels are described as latest MBPD best-hit annotations."), "Explicit top1 file retained only for content-identity audit."),
  stringsAsFactors = FALSE)

robust_raw <- dplyr::mutate(read_result("19_robust_ASV_response_strength_ranking.csv"), taxon_id = as.character(taxon_id),
  ASV_ID = as.character(ASV_ID), robust_direction = as.character(robust_direction), change_point_observed = as.numeric(change_point_observed),
  change_point_q05 = as.numeric(change_point_q05), change_point_median = as.numeric(change_point_median),
  change_point_q95 = as.numeric(change_point_q95), zscore_bootstrap_median = as.numeric(zscore_bootstrap_median),
  permutation_probability_observed = as.numeric(permutation_probability_observed), purity = as.numeric(purity),
  reliability = as.numeric(reliability))

if (nrow(robust_raw) != 181L || anyDuplicated(robust_raw$ASV_ID)) {
  stop("正式稳健ASV排名表必须含181个唯一ASV。")
}

robust_recomputed <- dplyr::mutate(dplyr::arrange(robust_raw, dplyr::desc(zscore_bootstrap_median), permutation_probability_observed,
  dplyr::desc(purity), dplyr::desc(reliability), ASV_ID), response_strength_rank_recomputed = dplyr::row_number())

if (!all(robust_recomputed$response_strength_rank_recomputed == as.integer(robust_recomputed$response_strength_rank))) {
  stop("Top20响应强度排序不能从正式稳健ASV表精确复现。")
}

candidate_screen <- read_result("05_all_candidate_ASV_input_screening.csv")

if (!setequal(candidate_screen$ASV_ID, annotation$ASV_ID)) {
  stop("最新注释与2755个候选致病菌ASV不能一一对应。")
}

if (!all(robust_recomputed$ASV_ID %in% annotation$ASV_ID)) {
  stop("最新注释没有覆盖全部181个稳健ASV。")
}

annotation_keep <- dplyr::select(annotation, ASV_ID, short_ASV_ID, dplyr::all_of(rank_columns), deepest_annotation_level,
  deepest_annotation, display_label, species_annotation_eligible, species_dedup_key, genus_species_consistent,
  genus_species_conflict, All_Filtered_Hit_Rows, Top_Bitscore_Hit_Rows, Top_Bitscore_Distinct_Taxonomies,
  Final_Numeric_Tie_Rows, Final_Numeric_Tie_Distinct_Taxonomies, Best_Hit_Rows, Best_Distinct_Taxonomies,
  Selected_sseqid, Selection_Rule, Best_Bitscore, Best_Pident, Best_Qcovs, Best_Alignment_Length, annotation_method_interpreted)

robust_annotated <- dplyr::left_join(robust_recomputed, annotation_keep, by = "ASV_ID")

if (anyNA(robust_annotated$display_label)) {
  stop("稳健ASV连接最新注释后出现缺失。")
}

species_annotated_robust <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(dplyr::arrange(dplyr::filter(robust_annotated,
  species_annotation_eligible), dplyr::desc(zscore_bootstrap_median), permutation_probability_observed,
  dplyr::desc(purity), dplyr::desc(reliability), ASV_ID), species_dedup_key), rank_within_species = dplyr::row_number(),
  selected_as_species_representative = rank_within_species == 1L, ASVs_with_same_species_annotation = dplyr::n()))

species_representatives <- dplyr::mutate(dplyr::arrange(dplyr::filter(species_annotated_robust, selected_as_species_representative),
  dplyr::desc(zscore_bootstrap_median), permutation_probability_observed, dplyr::desc(purity), dplyr::desc(reliability),
  ASV_ID), species_response_strength_rank = dplyr::row_number())

if (nrow(species_representatives) < DISPLAY_N) {
  stop("有有效种注释且去重后的稳健ASV不足", DISPLAY_N, "个。")
}

top20 <- dplyr::mutate(dplyr::filter(species_representatives, species_response_strength_rank <= DISPLAY_N),
  selected_for_top20_plot = TRUE, plot_label = species_epithet_label(species), response_direction_label = dplyr::recode(robust_direction,
    `negative_z-` = "Negative", `positive_z+` = "Positive"))

if (nrow(top20) != DISPLAY_N || anyDuplicated(top20$species_dedup_key) || anyDuplicated(top20$plot_label) ||
  anyDuplicated(top20$short_ASV_ID)) {
  stop("Top20必须含20个不同物种、唯一代表ASV和唯一显示标签。")
}

if (sum(top20$response_direction_label == "Negative") != 4L || sum(top20$response_direction_label ==
  "Positive") != 16L) {
  stop("Top20响应方向数量与正式结果不一致。")
}

bootstrap_all <- dplyr::transmute(read_result("17_ASV_bootstrap_metrics_long.csv"), taxon_id = as.character(taxon_id),
  ASV_ID = as.character(ASV_ID), bootstrap_replicate = as.integer(bootstrap_replicate), direction = as.character(direction),
  change_point = as.numeric(change_point), zscore = as.numeric(zscore), permutation_probability = as.numeric(permutation_probability))

bootstrap_top20 <- dplyr::filter(dplyr::inner_join(bootstrap_all, dplyr::select(top20, taxon_id, ASV_ID,
  robust_direction, response_direction_label, species_response_strength_rank, display_label, plot_label,
  deepest_annotation_level, genus_species_conflict), by = c("taxon_id", "ASV_ID")), direction == robust_direction,
  is.finite(change_point))

bootstrap_audit <- dplyr::summarise(dplyr::group_by(bootstrap_top20, taxon_id, ASV_ID, display_label,
  plot_label, response_direction_label, species_response_strength_rank), N_bootstrap_matching_direction = dplyr::n(),
  N_unique_change_points = dplyr::n_distinct(change_point), bootstrap_min = min(change_point), bootstrap_q05_recomputed = stats::quantile(change_point,
    0.05, names = FALSE), bootstrap_median_recomputed = stats::median(change_point), bootstrap_q95_recomputed = stats::quantile(change_point,
    0.95, names = FALSE), bootstrap_max = max(change_point), .groups = "drop")

missing_bootstrap <- dplyr::anti_join(top20, bootstrap_audit, by = "ASV_ID")

if (nrow(missing_bootstrap) > 0L || any(bootstrap_audit$N_bootstrap_matching_direction < 475L) || any(bootstrap_audit$N_bootstrap_matching_direction >
  500L)) {
  stop("每个Top20 ASV必须有至少95%的同方向bootstrap阈值。")
}

quantile_check <- dplyr::mutate(dplyr::inner_join(dplyr::select(top20, ASV_ID, change_point_q05, change_point_median,
  change_point_q95), bootstrap_audit, by = "ASV_ID"), q05_absolute_difference = abs(change_point_q05 -
  bootstrap_q05_recomputed), median_absolute_difference = abs(change_point_median - bootstrap_median_recomputed),
  q95_absolute_difference = abs(change_point_q95 - bootstrap_q95_recomputed))

if (any(!is.finite(as.matrix(quantile_check[, c("q05_absolute_difference", "median_absolute_difference",
  "q95_absolute_difference")])))) {
  stop("Top20 bootstrap分位数差异审计出现非有限值。")
}

if (any(bootstrap_top20$change_point < 149.67 - 1e-08) || any(bootstrap_top20$change_point > 170 + 1e-08)) {
  stop("Top20 bootstrap阈值超出实测水位范围。")
}

ridge_spacing <- 1.26

ridge_height <- 1.1

x_limits <- c(149.5, 170.5)

bubble_x <- x_limits[1] - 0.4

x_grid <- seq(x_limits[1], x_limits[2], length.out = 512L)

common_bandwidth_m <- 0.34

ridge_summary <- dplyr::mutate(dplyr::arrange(dplyr::inner_join(top20, bootstrap_audit, by = c("taxon_id",
  "ASV_ID", "display_label", "plot_label", "response_direction_label", "species_response_strength_rank")),
  change_point_median, species), threshold_interval_width_m = change_point_q95 - change_point_q05,
  threshold_precision_per_m = 1/threshold_interval_width_m, Ridge_index = dplyr::row_number() * ridge_spacing)

if (any(!is.finite(ridge_summary$threshold_precision_per_m)) || any(ridge_summary$threshold_precision_per_m <=
  0)) {
  stop("Top20阈值精度必须是有限正数。")
}

ridge_list <- lapply(seq_len(nrow(ridge_summary)), function(i) {
  row_i <- ridge_summary[i, , drop = FALSE]
  values_i <- bootstrap_top20$change_point[bootstrap_top20$ASV_ID == row_i$ASV_ID]
  density_i <- stats::density(values_i, bw = common_bandwidth_m, from = x_limits[1], to = x_limits[2],
    n = length(x_grid), cut = 0)
  normalized_density <- if (max(density_i$y) > 0) {
    density_i$y/max(density_i$y)
  }
  else {
    rep(0, length(density_i$y))
  }
  curve <- data.frame(taxon_id = row_i$taxon_id, ASV_ID = row_i$ASV_ID, short_ASV_ID = row_i$short_ASV_ID,
    display_label = row_i$display_label, plot_label = row_i$plot_label, deepest_annotation_level = row_i$deepest_annotation_level,
    response_strength_rank = row_i$species_response_strength_rank, robust_direction = row_i$robust_direction,
    response_direction_label = row_i$response_direction_label, zscore_bootstrap_median = row_i$zscore_bootstrap_median,
    permutation_probability_observed = row_i$permutation_probability_observed, purity = row_i$purity,
    reliability = row_i$reliability, threshold_interval_width_m = row_i$threshold_interval_width_m,
    threshold_precision_per_m = row_i$threshold_precision_per_m, Ridge_index = row_i$Ridge_index,
    Water_level = density_i$x, Density = density_i$y, Density_normalized = normalized_density, Ridge_y = row_i$Ridge_index +
      ridge_height * normalized_density, stringsAsFactors = FALSE)
  polygon <- dplyr::bind_rows(dplyr::mutate(curve[1L, , drop = FALSE], Water_level = x_limits[1], Density = 0,
    Density_normalized = 0, Ridge_y = row_i$Ridge_index), curve, dplyr::mutate(curve[nrow(curve),
    , drop = FALSE], Water_level = x_limits[2], Density = 0, Density_normalized = 0, Ridge_y = row_i$Ridge_index))
  list(curve = curve, polygon = polygon)
})

curve_data <- dplyr::bind_rows(lapply(ridge_list, `[[`, "curve"))

polygon_data <- dplyr::bind_rows(lapply(ridge_list, `[[`, "polygon"))

ridge_summary$Median_density_height <- vapply(seq_len(nrow(ridge_summary)), function(i) {
  curve_i <- curve_data[curve_data$ASV_ID == ridge_summary$ASV_ID[i], , drop = FALSE]
  stats::approx(x = curve_i$Water_level, y = curve_i$Ridge_y, xout = ridge_summary$change_point_median[i],
    rule = 2)$y
}, numeric(1))

write_csv(annotation_provenance, file.path(PLOTDATA_DIR, "01_annotation_provenance.csv"))

write_csv(robust_annotated, file.path(PLOTDATA_DIR, "02_all_181_robust_ASV_latest_annotation.csv"))

write_csv(species_annotated_robust, file.path(PLOTDATA_DIR, "02a_species_annotated_robust_ASV_dedup_audit.csv"))

write_csv(species_representatives, file.path(PLOTDATA_DIR, "02b_unique_species_representatives_response_ranking.csv"))

write_csv(top20, file.path(PLOTDATA_DIR, "03_selected_top20_unique_species_complete_results_annotation.csv"))

write_csv(bootstrap_top20, file.path(PLOTDATA_DIR, "04_top20_bootstrap_thresholds_long.csv"))

write_csv(bootstrap_audit, file.path(PLOTDATA_DIR, "05_top20_bootstrap_audit.csv"))

write_csv(quantile_check, file.path(PLOTDATA_DIR, "06_top20_bootstrap_quantile_check.csv"))

write_csv(ridge_summary, file.path(PLOTDATA_DIR, "07_top20_ridgeline_summary.csv"))

write_csv(curve_data, file.path(PLOTDATA_DIR, "08_top20_density_curve_plotdata.csv"))

write_csv(polygon_data, file.path(PLOTDATA_DIR, "09_top20_density_polygon_plotdata.csv"))

annotation_conflict_audit <- dplyr::arrange(dplyr::select(top20, species_response_strength_rank, response_strength_rank_recomputed,
  ASV_ID, short_ASV_ID, display_label, plot_label, deepest_annotation_level, genus, species, genus_species_consistent,
  genus_species_conflict, Top_Bitscore_Distinct_Taxonomies, Final_Numeric_Tie_Distinct_Taxonomies,
  Best_Distinct_Taxonomies, Selection_Rule, Best_Pident, Best_Qcovs, Selected_sseqid), species_response_strength_rank)

write_csv(annotation_conflict_audit, file.path(PLOTDATA_DIR, "10_top20_annotation_conflict_audit.csv"))

direction_colors <- c(Negative = "#2C7FB8", Positive = "#C62F43")

sunrise_colors <- colorRampPalette(c("#8075AD", "#F5EDF0", "#F5DA73", "#FFB800"))(100)

label_breaks <- ridge_summary$Ridge_index

label_values <- ridge_summary$plot_label

size_range <- range(ridge_summary$zscore_bootstrap_median, na.rm = TRUE)

size_breaks <- unique(seq(size_range[1L], size_range[2L], length.out = 4L))

precision_limits <- range(ridge_summary$threshold_precision_per_m, na.rm = TRUE)

ridge_plot <- ggplot() + geom_segment(data = ridge_summary, aes(x = x_limits[1L], xend = x_limits[2L],
  y = Ridge_index, yend = Ridge_index), color = "#D9D9D9", linewidth = 0.28) + geom_polygon(data = polygon_data,
  aes(x = Water_level, y = Ridge_y, group = ASV_ID, fill = threshold_precision_per_m), alpha = 0.82,
  color = NA) + geom_line(data = curve_data, aes(x = Water_level, y = Ridge_y, group = ASV_ID), color = "#5E4B44",
  linewidth = 0.58, lineend = "round") + geom_segment(data = ridge_summary, aes(x = change_point_median,
  xend = change_point_median, y = Ridge_index + 0.015, yend = Median_density_height), color = "white",
  linewidth = 0.72, lineend = "round") + geom_point(data = ridge_summary, aes(x = bubble_x, y = Ridge_index +
  0.04, size = zscore_bootstrap_median, color = response_direction_label), shape = 16, alpha = 0.98) +
  scale_color_manual(values = direction_colors, breaks = c("Negative", "Positive"), name = "Direction") +
  scale_fill_gradientn(colors = sunrise_colors, limits = precision_limits, name = expression("Threshold precision" ~
    (m^{
      -1
    })), guide = guide_colorbar(order = 1, title.position = "top", title.hjust = 0.5, barwidth = grid::unit(40,
    "mm"), barheight = grid::unit(3, "mm"))) + scale_size_continuous(range = c(1.4, 4.2), breaks = size_breaks,
  labels = function(x) format(round(x, 1), nsmall = 1, trim = TRUE), name = "Median z") + scale_x_continuous(breaks = seq(150,
  170, 5), expand = c(0, 0)) + scale_y_continuous(breaks = label_breaks, labels = label_values, expand = expansion(add = c(0.68,
  1.58))) + labs(x = "Water level (m)", y = NULL) + coord_cartesian(xlim = x_limits, clip = "off") +
  theme_classic(base_size = 11.8, base_family = "sans") + theme(axis.line = element_blank(), panel.border = element_rect(fill = NA,
  color = "black", linewidth = 0.55), axis.ticks.y = element_blank(), axis.ticks.x = element_line(color = "black",
  linewidth = 0.42), axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "#3F3F3F",
  size = 16, hjust = 1, margin = margin(r = 18)), axis.title.x = element_text(face = "bold", size = 14.5),
  panel.grid.major.x = element_line(color = "#E2E2E2", linewidth = 0.35), panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "bottom",
  legend.direction = "horizontal", legend.box = "vertical", legend.box.just = "center", legend.spacing.y = grid::unit(1.5,
    "pt"), legend.box.spacing = grid::unit(2, "pt"), legend.title = element_text(face = "bold", size = 11.5),
  legend.text = element_text(size = 11), legend.key = element_rect(fill = "transparent", color = NA),
  plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white",
    color = NA), plot.margin = margin(8, 10, 8, 8)) + guides(color = guide_legend(order = 2, nrow = 1,
  override.aes = list(size = 3.2)), size = guide_legend(order = 3, nrow = 1))

figure_width <- 8

figure_height <- 12.8

output_prefix <- file.path(FIGURE_DIR, "TITAN2_WL_Top20_species_ASV_ridge_sunrise_narrow")

save_four_formats(ridge_plot, output_prefix, width = figure_width, height = figure_height)

expected_figure_files <- paste0(output_prefix, c(".png", ".tiff", ".svg", ".pdf"))

qa <- data.frame(check = c("formal_500_permutations", "formal_500_bootstraps", "taxonomy_not_used_in_TITAN2",
  "annotation_rows_2755", "annotation_ASV_ID_unique", "candidate_annotation_coverage", "robust_annotation_coverage",
  "top20_annotation_coverage", "species_annotated_robust_ASVs_93", "unique_species_representatives_48",
  "top20_only_explicit_species_annotations", "top20_one_representative_per_species", "top20_representative_is_strongest_within_species",
  "top20_species_response_ranks_1_to_20", "top20_labels_unique", "top20_short_IDs_unique", "top20_direction_counts_4_negative_16_positive",
  "each_top20_ASV_has_at_least_95pct_same_direction_bootstraps", "same_direction_bootstrap_quantile_differences_audited",
  "bootstrap_thresholds_within_measured_gradient", "observed_permutation_P_all_equal_0.002", "sunrise_fill_equals_threshold_precision",
  "source_genus_species_conflicts_recorded", "figure_is_narrow_and_tall", "four_figure_formats_nonempty"),
  pass = c(as.integer(cfg[["num_permutations"]]) == 500L, as.integer(cfg[["n_bootstraps"]]) == 500L,
    toupper(cfg[["taxonomy_annotation_used"]]) == "FALSE", nrow(annotation) == 2755L, !anyDuplicated(annotation$ASV_ID),
    setequal(candidate_screen$ASV_ID, annotation$ASV_ID), all(robust_recomputed$ASV_ID %in% annotation$ASV_ID),
    all(top20$ASV_ID %in% annotation$ASV_ID), nrow(species_annotated_robust) == 93L, nrow(species_representatives) ==
      48L, all(top20$species_annotation_eligible), !anyDuplicated(top20$species_dedup_key), all(top20$selected_as_species_representative) &&
      all(top20$rank_within_species == 1L), identical(top20$species_response_strength_rank, 1:20),
    !anyDuplicated(top20$plot_label), !anyDuplicated(top20$short_ASV_ID), sum(top20$response_direction_label ==
      "Negative") == 4L && sum(top20$response_direction_label == "Positive") == 16L, all(bootstrap_audit$N_bootstrap_matching_direction >=
      475L) && all(bootstrap_audit$N_bootstrap_matching_direction <= 500L), all(is.finite(as.matrix(quantile_check[,
      c("q05_absolute_difference", "median_absolute_difference", "q95_absolute_difference")]))),
    min(bootstrap_top20$change_point) >= 149.67 - 1e-08 && max(bootstrap_top20$change_point) <= 170 +
      1e-08, length(unique(top20$permutation_probability_observed)) == 1L && abs(unique(top20$permutation_probability_observed) -
      0.002) <= 1e-12, all(abs(curve_data$threshold_precision_per_m - 1/curve_data$threshold_interval_width_m) <=
      1e-12), sum(top20$genus_species_conflict) == 4L, figure_height > figure_width, all(file.exists(expected_figure_files)) &&
      all(file.info(expected_figure_files)$size > 0)), expected_or_note = c("500", "500", "FALSE",
    "2755", "unique", "2755/2755", "181/181", "20/20", "93 explicit species-annotated robust ASVs",
    "48 species after exact annotation de-duplication", "no NA/sp./bacterium/oral-taxon placeholders",
    "20 distinct species", "z/P/purity/reliability/ASV tie-break", "1:20", "20", "20", "4 negative; 16 positive",
    "478-500 each", "formal all-bootstrap vs same-direction-only differences retained", "149.67-170.00 m",
    "all observed P = 0.002, so P cannot distinguish ridge colors", "fill = 1 / (bootstrap q95 - q05), in inverse metres",
    "4 conflicts marked with dagger and retained in audit", "8.0 x 12.8 inches", "PNG/TIFF/SVG/PDF"),
  stringsAsFactors = FALSE)

write_csv(qa, file.path(PLOTDATA_DIR, "11_figure_QA.csv"))

mapping_notes <- data.frame(visual_element = c("ridge_fill", "left_bubble_color", "left_bubble_size",
  "white_vertical_line", "ridge_height"), mapped_metric = c("threshold_precision_per_m", "robust_direction",
  "zscore_bootstrap_median", "change_point_median", "within_ASV_peak_normalized_bootstrap_density"),
  interpretation_cn = c(paste0("阈值精度=1/(bootstrap q95-q05)，单位1/m；值越大，", "90% bootstrap阈值区间越窄、阈值越稳定。Top20观测P均为0.002，",
    "所以不再用P填色。"), "蓝色为负响应，红色为正响应。", "气泡越大，bootstrap median z越大，响应越强。",
    "白线为同方向bootstrap阈值中位数。", "每个ASV单独峰值归一化；山脊高度不是丰度，也不能跨ASV比较总量。"),
  stringsAsFactors = FALSE)

write_csv(mapping_notes, file.path(PLOTDATA_DIR, "12_plot_visual_mapping_notes.csv"))

guide <- data.frame(file = c("01_annotation_provenance.csv", "02_all_181_robust_ASV_latest_annotation.csv",
  "02a_species_annotated_robust_ASV_dedup_audit.csv", "02b_unique_species_representatives_response_ranking.csv",
  "03_selected_top20_unique_species_complete_results_annotation.csv", "04_top20_bootstrap_thresholds_long.csv",
  "05_top20_bootstrap_audit.csv", "06_top20_bootstrap_quantile_check.csv", "07_top20_ridgeline_summary.csv",
  "08_top20_density_curve_plotdata.csv", "09_top20_density_polygon_plotdata.csv", "10_top20_annotation_conflict_audit.csv",
  "11_figure_QA.csv", "12_plot_visual_mapping_notes.csv"), purpose_cn = c("最新MBPD注释的绝对路径、时间、MD5及与top1文件的同一性审计",
  "全部181个稳健ASV的正式TITAN结果、强度排名和最新注释", "93个明确种注释稳健ASV；含同种内排名、代表ASV选择和被去重记录",
  "48个不同种的代表稳健ASV及重新计算的响应强度排名", "图中Top20的全部阈值、方向、z、P、稳健性、分类和显示标签",
  "Top20每个ASV的500次bootstrap阈值、z和P", "每个Top20 ASV的bootstrap数量、唯一阈值点、范围和重算分位数",
  "bootstrap重算5%、中位数、95%与正式TITAN2结果的逐ASV核对", "每条山脊的纵轴位置、显示标签、阈值、排名、z、P和稳健指标",
  "山脊KDE曲线逐点数据；Density_normalized为ASV内峰值归一化", "闭合山脊填充多边形逐点数据，可直接重画",
  "Top20原始属种字段、层级一致性、命中冲突和best-hit质量审计", "正式结果、注释连接、Top20、bootstrap、阈值、标签和四格式图件QA",
  "日出填色、方向点、Median z气泡、白线和山脊高度的逐项含义"), stringsAsFactors = FALSE)

write_csv(guide, file.path(PLOTDATA_DIR, "00_plotdata_guide.csv"))

capture.output(sessionInfo(), file = file.path(PLOTDATA_DIR, "sessionInfo_plot.txt"))

if (any(!qa$pass)) {
  stop("Top20稳健ASV山脊图QA未全部通过；请查看11_figure_QA.csv。")
}

message("完成Top20不同种代表稳健ASV窄高山脊图：", "181个稳健ASV中93个有明确种注释，去重后48种；",
  "Top20为4个负响应、16个正响应；物种标签不显示冲突符号；", "sunrise填色已改为阈值精度，四种格式和全部QA均通过。")
