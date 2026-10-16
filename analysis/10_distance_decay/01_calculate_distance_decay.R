rm(list = ls())

set.seed(123)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- Sys.getenv("DJK16S_ACTIVE_ROOT", unset = file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA/decay"))

project_dir <- dirname(script_dir)

tpm_file <- file.path(script_dir, "tpm.csv")

env_file <- file.path(project_dir, "env.csv")

group_size_file <- file.path(script_dir, "groupsize.csv")

group_month_file <- file.path(script_dir, "groupmonth.csv")

coord_file <- file.path(script_dir, "1.csv")

output_dir <- file.path(script_dir, "outputs", "decay_analysis")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

MIN_TOTAL_ABUNDANCE <- 0.1

N_PERM <- 9999

MANTEL_METHOD <- "spearman"

RUN_PARTIAL_MANTEL <- TRUE

pkgs <- c("ggplot2", "dplyr", "tidyr", "stringr", "vegan", "geosphere", "scales", "patchwork", "svglite",
  "ragg")

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("缺少 R 包：", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(vegan)
  library(geosphere)
  library(scales)
  library(patchwork)
})

read_csv_smart <- function(file, ...) {
  out <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM",
    ...), error = function(e) NULL)
  if (is.null(out)) {
    out <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "GBK", ...)
  }
  out
}

parse_month_index <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0$", "", x)
  x <- stringr::str_pad(x, width = 4, side = "left", pad = "0")
  yy <- as.integer(substr(x, 1, 2))
  mm <- as.integer(substr(x, 3, 4))
  year <- ifelse(yy < 50, 2000 + yy, 1900 + yy)
  year * 12 + mm
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

safe_rescale_01 <- function(x) {
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1) {
    return(rep(0.5, length(x)))
  }
  scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
}

save_plot_all_formats <- function(p, file_prefix, width, height) {
  png_file <- file.path(output_dir, paste0(file_prefix, ".png"))
  svg_file <- file.path(output_dir, paste0(file_prefix, ".svg"))
  pdf_file <- file.path(output_dir, paste0(file_prefix, ".pdf"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_file, p, width = width, height = height, dpi = 300, device = ragg::agg_png)
  }
  else {
    ggsave(png_file, p, width = width, height = height, dpi = 300)
  }
  ggsave(svg_file, p, width = width, height = height, device = "svg")
  ggsave(pdf_file, p, width = width, height = height, device = cairo_pdf)
}

env <- read_csv_smart(env_file)

group_size <- read_csv_smart(group_size_file)

group_month <- read_csv_smart(group_month_file)

coord_raw <- read_csv_smart(coord_file)

names(group_size)[1:2] <- c("sample", "size_group")

names(group_month)[1:2] <- c("sample", "month_group")

env$sample <- trimws(as.character(env$sample))

group_size$sample <- trimws(as.character(group_size$sample))

group_month$sample <- trimws(as.character(group_month$sample))

group_size$size_group <- tolower(trimws(as.character(group_size$size_group)))

site_col <- names(coord_raw)[1]

lon_col <- names(coord_raw)[grepl("经度|longitude|lon", names(coord_raw), ignore.case = TRUE)][1]

lat_col <- names(coord_raw)[grepl("纬度|latitude|lat", names(coord_raw), ignore.case = TRUE)][1]

if (is.na(lon_col) || is.na(lat_col)) {
  stop("经纬度文件中没有识别到经度/纬度列，请检查列名是否包含 经度/纬度 或 longitude/latitude。")
}

coord <- coord_raw %>%
  transmute(site = trimws(as.character(.data[[site_col]])), longitude = as.numeric(.data[[lon_col]]),
    latitude = as.numeric(.data[[lat_col]])) %>%
  filter(!is.na(site), site != "", !is.na(longitude), !is.na(latitude))

tpm_raw <- read_csv_smart(tpm_file)

tax_col <- names(tpm_raw)[1]

tax_id <- make.unique(as.character(tpm_raw[[tax_col]]))

tpm_raw[[tax_col]] <- NULL

tpm_mat <- as.matrix(tpm_raw)

mode(tpm_mat) <- "numeric"

rownames(tpm_mat) <- tax_id

tpm_mat <- tpm_mat[rowSums(tpm_mat, na.rm = TRUE) > MIN_TOTAL_ABUNDANCE, , drop = FALSE]

abund <- as.data.frame(t(tpm_mat), check.names = FALSE)

abund$sample <- rownames(abund)

meta <- env %>%
  inner_join(group_month, by = "sample") %>%
  inner_join(group_size, by = "sample") %>%
  mutate(site = stringr::str_remove(sample, "\\d+$"), month_index = parse_month_index(month_group)) %>%
  left_join(coord, by = "site")

if (any(is.na(meta$longitude)) || any(is.na(meta$latitude))) {
  missing_site <- meta %>%
    filter(is.na(longitude) | is.na(latitude)) %>%
    distinct(site) %>%
    pull(site)
  stop(paste0("以下站点没有在经纬度文件中匹配到坐标：", paste(missing_site, collapse = ", ")))
}

common_samples <- Reduce(intersect, list(meta$sample, abund$sample))

if (length(common_samples) < 3) {
  stop("可用于分析的共同样本数少于 3，请检查 tpm/env/group 文件中的 sample 是否一致。")
}

meta <- meta %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

abund <- abund %>%
  filter(sample %in% meta$sample) %>%
  arrange(match(sample, meta$sample))

rownames(abund) <- abund$sample

abund$sample <- NULL

env_vars <- setdiff(names(env), c("sample", "Longitude", "Latitude", "Q"))

env_vars <- env_vars[sapply(env[env_vars], function(x) suppressWarnings(!all(is.na(as.numeric(x)))))]

message("使用的环境因子：", paste(env_vars, collapse = ", "))

calc_decay_one_group <- function(group_name) {
  if (group_name == "all") {
    sub_meta <- meta
  }
  else {
    sub_meta <- meta %>%
      filter(size_group == group_name)
  }
  sub_abund <- abund[sub_meta$sample, , drop = FALSE]
  sub_abund <- sub_abund[, colSums(sub_abund, na.rm = TRUE) > 0, drop = FALSE]
  if (nrow(sub_abund) < 3) {
    stop(paste0(group_name, " 组样本数少于 3，无法进行 Mantel test。"))
  }
  if (ncol(sub_abund) < 2) {
    stop(paste0(group_name, " 组有效类群数少于 2，无法计算 Bray-Curtis。"))
  }
  dist_abund <- vegdist(sub_abund, method = "bray")
  env_sub <- sub_meta[, env_vars, drop = FALSE] %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))
  keep_env <- names(env_sub)[sapply(env_sub, function(x) {
    sum(!is.na(x)) >= 3 && sd(x, na.rm = TRUE) > 0
  })]
  env_sub <- env_sub[, keep_env, drop = FALSE]
  env_sub <- env_sub %>%
    mutate(across(everything(), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))
  scale_env <- scale(env_sub, center = TRUE, scale = TRUE)
  dist_env <- dist(scale_env, method = "euclidean")
  geo_xy <- as.matrix(sub_meta[, c("longitude", "latitude")])
  dist_geo_m <- as.dist(geosphere::distm(geo_xy, fun = geosphere::distHaversine))
  dist_time <- as.dist(abs(outer(sub_meta$month_index, sub_meta$month_index, "-")))
  mantel_env <- vegan::mantel(dist_abund, dist_env, method = MANTEL_METHOD, permutations = N_PERM,
    na.rm = TRUE)
  mantel_geo <- vegan::mantel(dist_abund, dist_geo_m, method = MANTEL_METHOD, permutations = N_PERM,
    na.rm = TRUE)
  mantel_time <- vegan::mantel(dist_abund, dist_time, method = MANTEL_METHOD, permutations = N_PERM,
    na.rm = TRUE)
  mantel_tbl <- tibble(Group = group_name, N_sample = nrow(sub_meta), N_taxa = ncol(sub_abund), Distance_code = c("Environmental_distance",
    "Geographic_distance_km", "Temporal_distance_month"), Distance_type = c("Environmental distance",
    "Geographic distance", "Temporal distance"), Mantel_r = c(unname(mantel_env$statistic), unname(mantel_geo$statistic),
    unname(mantel_time$statistic)), P_value = c(mantel_env$signif, mantel_geo$signif, mantel_time$signif),
    P_label = fmt_p(c(mantel_env$signif, mantel_geo$signif, mantel_time$signif)), Env_vars_used = paste(keep_env,
      collapse = "; "))
  pair_tbl <- tibble(Group = group_name, Pair_id = seq_along(as.vector(dist_abund)), Bray_Curtis = as.vector(dist_abund),
    Environmental_distance = as.vector(dist_env), Geographic_distance_km = as.vector(dist_geo_m)/1000,
    Temporal_distance_month = as.vector(dist_time))
  partial_tbl <- NULL
  if (RUN_PARTIAL_MANTEL) {
    pe_t <- vegan::mantel.partial(dist_abund, dist_env, dist_time, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    pt_e <- vegan::mantel.partial(dist_abund, dist_time, dist_env, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    pe_g <- vegan::mantel.partial(dist_abund, dist_env, dist_geo_m, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    pg_e <- vegan::mantel.partial(dist_abund, dist_geo_m, dist_env, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    pt_g <- vegan::mantel.partial(dist_abund, dist_time, dist_geo_m, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    pg_t <- vegan::mantel.partial(dist_abund, dist_geo_m, dist_time, method = MANTEL_METHOD, permutations = N_PERM,
      na.rm = TRUE)
    partial_tbl <- tibble(Group = group_name, Test = c("Environmental | Temporal", "Temporal | Environmental",
      "Environmental | Geographic", "Geographic | Environmental", "Temporal | Geographic", "Geographic | Temporal"),
      Partial_Mantel_r = c(unname(pe_t$statistic), unname(pt_e$statistic), unname(pe_g$statistic),
        unname(pg_e$statistic), unname(pt_g$statistic), unname(pg_t$statistic)), P_value = c(pe_t$signif,
        pt_e$signif, pe_g$signif, pg_e$signif, pt_g$signif, pg_t$signif), P_label = fmt_p(c(pe_t$signif,
        pt_e$signif, pe_g$signif, pg_e$signif, pt_g$signif, pg_t$signif)))
  }
  list(mantel = mantel_tbl, pair = pair_tbl, partial = partial_tbl)
}

group_order <- c("large", "small", "all")

res_list <- lapply(group_order, calc_decay_one_group)

mantel_results <- bind_rows(lapply(res_list, `[[`, "mantel")) %>%
  mutate(Group = factor(Group, levels = group_order), Distance_type = factor(Distance_type, levels = c("Environmental distance",
    "Geographic distance", "Temporal distance"))) %>%
  group_by(Group) %>%
  mutate(P_adj_BH = p.adjust(P_value, method = "BH"), P_label_BH = fmt_p(P_adj_BH)) %>%
  ungroup() %>%
  arrange(Group, desc(Mantel_r))

pair_wide <- bind_rows(lapply(res_list, `[[`, "pair")) %>%
  mutate(Group = factor(Group, levels = group_order))

partial_results <- bind_rows(lapply(res_list, `[[`, "partial")) %>%
  mutate(Group = factor(Group, levels = group_order)) %>%
  group_by(Group) %>%
  mutate(P_adj_BH = p.adjust(P_value, method = "BH")) %>%
  ungroup()

write.csv(mantel_results, file.path(output_dir, "mantel_decay_results.csv"), row.names = FALSE)

write.csv(pair_wide, file.path(output_dir, "decay_pairwise_plot_data.csv"), row.names = FALSE)

if (RUN_PARTIAL_MANTEL) {
  write.csv(partial_results, file.path(output_dir, "partial_mantel_results.csv"), row.names = FALSE)
}

driver_rank <- mantel_results %>%
  group_by(Group) %>%
  arrange(desc(Mantel_r), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

write.csv(driver_rank, file.path(output_dir, "driver_rank_by_mantel_r.csv"), row.names = FALSE)

message("\nMantel 结果：")

print(mantel_results)

message("\n按 Mantel r 从大到小排序：")

print(driver_rank %>%
  select(Group, Distance_type, Mantel_r, P_value, P_adj_BH, Rank))

plot_long <- pair_wide %>%
  pivot_longer(cols = c(Environmental_distance, Geographic_distance_km, Temporal_distance_month), names_to = "Distance_code",
    values_to = "Distance") %>%
  mutate(Distance_type = recode(Distance_code, Environmental_distance = "Environmental distance", Geographic_distance_km = "Geographic distance",
    Temporal_distance_month = "Temporal distance"), Distance_type = factor(Distance_type, levels = c("Environmental distance",
    "Geographic distance", "Temporal distance"))) %>%
  group_by(Group, Distance_type) %>%
  mutate(Distance_scaled = safe_rescale_01(Distance)) %>%
  ungroup()

write.csv(plot_long, file.path(output_dir, "decay_plot_long_data.csv"), row.names = FALSE)

anno_df <- plot_long %>%
  group_by(Group, Distance_type) %>%
  summarise(x_pos = min(Distance, na.rm = TRUE) + 0.97 * (max(Distance, na.rm = TRUE) - min(Distance,
    na.rm = TRUE)), y_pos = min(Bray_Curtis, na.rm = TRUE) + 0.97 * (max(Bray_Curtis, na.rm = TRUE) -
    min(Bray_Curtis, na.rm = TRUE)), .groups = "drop") %>%
  left_join(mantel_results %>%
    select(Group, Distance_type, Mantel_r, P_value, P_adj_BH, P_label_BH), by = c("Group", "Distance_type")) %>%
  mutate(label = paste0("Mantel r = ", sprintf("%.3f", Mantel_r), "\nBH P = ", P_label_BH))

group_labels <- c(large = "Large group", small = "Small group", all = "All samples")

distance_labels <- c(`Environmental distance` = "Environmental distance\nstandardized Euclidean", `Geographic distance` = "Geographic distance\nkm",
  `Temporal distance` = "Temporal distance\nmonths")

base_theme <- theme_bw(base_size = 12) + theme(panel.grid = element_blank(), panel.border = element_rect(color = "black",
  linewidth = 0.5), strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  strip.text = element_text(face = "bold", color = "black"), axis.text = element_text(color = "black"),
  axis.title = element_text(face = "bold", color = "black"), legend.position = "top", legend.title = element_text(face = "bold"),
  plot.title = element_text(face = "bold", hjust = 0.5))

p_all_decay <- ggplot(plot_long, aes(x = Distance, y = Bray_Curtis)) + geom_point(aes(fill = Distance_scaled),
  shape = 21, color = "black", size = 2.4, alpha = 0.7, stroke = 0.2) + geom_smooth(method = "lm",
  formula = y ~ x, se = TRUE, color = "black", fill = "grey80", linewidth = 0.75, alpha = 0.35) + geom_text(data = anno_df,
  aes(x = x_pos, y = y_pos, label = label), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3.6,
  fontface = "bold", color = "black") + facet_grid(Group ~ Distance_type, scales = "free_x", labeller = labeller(Group = group_labels,
  Distance_type = distance_labels)) + scale_fill_gradient(low = "#9ECAE1", high = "#CB181D", name = "Distance\ngradient") +
  labs(x = NULL, y = "Bray-Curtis dissimilarity", title = "Distance-decay relationships of microbial community") +
  base_theme

save_plot_all_formats(p_all_decay, "decay_3groups_3distances", width = 12, height = 9)

for (g in group_order) {
  dat_g <- plot_long %>%
    filter(Group == g)
  anno_g <- anno_df %>%
    filter(Group == g)
  p_g <- ggplot(dat_g, aes(x = Distance, y = Bray_Curtis)) + geom_point(aes(fill = Distance_scaled),
    shape = 21, color = "black", size = 2.8, alpha = 0.72, stroke = 0.22) + geom_smooth(method = "lm",
    formula = y ~ x, se = TRUE, color = "black", fill = "grey80", linewidth = 0.8, alpha = 0.35) +
    geom_text(data = anno_g, aes(x = x_pos, y = y_pos, label = label), inherit.aes = FALSE, hjust = 1,
      vjust = 1, size = 4, fontface = "bold", color = "black") + facet_wrap(~Distance_type, scales = "free_x",
    nrow = 1, labeller = labeller(Distance_type = distance_labels)) + scale_fill_gradient(low = "#9ECAE1",
    high = "#CB181D", name = "Distance\ngradient") + labs(x = NULL, y = "Bray-Curtis dissimilarity",
    title = paste0(group_labels[[g]], ": distance-decay relationships")) + base_theme
  save_plot_all_formats(p_g, paste0("decay_", g), width = 12, height = 4.2)
}

message("\n分析完成！输出目录：", output_dir)

message("建议优先查看：")

message("1) decay_3groups_3distances.png")

message("2) mantel_decay_results.csv")

message("3) driver_rank_by_mantel_r.csv")

if (RUN_PARTIAL_MANTEL) {
  message("4) partial_mantel_results.csv")
}
