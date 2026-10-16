set.seed(123)

project_root <- Sys.getenv("DJK16S_PROJECT_ROOT", unset = Sys.getenv("DJK_PROJECT_ROOT", unset = "."))

.libPaths(unique(c(file.path(project_root, "work", "r_packages_extra"), file.path(project_root, "work",
  "figure2_hydro_core", "Rlib"), .libPaths())))

suppressMessages(suppressWarnings(library(dplyr)))

suppressMessages(suppressWarnings(library(vegan)))

suppressMessages(suppressWarnings(library(linkET)))

suppressMessages(suppressWarnings(library(ggplot2)))

rcorr_pearson <- function(x) {
  x <- as.matrix(x)
  p <- ncol(x)
  r <- cor(x, use = "pairwise.complete.obs", method = "pearson")
  n <- matrix(0L, p, p, dimnames = dimnames(r))
  P <- matrix(NA_real_, p, p, dimnames = dimnames(r))
  for (i in seq_len(p)) for (j in seq_len(p)) {
    ok <- is.finite(x[, i]) & is.finite(x[, j])
    n[i, j] <- sum(ok)
    if (i != j && n[i, j] >= 3L) {
      P[i, j] <- suppressWarnings(cor.test(x[ok, i], x[ok, j], method = "pearson")$p.value)
    }
  }
  structure(list(r = r, n = n, P = P), class = "rcorr")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- Sys.getenv("DJK16S_ACTIVE_ROOT", unset = "")

if (!nzchar(script_dir)) {
  script_dir <- file.path(project_root, "work", "mantel_nfr_hrt_airt_inputs", "whole")
}

project_dir <- dirname(script_dir)

tpm_file <- Sys.getenv("DJK16S_WHOLE_TPM_FILE", unset = file.path(script_dir, "tpm_rarefied_21000.csv"))

env_file <- Sys.getenv("DJK16S_ENV_FILE", unset = file.path(project_root, "outputs", "20260825_env_hydro_air_update",
  "tables", "environment_analysis_120.csv"))

group_file <- file.path(script_dir, "groupsize.csv")

output_dir <- Sys.getenv("DJK16S_OUTPUT_DIR", unset = file.path(project_root, "outputs", "20260825_env_hydro_air_update",
  "analyses", "AllCommunity", "mantel_link_triangle"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tpm <- read.csv(tpm_file, row.names = 1, check.names = FALSE)

env <- read.csv(env_file, row.names = 1, check.names = FALSE)

groupsize <- read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE)

rarefied_depth <- colSums(tpm, na.rm = TRUE)

if (ncol(tpm) != 120L || any(rarefied_depth != 21000)) {
  stop("Rarefaction audit failed: expected 120 samples and exactly ", "21,000 reads per sample.")
}

message("Rarefied community table verified: ", ncol(tpm), " samples; depth = ", unique(rarefied_depth),
  " reads/sample")

if (!all(c("sample", "group") %in% colnames(groupsize))) {
  stop("groupsize(2).csv 必须包含 sample 和 group 两列。")
}

spec <- as.data.frame(t(tpm), check.names = FALSE)

spec[] <- lapply(spec, function(x) as.numeric(as.character(x)))

numeric_candidate_cols <- setdiff(colnames(env), c("MonthGroup", "SizeGroup"))

env_numeric <- as.data.frame(lapply(env[, numeric_candidate_cols, drop = FALSE], function(x) as.numeric(as.character(x))),
  row.names = rownames(env), check.names = FALSE)

original_env_vars <- c("WT", "pH", "DO", "CODmn", "BOD5", "NH3N", "TN", "TP", "Cu", "Zn", "F", "Se",
  "As", "Pb", "Chla", "WL")

added_env_vars <- c("NFR", "HRT", "AIRT")

mantel_env_vars <- c(original_env_vars, added_env_vars)

plot_env_vars <- c("WT", "AIRT", setdiff(mantel_env_vars, c("WT", "AIRT")))

missing_env_vars <- setdiff(mantel_env_vars, colnames(env_numeric))

if (length(missing_env_vars)) {
  stop("Updated environment table is missing: ", paste(missing_env_vars, collapse = ", "))
}

env_numeric <- env_numeric[, mantel_env_vars, drop = FALSE]

message("Environmental variables retained: ", paste(colnames(env_numeric), collapse = ", "))

env_numeric <- env_numeric[, colSums(!is.na(env_numeric)) > 0, drop = FALSE]

env_numeric <- env_numeric[, sapply(env_numeric, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

spec <- spec[, colSums(spec, na.rm = TRUE) > 0, drop = FALSE]

common_samples <- Reduce(intersect, list(rownames(spec), rownames(env_numeric), groupsize$sample))

if (length(common_samples) == 0) {
  stop("错误：tpm.csv、env(1).csv 和 groupsize(2).csv 之间没有匹配的样本名。")
}

spec <- spec[common_samples, , drop = FALSE]

env_numeric <- env_numeric[common_samples, , drop = FALSE]

groupsize <- groupsize[match(common_samples, groupsize$sample), , drop = FALSE]

groupsize$group <- tolower(groupsize$group)

if (!all(groupsize$group %in% c("small", "large"))) {
  stop("groupsize(2).csv 的 group 列只能包含 small 和 large。")
}

message("样本数：all = ", length(common_samples), "；small = ", sum(groupsize$group == "small"),
  "；large = ", sum(groupsize$group == "large"))

complete_env_samples <- rownames(env_numeric)[complete.cases(env_numeric)]

spec <- spec[complete_env_samples, , drop = FALSE]

env_numeric <- env_numeric[complete_env_samples, , drop = FALSE]

groupsize <- groupsize[groupsize$sample %in% complete_env_samples, , drop = FALSE]

run_mantel_by_sample_group <- function(group_name, sample_ids) {
  sample_ids <- intersect(sample_ids, rownames(spec))
  if (length(sample_ids) < 4) {
    stop(paste0(group_name, " 组样本数少于 4，无法稳定进行 Mantel test。"))
  }
  spec_sub <- spec[sample_ids, , drop = FALSE]
  env_sub <- env_numeric[sample_ids, , drop = FALSE]
  spec_sub <- spec_sub[, colSums(spec_sub, na.rm = TRUE) > 0, drop = FALSE]
  mantel_test(spec_sub, env_sub, mantel_fun = "mantel", spec_dist = "bray", env_dist = "euclidean",
    permutations = 9999) %>%
    mutate(spec = group_name)
}

mantel <- bind_rows(run_mantel_by_sample_group("all", rownames(spec)), run_mantel_by_sample_group("small",
  groupsize$sample[groupsize$group == "small"]), run_mantel_by_sample_group("large", groupsize$sample[groupsize$group ==
  "large"])) %>%
  mutate(spec = factor(spec, levels = c("all", "small", "large"), labels = c("All", "USF", "ConF"))) %>%
  group_by(spec) %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH"), r_abs = abs(r), r_level = cut(r_abs, breaks = c(-Inf,
    0.1, 0.2, 0.3, 0.4, Inf), labels = c("|r| < 0.10", "0.10 <= |r| < 0.20", "0.20 <= |r| < 0.30",
    "0.30 <= |r| < 0.40", "|r| >= 0.40"), right = FALSE), p_level = cut(p_adj_BH, breaks = c(-Inf,
    0.005, 0.01, 0.05, 0.1, Inf), labels = c("p < 0.005", "0.005 <= p < 0.010", "0.010 <= p < 0.050",
    "0.050 <= p < 0.100", "p >= 0.100"), right = FALSE)) %>%
  ungroup()

write.csv(mantel, file.path(output_dir, "mantel_all_small_large.csv"), quote = FALSE, row.names = FALSE)

write.csv(dplyr::filter(mantel, env %in% added_env_vars), file.path(output_dir, "mantel_added_NFR_HRT_AIRT.csv"),
  quote = FALSE, row.names = FALSE)

env_cor <- rcorr_pearson(env_numeric[, plot_env_vars, drop = FALSE])

write.csv(env_cor$r, file.path(output_dir, "environment_pearson_r.csv"))

write.csv(env_cor$P, file.path(output_dir, "environment_pearson_p.csv"))

write.csv(data.frame(community_input = normalizePath(tpm_file, winslash = "/"), rarefaction_depth = unique(rarefied_depth),
  n_samples = ncol(tpm), community_distance = "Bray-Curtis on rarefied counts", environmental_variables = paste(colnames(env_numeric),
    collapse = "; "), plot_environmental_variables = paste(plot_env_vars, collapse = "; "), original_environmental_variables = paste(original_env_vars,
    collapse = "; "), added_environmental_variables = paste(added_env_vars, collapse = "; "), excluded_new_hydrology = "MWL; RSV; Qin; Qout; Qnet",
  excluded_metadata_or_legacy = "Q; Longitude; Latitude; MonthGroup; SizeGroup", mantel_permutations = 9999),
  file.path(output_dir, "mantel_input_audit.csv"), row.names = FALSE)

gg <- qcorrplot(env_cor, type = "upper", diag = FALSE, grid_size = 0.15, grid_col = "grey50") + geom_square(aes(fill = r),
  colour = "black", size = 0.2) + scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
  midpoint = 0, limits = c(-1, 1), name = "Pearson's r") + geom_couple(aes(colour = p_level, size = r_level),
  data = mantel, curvature = 0.15, nudge_x = 0.6, offset_x = c(All = 3.8, USF = 3.8 - 0.275 * ncol(env_numeric),
    ConF = 3.8 - 0.55 * ncol(env_numeric)), offset_y = c(All = 3.9, USF = 0.1, ConF = -3.7), node.color = "grey20",
  node.fill = "grey80", alpha = 0.88, drop = TRUE) + scale_size_manual(values = c(0.25, 0.55, 0.9,
  1.35, 1.9), drop = FALSE, name = "Mantel's |r|") + scale_colour_manual(values = c(`p < 0.005` = "#6A3D9A",
  `0.005 <= p < 0.010` = "#0072B2", `0.010 <= p < 0.050` = "#009E73", `0.050 <= p < 0.100` = "#E69F00",
  `p >= 0.100` = "#C7C7C7"), drop = FALSE, name = "Mantel's BH-adjusted p") + guides(colour = guide_legend(order = 1,
  title.position = "top", direction = "vertical", ncol = 1, override.aes = list(size = 3, alpha = 1)),
  size = guide_legend(order = 2, title.position = "top", direction = "vertical", ncol = 1, override.aes = list(colour = "grey35",
    alpha = 1)), fill = guide_colorbar(order = 3, title.position = "top", direction = "vertical",
    barwidth = grid::unit(0.38, "cm"), barheight = grid::unit(3.2, "cm"))) + geom_mark(only_mark = TRUE,
  size = 3, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, colour = "white") + coord_fixed(clip = "off") +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, color = "black", size = 11.5),
    axis.text.y = element_text(color = "black", size = 11.5), legend.title = element_text(size = 11),
    legend.text = element_text(size = 10.5), legend.position = "bottom", legend.justification = "right",
    legend.location = "plot", legend.box = "horizontal", legend.box.just = "right", legend.box.margin = margin(t = 48,
      r = 0, b = 0, l = 0), legend.spacing.x = grid::unit(0.25, "cm"), legend.spacing.y = grid::unit(0.1,
      "cm"), plot.margin = margin(t = 25, r = 25, b = 12, l = 55))

ggsave(file.path(output_dir, "mantel_all_small_large.png"), plot = gg, width = 10, height = 10.2, dpi = 600)

ggsave(file.path(output_dir, "mantel_all_small_large.pdf"), plot = gg, width = 10, height = 10.2)

ggsave(file.path(output_dir, "mantel_all_small_large.svg"), plot = gg, width = 10, height = 10.2)

print(gg)
