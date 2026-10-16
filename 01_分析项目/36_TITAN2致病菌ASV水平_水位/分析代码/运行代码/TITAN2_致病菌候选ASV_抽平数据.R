options(stringsAsFactors = FALSE, scipen = 999)

try(Sys.setlocale("LC_ALL", ""), silent = TRUE)

ROOT_DIR <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT_DIR <- file.path(ROOT_DIR, "01_分析项目", "36_TITAN2致病菌ASV水平_水位")

args <- commandArgs(trailingOnly = TRUE)

ENV_VAR <- if (length(args) >= 1L) args[[1L]] else "WL"

RUN_MODE <- if (length(args) >= 2L) args[[2L]] else "formal"

ALLOWED_ENV <- c("WL", "TN", "Chla", "WT", "DO")

if (!ENV_VAR %in% ALLOWED_ENV) {
  stop("环境变量只允许：", paste(ALLOWED_ENV, collapse = ", "))
}

if (!RUN_MODE %in% c("formal", "smoke")) {
  stop("运行模式只允许 formal 或 smoke。")
}

ENV_LABELS <- c(WL = "Water level (m)", TN = "Total nitrogen", Chla = "Chlorophyll a", WT = "Water temperature",
  DO = "Dissolved oxygen")

ENV_LABELS_CN <- c(WL = "水位", TN = "总氮", Chla = "叶绿素a", WT = "水温", DO = "溶解氧")

ENV_LABEL <- unname(ENV_LABELS[[ENV_VAR]])

RUN_TAG <- if (RUN_MODE == "formal") "ASV" else "ASV_smoke"

TABLE_DIR <- file.path(PROJECT_DIR, "结果表格", ENV_VAR, RUN_TAG)

FIGURE_DIR <- file.path(PROJECT_DIR, "结果图片", ENV_VAR, RUN_TAG)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

FULL_RAREFIED_FILE <- file.path(ROOT_DIR, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv")

PATHOGEN_RAREFIED_FILE <- file.path(ROOT_DIR, "00_输入文件", "05_致病菌候选", "完整群落抽平21000后候选ASV子表.csv")

ENV_FILE <- file.path(ROOT_DIR, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv")

RAREFACTION_DEPTH <- 21000L

MIN_OCCURRENCES <- 5L

MIN_SPLIT <- 5L

PURITY_CUTOFF <- 0.95

RELIABILITY_CUTOFF <- 0.95

IMAX <- FALSE

IV_TOTAL <- FALSE

MEMORY_MODE <- TRUE

RANDOM_SEED <- 20260822L

MODE_SETTINGS <- list(formal = list(num_perm = 500L, n_boot = 500L), smoke = list(num_perm = 25L, n_boot = 20L))

NUM_PERM <- MODE_SETTINGS[[RUN_MODE]]$num_perm

N_BOOT <- MODE_SETTINGS[[RUN_MODE]]$n_boot

MAX_CPUS <- suppressWarnings(as.integer(Sys.getenv("TITAN_MAX_CPUS", "4")))

if (is.na(MAX_CPUS) || MAX_CPUS < 1L) MAX_CPUS <- 1L

.libPaths(unique(c(file.path(ROOT_DIR, "99_运行环境", "R_library"), .libPaths())))

required_packages <- c("TITAN2", "dplyr", "tidyr", "cli")

available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (any(!available)) {
  stop("缺少 R 包：", paste(required_packages[!available], collapse = ", "))
}

suppressPackageStartupMessages(library(TITAN2))

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

md5 <- function(path) unname(tools::md5sum(path))

run_log <- file.path(TABLE_DIR, "run.log")

if (file.exists(run_log)) file.remove(run_log)

write_log <- function(...) {
  msg <- paste0(...)
  message(msg)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n", file = run_log, append = TRUE, sep = "")
}

detect_physical_cores <- function() {
  x <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(x) || x < 1L)
    x <- suppressWarnings(parallel::detectCores())
  if (is.na(x) || x < 1L)
    x <- 1L
  as.integer(x)
}

N_CPUS <- max(1L, min(MAX_CPUS, detect_physical_cores(), N_BOOT))

patch_titan2_windows_parallel <- function() {
  if (.Platform$OS.type != "windows")
    return(invisible(FALSE))
  patched_boot_titan <- function(env, taxa, ivTot = ivTot, boot = boot, ncpus = ncpus, nBoot = nBoot,
    minSplt = minSplt, nPerm = 250, memory = memory, imax = imax, numUnit = numUnit) {
    if (ncpus > 1) {
      cli::cli_alert_info("Bootstrap resampling in parallel using {ncpus} CPUs.")
      cl <- parallel::makeCluster(rep("localhost", ncpus))
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      worker_libpaths <- .libPaths()
      parallel::clusterCall(cl, function(paths) {
        try(Sys.setlocale("LC_ALL", ""), silent = TRUE)
        .libPaths(paths)
        suppressPackageStartupMessages(library(TITAN2))
        NULL
      }, worker_libpaths)
      parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1L))
      b_seq <- parallel::clusterSplit(cl, seq_len(nBoot))
      out <- parallel::clusterApply(cl, b_seq, TITAN2:::tboot, env = env, taxa = taxa, ivTot = ivTot,
        minSplt = minSplt, nPerm = nPerm, memory = memory, imax = imax)
    }
    else {
      b_seq <- 0
      out <- vector("list", nBoot)
      for (i in seq_len(nBoot)) {
        out[[i]] <- TITAN2:::tboot(i, env = env, taxa = taxa, ivTot = ivTot, minSplt = minSplt,
          nPerm = nPerm, memory = memory, imax = imax)
      }
    }
    list(b_seq, out)
  }
  assignInNamespace("boot.titan", patched_boot_titan, ns = "TITAN2")
  invisible(TRUE)
}

patch_titan2_windows_parallel()

assert_count_matrix <- function(x, label) {
  if (!is.numeric(x))
    stop(label, " 不是数值矩阵。")
  if (any(!is.finite(x)))
    stop(label, " 含 NA/Inf。")
  if (any(x < 0))
    stop(label, " 含负数。")
  if (any(abs(x - round(x)) > 1e-08))
    stop(label, " 含非整数计数。")
  invisible(TRUE)
}

extract_taxon_results <- function(result, map) {
  s <- as.data.frame(result$sppmax, check.names = FALSE)
  s$taxon_id <- rownames(s)
  mi <- match(s$taxon_id, map$taxon_id)
  if (anyNA(mi))
    stop("TITAN2结果中的ASV ID无法映射回输入。")
  data.frame(taxon_level = "ASV", taxonomy_annotation_used = FALSE, pathogen_scope = "candidate_pathogen_ASVs",
    taxon_id = s$taxon_id, taxon_name = map$taxon_name[mi], ASV_ID = map$ASV_ID[mi], direction_code = s$maxgrp,
    direction = ifelse(s$maxgrp == 1, "negative_z-", "positive_z+"), robust_indicator = s$filter %in%
      c(1, 2), robust_direction = ifelse(s$filter == 1, "negative_z-", ifelse(s$filter == 2, "positive_z+",
      "not_robust")), occurrence = s$freq, change_point_observed = s$zenv.cp, change_point_q05 = s[["5%"]],
    change_point_q10 = s[["10%"]], change_point_median = s[["50%"]], change_point_q90 = s[["90%"]],
    change_point_q95 = s[["95%"]], change_point_uncertainty_width_90 = s[["95%"]] - s[["5%"]], IndVal_observed = s$IndVal,
    permutation_probability_observed = s$obsiv.prob, zscore_observed = s$zscore, zscore_bootstrap_median = s$z.median,
    purity = s$purity, reliability = s$reliability, purity_pass = s$purity >= PURITY_CUTOFF, reliability_pass = s$reliability >=
      RELIABILITY_CUTOFF, filter_code = s$filter, stringsAsFactors = FALSE)
}

extract_candidate_metrics <- function(result, map) {
  scores <- result$ivzScores
  n_taxa <- nrow(map)
  n_split <- ncol(scores)
  srt <- as.numeric(result$srtEnv)
  split_index <- seq_len(n_split)
  left_index <- MIN_SPLIT + split_index - 1L
  split_info <- data.frame(split_index = split_index, environment_left = srt[left_index], environment_right = srt[left_index +
    1L], change_point = (srt[left_index] + srt[left_index + 1L])/2, n_samples_low_side = left_index,
    n_samples_high_side = length(srt) - left_index, stringsAsFactors = FALSE)
  ti <- rep(seq_len(n_taxa), each = n_split)
  si <- rep(seq_len(n_split), times = n_taxa)
  direction <- as.vector(t(scores[seq_len(n_taxa), , drop = FALSE]))
  data.frame(taxon_level = "ASV", taxon_id = map$taxon_id[ti], taxon_name = map$taxon_name[ti], ASV_ID = map$ASV_ID[ti],
    split_info[si, , drop = FALSE], direction_code = direction, direction = ifelse(direction == 1,
      "negative_z-", "positive_z+"), zscore = as.vector(t(scores[n_taxa + seq_len(n_taxa), , drop = FALSE])),
    IndVal = as.vector(t(scores[2L * n_taxa + seq_len(n_taxa), , drop = FALSE])), permutation_probability = as.vector(t(scores[3L *
      n_taxa + seq_len(n_taxa), , drop = FALSE])), stringsAsFactors = FALSE, row.names = NULL)
}

extract_bootstrap_taxa <- function(result, map, taxon_result) {
  a <- result$metricArray
  n_taxa <- dim(a)[1L]
  n_boot <- dim(a)[3L]
  if (n_taxa != nrow(map))
    stop("bootstrap ASV维度与输入不一致。")
  ti <- rep(seq_len(n_taxa), each = n_boot)
  bi <- rep(seq_len(n_boot), times = n_taxa)
  metric <- function(k) as.vector(t(a[, k, , drop = FALSE][, 1L, ]))
  direction <- metric(1L)
  ri <- match(map$taxon_id, taxon_result$taxon_id)
  if (anyNA(ri))
    stop("bootstrap ASV无法映射回TITAN2结果。")
  data.frame(taxon_level = "ASV", taxon_id = map$taxon_id[ti], taxon_name = map$taxon_name[ti], ASV_ID = map$ASV_ID[ti],
    robust_indicator = rep(taxon_result$robust_indicator[ri], each = n_boot), robust_direction = rep(taxon_result$robust_direction[ri],
      each = n_boot), bootstrap_replicate = bi, direction_code = direction, direction = ifelse(direction ==
      1, "negative_z-", "positive_z+"), change_point = metric(2L), zscore = metric(3L), permutation_probability = metric(4L),
    stringsAsFactors = FALSE)
}

extract_community <- function(result, robust) {
  cp <- as.data.frame(result$sumz.cp, check.names = FALSE)
  cp$threshold_type <- rownames(cp)
  rownames(cp) <- NULL
  data.frame(taxon_level = "ASV", threshold_type = cp$threshold_type, taxa_scope = ifelse(grepl("^f",
    cp$threshold_type), "robust_filtered", "all_taxa"), direction = ifelse(grepl("-$", cp$threshold_type),
    "negative_z-", "positive_z+"), change_point_observed = cp$cp, change_point_q05 = cp[["0.05"]],
    change_point_q10 = cp[["0.10"]], change_point_median = cp[["0.50"]], change_point_q90 = cp[["0.90"]],
    change_point_q95 = cp[["0.95"]], robust_taxa_n = ifelse(grepl("-$", cp$threshold_type), sum(robust$robust_direction ==
      "negative_z-"), sum(robust$robust_direction == "positive_z+")), stringsAsFactors = FALSE)
}

extract_sumz_profile <- function(result) {
  n <- nrow(result$ivz)
  srt <- as.numeric(result$srtEnv)
  j <- seq_len(n)
  left <- MIN_SPLIT + j - 1L
  data.frame(taxon_level = "ASV", split_index = j, environment_left = srt[left], environment_right = srt[left +
    1L], change_point = (srt[left] + srt[left + 1L])/2, n_samples_low_side = left, n_samples_high_side = length(srt) -
    left, sumz_negative_all_taxa = result$ivz[, 1L], sumz_positive_all_taxa = result$ivz[, 2L], sumz_negative_robust_filtered = result$ivz.f[,
    1L], sumz_positive_robust_filtered = result$ivz.f[, 2L], stringsAsFactors = FALSE)
}

extract_bootstrap_community <- function(result) {
  data.frame(taxon_level = "ASV", bootstrap_replicate = seq_len(nrow(result$maxSumz)), sumz_negative_all_taxa_cp = result$maxSumz[,
    1L], sumz_positive_all_taxa_cp = result$maxSumz[, 2L], fsumz_negative_robust_cp = result$maxFsumz[,
    1L], fsumz_positive_robust_cp = result$maxFsumz[, 2L], stringsAsFactors = FALSE)
}

write_log("读取并核验完整抽平表、候选致病菌ASV子表和环境表；不读取任何分类注释。")

full <- read.csv(FULL_RAREFIED_FILE, check.names = FALSE)

candidate <- read.csv(PATHOGEN_RAREFIED_FILE, check.names = FALSE)

names(full)[1L] <- names(candidate)[1L] <- "ASV_ID"

if (anyNA(full$ASV_ID) || anyNA(candidate$ASV_ID) || any(trimws(full$ASV_ID) == "") || any(trimws(candidate$ASV_ID) ==
  "")) {
  stop("ASV ID存在缺失或空值。")
}

if (anyDuplicated(full$ASV_ID) || anyDuplicated(candidate$ASV_ID)) {
  stop("ASV ID重复。")
}

sample_ids <- names(full)[-1L]

if (anyDuplicated(sample_ids) || any(sample_ids == "")) {
  stop("样本列名重复或为空。")
}

if (!identical(names(candidate)[-1L], sample_ids)) {
  stop("候选子表与完整抽平表的样本列不一致。")
}

full_matrix <- as.matrix(full[, -1L, drop = FALSE])

candidate_matrix <- as.matrix(candidate[, -1L, drop = FALSE])

storage.mode(full_matrix) <- storage.mode(candidate_matrix) <- "numeric"

assert_count_matrix(full_matrix, "完整抽平表")

assert_count_matrix(candidate_matrix, "候选致病菌ASV子表")

if (any(colSums(full_matrix) != RAREFACTION_DEPTH)) {
  stop("完整群落不是严格抽平至21,000 reads。")
}

full_idx <- match(candidate$ASV_ID, full$ASV_ID)

if (anyNA(full_idx) || !isTRUE(all.equal(candidate_matrix, full_matrix[full_idx, , drop = FALSE], tolerance = 0))) {
  stop("候选致病菌ASV子表不能由完整抽平表精确复现。")
}

safe_taxon_id <- sprintf("A%04d", seq_len(nrow(candidate)))

rownames(candidate_matrix) <- safe_taxon_id

env <- read.csv(ENV_FILE, check.names = FALSE)

if (!all(c("sample", ENV_VAR) %in% names(env))) {
  stop("环境表缺少 sample 或 ", ENV_VAR)
}

if (anyDuplicated(env$sample)) stop("环境表样本ID重复。")

env_index <- match(sample_ids, env$sample)

if (anyNA(env_index)) stop("候选ASV表中的样本未全部匹配环境表。")

environment_value <- suppressWarnings(as.numeric(env[[ENV_VAR]][env_index]))

if (any(!is.finite(environment_value))) {
  stop("120个样本并非全部具有有效环境值。")
}

occurrence <- rowSums(candidate_matrix > 0)

variable_abundance <- apply(candidate_matrix, 1L, function(x) length(unique(x)) > 1L)

included <- occurrence >= MIN_OCCURRENCES & variable_abundance

total_by_asv <- rowSums(candidate_matrix)

screened_reads <- sum(total_by_asv[included])

all_candidate_reads <- sum(total_by_asv)

screening <- data.frame(taxon_level = "ASV", taxonomy_annotation_used = FALSE, pathogen_scope = "candidate_pathogen_ASVs",
  taxon_id = safe_taxon_id, taxon_name = candidate$ASV_ID, ASV_ID = candidate$ASV_ID, occurrence = occurrence,
  prevalence = occurrence/ncol(candidate_matrix), total_rarefied_count = total_by_asv, mean_rarefied_count = rowMeans(candidate_matrix),
  sd_rarefied_count = apply(candidate_matrix, 1L, stats::sd), median_rarefied_count = apply(candidate_matrix,
    1L, stats::median), min_rarefied_count = apply(candidate_matrix, 1L, min), max_rarefied_count = apply(candidate_matrix,
    1L, max), variable_abundance = variable_abundance, required_occurrence = MIN_OCCURRENCES, included_in_TITAN2 = included,
  exclusion_reason = ifelse(included, "", ifelse(occurrence == 0L, "all_zero", ifelse(occurrence <
    MIN_OCCURRENCES, "occurrence_below_5_samples", "constant_abundance"))), stringsAsFactors = FALSE)

if (nrow(candidate) != 2755L || length(sample_ids) != 120L) {
  stop("候选ASV数或样本数与审计基线不一致。")
}

if (sum(included) != 476L || all_candidate_reads != 352425 || screened_reads != 294691) {
  stop("ASV筛选数或reads覆盖率与审计基线不一致。")
}

counts_sample_by_asv <- t(candidate_matrix)

counts_sample_by_asv <- counts_sample_by_asv[sample_ids, , drop = FALSE]

screened_counts <- counts_sample_by_asv[, included, drop = FALSE]

map_for_titan <- screening[included, c("taxon_id", "taxon_name", "ASV_ID"), drop = FALSE]

metadata <- data.frame(sample = sample_ids, environment_variable = ENV_VAR, environment_label_cn = unname(ENV_LABELS_CN[[ENV_VAR]]),
  environment_value = environment_value, full_community_rarefied_reads = as.numeric(colSums(full_matrix)),
  all_candidate_pathogen_ASV_reads = as.numeric(colSums(candidate_matrix)), screened_candidate_pathogen_ASV_reads = as.numeric(rowSums(screened_counts)),
  candidate_reads_percent_of_total_community = as.numeric(colSums(candidate_matrix))/RAREFACTION_DEPTH *
    100, screened_reads_percent_of_total_community = as.numeric(rowSums(screened_counts))/RAREFACTION_DEPTH *
    100, stringsAsFactors = FALSE)

write_csv(metadata, file.path(TABLE_DIR, "01_sample_alignment_depth_and_ASV_audit.csv"))

write_csv(as.data.frame(table(environment_value = metadata$environment_value)), file.path(TABLE_DIR,
  "02_environment_value_replication.csv"))

input_paths <- c(full_rarefied_ASV = FULL_RAREFIED_FILE, candidate_pathogen_rarefied_ASV = PATHOGEN_RAREFIED_FILE,
  environment = ENV_FILE)

write_csv(data.frame(role = names(input_paths), absolute_path = unname(input_paths), bytes = file.info(input_paths)$size,
  md5 = vapply(input_paths, md5, character(1)), stringsAsFactors = FALSE), file.path(TABLE_DIR, "03_input_file_provenance.csv"))

config <- data.frame(parameter = c("environment_variable", "environment_label_cn", "taxon_level", "taxonomy_annotation_used",
  "analysis_pool", "candidate_ASV_total", "screened_ASV_entering_TITAN2", "screening_rule", "abundance_input",
  "relative_abundance_transform", "rarefaction_depth_full_community", "candidate_reads_total", "screened_reads_total",
  "screened_reads_percent_of_candidate_reads", "sample_scope", "min_occurrences", "min_split", "num_permutations",
  "n_bootstraps", "purity_cutoff", "reliability_cutoff", "imax", "iv_total", "bootstrap_disk_memory_mode",
  "random_seed", "run_mode", "n_cpus", "TITAN2_version", "R_version"), value = c(ENV_VAR, unname(ENV_LABELS_CN[[ENV_VAR]]),
  "ASV", "FALSE", "all_2755_candidate_pathogen_ASVs_before_screening", nrow(candidate), sum(included),
  "occurrence_at_least_5_samples_and_abundance_not_constant", "rarefied_candidate_ASV_counts", "none",
  RAREFACTION_DEPTH, all_candidate_reads, screened_reads, sprintf("%.6f", 100 * screened_reads/all_candidate_reads),
  "all_120_samples", MIN_OCCURRENCES, MIN_SPLIT, NUM_PERM, N_BOOT, PURITY_CUTOFF, RELIABILITY_CUTOFF,
  IMAX, IV_TOTAL, MEMORY_MODE, RANDOM_SEED, RUN_MODE, N_CPUS, as.character(packageVersion("TITAN2")),
  R.version.string), stringsAsFactors = FALSE)

write_csv(config, file.path(TABLE_DIR, "04_analysis_config.csv"))

write_csv(screening, file.path(TABLE_DIR, "05_all_candidate_ASV_input_screening.csv"))

write_csv(screening[, c("taxon_id", "taxon_name", "ASV_ID", "taxonomy_annotation_used", "included_in_TITAN2")],
  file.path(TABLE_DIR, "05a_ASV_internal_ID_map.csv"))

n_sample <- nrow(metadata)

n_asv_all <- ncol(counts_sample_by_asv)

all_counts_wide <- cbind(metadata[, c("sample", "environment_value", "full_community_rarefied_reads",
  "all_candidate_pathogen_ASV_reads", "screened_candidate_pathogen_ASV_reads")], as.data.frame(counts_sample_by_asv,
  check.names = FALSE))

names(all_counts_wide)[(ncol(all_counts_wide) - n_asv_all + 1L):ncol(all_counts_wide)] <- screening$ASV_ID

write_csv(all_counts_wide, file.path(TABLE_DIR, "06_all_candidate_ASV_counts_sample_wide.csv"))

rm(all_counts_wide)

long <- data.frame(taxon_level = "ASV", taxonomy_annotation_used = FALSE, pathogen_scope = "candidate_pathogen_ASVs",
  sample = rep(metadata$sample, times = n_asv_all), environment_variable = ENV_VAR, environment_value = rep(metadata$environment_value,
    times = n_asv_all), taxon_id = rep(screening$taxon_id, each = n_sample), taxon_name = rep(screening$taxon_name,
    each = n_sample), ASV_ID = rep(screening$ASV_ID, each = n_sample), rarefied_count = as.vector(counts_sample_by_asv),
  relative_abundance_percent_of_total_community = as.vector(counts_sample_by_asv)/RAREFACTION_DEPTH *
    100, stringsAsFactors = FALSE)

long$present <- long$rarefied_count > 0

long$included_in_TITAN2 <- rep(included, each = n_sample)

write_csv(long, file.path(TABLE_DIR, "07_all_candidate_ASV_abundance_sample_long.csv"))

env_summary <- dplyr::summarise(dplyr::group_by(long, taxon_level, taxon_id, taxon_name, ASV_ID, included_in_TITAN2,
  environment_variable, environment_value), n_samples = dplyr::n(), occurrence = sum(present), prevalence = mean(present),
  mean_count = mean(rarefied_count), sd_count = stats::sd(rarefied_count), se_count = sd_count/sqrt(n_samples),
  median_count = stats::median(rarefied_count), min_count = min(rarefied_count), max_count = max(rarefied_count),
  mean_relative_abundance_percent_of_total_community = mean(relative_abundance_percent_of_total_community),
  sd_relative_abundance_percent_of_total_community = stats::sd(relative_abundance_percent_of_total_community),
  se_relative_abundance_percent_of_total_community = sd_relative_abundance_percent_of_total_community/sqrt(n_samples),
  median_relative_abundance_percent_of_total_community = stats::median(relative_abundance_percent_of_total_community),
  .groups = "drop")

write_csv(env_summary, file.path(TABLE_DIR, "08_all_candidate_ASV_summary_by_environment_value.csv"))

rm(long, env_summary)

screened_counts_wide <- cbind(metadata[, c("sample", "environment_value")], as.data.frame(screened_counts,
  check.names = FALSE))

names(screened_counts_wide)[3L:ncol(screened_counts_wide)] <- screening$ASV_ID[included]

write_csv(screened_counts_wide, file.path(TABLE_DIR, "09_screened_ASV_counts_sample_wide.csv"))

rm(screened_counts_wide, full, full_matrix, candidate_matrix)

invisible(gc(verbose = FALSE))

setwd(TABLE_DIR)

set.seed(RANDOM_SEED + 100L * match(ENV_VAR, ALLOWED_ENV) + ifelse(RUN_MODE == "smoke", 10000L, 0L))

write_log("开始 ", ENV_VAR, " 的致病菌候选ASV水平TITAN2；总候选ASV ", nrow(candidate),
  "；筛选后进入分析 ", sum(included), "；分类注释未读取。")

titan_result <- TITAN2::titan(env = metadata$environment_value, txa = screened_counts, minSplt = MIN_SPLIT,
  numPerm = NUM_PERM, boot = TRUE, nBoot = N_BOOT, imax = IMAX, ivTot = IV_TOTAL, pur.cut = PURITY_CUTOFF,
  rel.cut = RELIABILITY_CUTOFF, ncpus = N_CPUS, memory = MEMORY_MODE, messaging = TRUE)

saveRDS(titan_result, file.path(TABLE_DIR, "TITAN2_result_full.rds"), compress = "xz")

taxa_result <- extract_taxon_results(titan_result, map_for_titan)

screen_index <- match(taxa_result$taxon_id, screening$taxon_id)

if (anyNA(screen_index) || !all(taxa_result$occurrence == screening$occurrence[screen_index])) {
  stop("TITAN2输出的ASV occurrence与筛选表不一致。")
}

if (!all(taxa_result$robust_indicator == (taxa_result$purity_pass & taxa_result$reliability_pass))) {
  stop("TITAN2 filter与purity/reliability稳健性判定不一致。")
}

taxa_result <- taxa_result[order(!taxa_result$robust_indicator, taxa_result$robust_direction, taxa_result$change_point_median,
  taxa_result$ASV_ID), , drop = FALSE]

robust <- taxa_result[taxa_result$robust_indicator, , drop = FALSE]

negative <- robust[robust$robust_direction == "negative_z-", , drop = FALSE]

positive <- robust[robust$robust_direction == "positive_z+", , drop = FALSE]

candidate_metrics <- extract_candidate_metrics(titan_result, map_for_titan)

sumz_profile <- extract_sumz_profile(titan_result)

community <- extract_community(titan_result, robust)

bootstrap_taxa <- extract_bootstrap_taxa(titan_result, map_for_titan, taxa_result)

bootstrap_community <- extract_bootstrap_community(titan_result)

if (nrow(bootstrap_community) != N_BOOT || nrow(bootstrap_taxa) != nrow(map_for_titan) * N_BOOT) {
  stop("bootstrap输出行数与设定重复次数不一致。")
}

filtered_thresholds <- community[community$taxa_scope == "robust_filtered" & community$robust_taxa_n >
  0, , drop = FALSE]

for (i in seq_len(nrow(filtered_thresholds))) {
  direction_i <- filtered_thresholds$direction[i]
  profile_column <- if (direction_i == "negative_z-") {
    "sumz_negative_robust_filtered"
  }
  else {
    "sumz_positive_robust_filtered"
  }
  raw_y <- sumz_profile[[profile_column]]
  peak_x <- sumz_profile$change_point[abs(raw_y - max(raw_y, na.rm = TRUE)) < 1e-10]
  if (!any(abs(peak_x - filtered_thresholds$change_point_observed[i]) < 1e-08)) {
    stop("观察稳健总体阈值不是原始filtered sum(z)峰。")
  }
}

strength_ranking <- dplyr::mutate(dplyr::arrange(robust, dplyr::desc(zscore_bootstrap_median), permutation_probability_observed,
  dplyr::desc(purity), dplyr::desc(reliability), ASV_ID), response_strength_rank = dplyr::row_number(),
  selected_for_top15_summary = response_strength_rank <= 15L, ranking_primary_metric = "zscore_bootstrap_median_descending",
  ranking_tie_breakers = paste("permutation_probability_observed_ascending;", "purity_descending;reliability_descending;ASV_ID_alphabetical"))

top15 <- strength_ranking[strength_ranking$selected_for_top15_summary, , drop = FALSE]

write_csv(taxa_result, file.path(TABLE_DIR, "10_all_analyzed_ASV_TITAN2_results.csv"))

write_csv(robust, file.path(TABLE_DIR, "11_robust_indicator_ASVs.csv"))

write_csv(negative, file.path(TABLE_DIR, "12_robust_negative_ASVs.csv"))

write_csv(positive, file.path(TABLE_DIR, "13_robust_positive_ASVs.csv"))

write_csv(candidate_metrics, file.path(TABLE_DIR, "14_ASV_by_candidate_split_metrics_long.csv"))

write_csv(sumz_profile, file.path(TABLE_DIR, "15_community_sumz_profile_by_split.csv"))

write_csv(community, file.path(TABLE_DIR, "16_community_thresholds.csv"))

write_csv(bootstrap_taxa, file.path(TABLE_DIR, "17_ASV_bootstrap_metrics_long.csv"))

write_csv(bootstrap_community, file.path(TABLE_DIR, "18_community_bootstrap_thresholds.csv"))

write_csv(strength_ranking, file.path(TABLE_DIR, "19_robust_ASV_response_strength_ranking.csv"))

write_csv(top15, file.path(TABLE_DIR, "20_top15_most_pronounced_robust_ASVs.csv"))

summary_table <- data.frame(metric = c("candidate_ASV_total", "ASV_entering_TITAN2", "ASV_excluded_all_zero",
  "ASV_excluded_occurrence_1_to_4", "candidate_reads_total", "screened_reads_total", "screened_reads_percent_of_candidate_reads",
  "robust_ASV_total", "robust_negative_ASV", "robust_positive_ASV"), value = c(nrow(candidate), sum(included),
  sum(occurrence == 0L), sum(occurrence > 0L & occurrence < MIN_OCCURRENCES), all_candidate_reads,
  screened_reads, 100 * screened_reads/all_candidate_reads, nrow(robust), nrow(negative), nrow(positive)),
  stringsAsFactors = FALSE)

write_csv(summary_table, file.path(TABLE_DIR, "21_analysis_summary.csv"))

guide <- data.frame(file = c("01_sample_alignment_depth_and_ASV_audit.csv", "02_environment_value_replication.csv",
  "03_input_file_provenance.csv", "04_analysis_config.csv", "05_all_candidate_ASV_input_screening.csv",
  "05a_ASV_internal_ID_map.csv", "06_all_candidate_ASV_counts_sample_wide.csv", "07_all_candidate_ASV_abundance_sample_long.csv",
  "08_all_candidate_ASV_summary_by_environment_value.csv", "09_screened_ASV_counts_sample_wide.csv",
  "10_all_analyzed_ASV_TITAN2_results.csv", "11_robust_indicator_ASVs.csv", "12_robust_negative_ASVs.csv",
  "13_robust_positive_ASVs.csv", "14_ASV_by_candidate_split_metrics_long.csv", "15_community_sumz_profile_by_split.csv",
  "16_community_thresholds.csv", "17_ASV_bootstrap_metrics_long.csv", "18_community_bootstrap_thresholds.csv",
  "19_robust_ASV_response_strength_ranking.csv", "20_top15_most_pronounced_robust_ASVs.csv", "21_analysis_summary.csv",
  "TITAN2_result_full.rds"), purpose_cn = c("120个样本的环境匹配、完整群落抽平深度、筛选前后候选ASV reads",
  "各水位的重复样本数", "所有输入文件的绝对路径、大小和MD5", "全部分析参数；明确分类注释未使用",
  "全部2755个候选ASV的出现次数、丰度统计、是否入模和排除原因", "TITAN2安全内部编号与原始32位ASV哈希的一对一映射",
  "120个样本×全部2755个候选ASV的抽平计数宽表", "样本×全部候选ASV长表，含水位、计数、总群落相对百分比和入模状态",
  "每个水位值下全部候选ASV的计数、出现率及相对丰度汇总", "120个样本×筛选后476个入模ASV的抽平计数宽表",
  "全部476个入模ASV的TITAN2阈值、方向、z、P、purity和reliability", "purity与reliability均达到0.95的全部稳健ASV",
  "稳健负响应ASV", "稳健正响应ASV", "每个入模ASV×每个候选切点的方向、z、IndVal和置换概率",
  "第一张图四条community sum(z)曲线的原始候选切点数据", "全部入模ASV和稳健ASV的正负总体阈值及bootstrap分位数",
  "每个入模ASV×全部bootstrap的方向、阈值、z和概率", "每次bootstrap的四类总体阈值",
  "全部稳健ASV按bootstrap median z从高到低的响应强度排名", "响应最明显的前15个稳健ASV汇总（当前不用于第一张图筛选）",
  "筛选规模、reads覆盖率和稳健ASV数量摘要", "完整TITAN2对象，供复核或重新提取结果"),
  stringsAsFactors = FALSE)

write_csv(guide, file.path(TABLE_DIR, "00_output_table_guide.csv"))

capture.output(sessionInfo(), file = file.path(TABLE_DIR, "sessionInfo.txt"))

write_log("完成：候选ASV ", nrow(candidate), "；进入TITAN2 ", sum(included), "；稳健z- ",
  nrow(negative), "；稳健z+ ", nrow(positive), "；筛选后保留候选reads ", sprintf("%.3f%%",
    100 * screened_reads/all_candidate_reads), "。")
