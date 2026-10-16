options(stringsAsFactors = FALSE, scipen = 999)

try(Sys.setlocale("LC_ALL", ""), silent = TRUE)

ROOT_DIR <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT_DIR <- file.path(ROOT_DIR, "01_分析项目", "36_TITAN2致病菌ASV水平_水位")

TABLE_DIR <- file.path(PROJECT_DIR, "结果表格", "WL", "ASV")

FIGURE_DIR <- file.path(PROJECT_DIR, "结果图片", "WL", "ASV", "01_community_sumz_ASV")

read_csv <- function(name) {
  path <- file.path(TABLE_DIR, name)
  if (!file.exists(path))
    stop("缺少正式结果文件：", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

checks <- list()

add_check <- function(check, pass, expected, observed, detail_cn = "") {
  checks[[length(checks) + 1L]] <<- data.frame(check = check, pass = isTRUE(pass), expected = as.character(expected),
    observed = as.character(observed), detail_cn = detail_cn, stringsAsFactors = FALSE)
}

config <- read_csv("04_analysis_config.csv")

cfg <- stats::setNames(as.character(config$value), config$parameter)

metadata <- read_csv("01_sample_alignment_depth_and_ASV_audit.csv")

water_rep <- read_csv("02_environment_value_replication.csv")

screening <- read_csv("05_all_candidate_ASV_input_screening.csv")

id_map <- read_csv("05a_ASV_internal_ID_map.csv")

analyzed <- read_csv("10_all_analyzed_ASV_TITAN2_results.csv")

robust <- read_csv("11_robust_indicator_ASVs.csv")

negative <- read_csv("12_robust_negative_ASVs.csv")

positive <- read_csv("13_robust_positive_ASVs.csv")

profile <- read_csv("15_community_sumz_profile_by_split.csv")

community <- read_csv("16_community_thresholds.csv")

bootstrap_taxa <- read_csv("17_ASV_bootstrap_metrics_long.csv")

bootstrap_community <- read_csv("18_community_bootstrap_thresholds.csv")

ranking <- read_csv("19_robust_ASV_response_strength_ranking.csv")

add_check("formal_run_mode", cfg[["run_mode"]] == "formal", "formal", cfg[["run_mode"]], "拒绝把smoke结果当正式结果")

add_check("num_permutations", as.integer(cfg[["num_permutations"]]) == 500L, 500, cfg[["num_permutations"]])

add_check("n_bootstraps", as.integer(cfg[["n_bootstraps"]]) == 500L, 500, cfg[["n_bootstraps"]])

add_check("taxon_level", cfg[["taxon_level"]] == "ASV", "ASV", cfg[["taxon_level"]])

add_check("taxonomy_annotation_not_used", toupper(cfg[["taxonomy_annotation_used"]]) == "FALSE", FALSE,
  cfg[["taxonomy_annotation_used"]], "未读取、未筛选、未聚合任何分类注释")

add_check("sample_count", nrow(metadata) == 120L, 120, nrow(metadata))

add_check("full_rarefaction_depth", all(as.numeric(metadata$full_community_rarefied_reads) == 21000),
  "all samples = 21000", paste(range(metadata$full_community_rarefied_reads), collapse = "-"))

add_check("unique_water_levels", nrow(water_rep) == 15L, 15, nrow(water_rep))

add_check("replicates_per_water_level", all(as.integer(water_rep$Freq) == 8L), "all levels = 8", paste(range(as.integer(water_rep$Freq)),
  collapse = "-"))

add_check("candidate_ASV_screening_pool", nrow(screening) == 2755L, 2755, nrow(screening))

included <- as.logical(screening$included_in_TITAN2)

add_check("ASV_entering_TITAN2", sum(included) == 476L, 476, sum(included))

add_check("excluded_all_zero", sum(screening$occurrence == 0L) == 258L, 258, sum(screening$occurrence ==
  0L))

add_check("excluded_occurrence_1_to_4", sum(screening$occurrence > 0L & screening$occurrence < 5L) ==
  2021L, 2021, sum(screening$occurrence > 0L & screening$occurrence < 5L))

add_check("ASV_internal_ID_map_unique", nrow(id_map) == 2755L && !anyDuplicated(id_map$taxon_id) && !anyDuplicated(id_map$ASV_ID),
  "2755 one-to-one IDs", paste(nrow(id_map), anyDuplicated(id_map$taxon_id), anyDuplicated(id_map$ASV_ID),
    sep = ";"))

add_check("analyzed_result_rows", nrow(analyzed) == 476L, 476, nrow(analyzed))

screen_index <- match(analyzed$taxon_id, screening$taxon_id)

add_check("analyzed_IDs_equal_screened_IDs", !anyNA(screen_index) && setequal(analyzed$taxon_id, screening$taxon_id[included]),
  "same 476 IDs", paste("unmatched", sum(is.na(screen_index))))

add_check("occurrence_preserved", !anyNA(screen_index) && all(as.numeric(analyzed$occurrence) == as.numeric(screening$occurrence[screen_index])),
  "TITAN2 frequency = screening occurrence", if (anyNA(screen_index)) "unmatched IDs" else sum(as.numeric(analyzed$occurrence) !=
    as.numeric(screening$occurrence[screen_index])))

robust_flag <- as.logical(analyzed$robust_indicator)

purity_pass <- as.numeric(analyzed$purity) >= 0.95

reliability_pass <- as.numeric(analyzed$reliability) >= 0.95

add_check("robust_filter_definition", all(robust_flag == (purity_pass & reliability_pass)), "robust iff purity>=0.95 and reliability>=0.95",
  sum(robust_flag != (purity_pass & reliability_pass)))

add_check("robust_ID_set", setequal(robust$taxon_id, analyzed$taxon_id[robust_flag]), "same robust IDs",
  paste(nrow(robust), sum(robust_flag), sep = ";"))

add_check("robust_negative_ID_set", setequal(negative$taxon_id, robust$taxon_id[robust$robust_direction ==
  "negative_z-"]), "same negative IDs", nrow(negative))

add_check("robust_positive_ID_set", setequal(positive$taxon_id, robust$taxon_id[robust$robust_direction ==
  "positive_z+"]), "same positive IDs", nrow(positive))

add_check("taxon_bootstrap_rows", nrow(bootstrap_taxa) == 476L * 500L, 476L * 500L, nrow(bootstrap_taxa))

add_check("taxon_bootstrap_replicates", length(unique(bootstrap_taxa$bootstrap_replicate)) == 500L, 500,
  length(unique(bootstrap_taxa$bootstrap_replicate)))

add_check("community_bootstrap_rows", nrow(bootstrap_community) == 500L, 500, nrow(bootstrap_community))

add_check("strength_ranking_rows", nrow(ranking) == nrow(robust) && all(sort(as.integer(ranking$response_strength_rank)) ==
  seq_len(nrow(robust))), nrow(robust), nrow(ranking))

filtered <- community[community$threshold_type %in% c("fsumz-", "fsumz+") & community$taxa_scope == "robust_filtered",
  , drop = FALSE]

add_check("filtered_community_threshold_rows", nrow(filtered) == 2L && !anyDuplicated(filtered$direction),
  "one row per direction", nrow(filtered))

for (direction_i in c("negative_z-", "positive_z+")) {
  robust_n_i <- sum(robust$robust_direction == direction_i)
  threshold_i <- filtered[filtered$direction == direction_i, , drop = FALSE]
  add_check(paste0("robust_count_in_threshold_", direction_i), nrow(threshold_i) == 1L && as.integer(threshold_i$robust_taxa_n) ==
    robust_n_i, robust_n_i, if (nrow(threshold_i) == 1L)
    threshold_i$robust_taxa_n
  else "missing")
  if (nrow(threshold_i) == 1L && robust_n_i > 0L) {
    profile_column <- if (direction_i == "negative_z-") {
      "sumz_negative_robust_filtered"
    }
    else {
      "sumz_positive_robust_filtered"
    }
    raw_y <- as.numeric(profile[[profile_column]])
    peak_x <- as.numeric(profile$change_point)[abs(raw_y - max(raw_y, na.rm = TRUE)) < 1e-10]
    add_check(paste0("observed_filtered_peak_", direction_i), any(abs(peak_x - as.numeric(threshold_i$change_point_observed)) <
      1e-08), paste(peak_x, collapse = ";"), threshold_i$change_point_observed, "观察阈值必须来自原始filtered sum(z)峰，不能来自平滑曲线")
  }
}

figure_prefix <- file.path(FIGURE_DIR, "TITAN2_WL_ASV_community_sumz_smooth")

figure_files <- paste0(figure_prefix, c(".png", ".tiff", ".svg", ".pdf"))

figure_ok <- all(file.exists(figure_files)) && all(file.info(figure_files)$size > 0)

add_check("first_figure_four_formats", figure_ok, "PNG/TIFF/SVG/PDF all non-empty", paste(basename(figure_files[file.exists(figure_files)]),
  collapse = ";"))

qa <- do.call(rbind, checks)

write_csv(qa, file.path(TABLE_DIR, "22_formal_QA_checks.csv"))

key_results <- data.frame(metric = c("candidate_ASV_screening_pool", "ASV_entering_TITAN2", "screened_reads_percent_of_candidate_reads",
  "robust_ASV_total", "robust_negative_ASV", "robust_positive_ASV", "negative_observed_community_threshold_m",
  "negative_bootstrap_median_threshold_m", "negative_bootstrap_q05_threshold_m", "negative_bootstrap_q95_threshold_m",
  "positive_observed_community_threshold_m", "positive_bootstrap_median_threshold_m", "positive_bootstrap_q05_threshold_m",
  "positive_bootstrap_q95_threshold_m"), value = c(nrow(screening), sum(included), as.numeric(cfg[["screened_reads_percent_of_candidate_reads"]]),
  nrow(robust), nrow(negative), nrow(positive), filtered$change_point_observed[filtered$direction ==
    "negative_z-"], filtered$change_point_median[filtered$direction == "negative_z-"], filtered$change_point_q05[filtered$direction ==
    "negative_z-"], filtered$change_point_q95[filtered$direction == "negative_z-"], filtered$change_point_observed[filtered$direction ==
    "positive_z+"], filtered$change_point_median[filtered$direction == "positive_z+"], filtered$change_point_q05[filtered$direction ==
    "positive_z+"], filtered$change_point_q95[filtered$direction == "positive_z+"]), stringsAsFactors = FALSE)

write_csv(key_results, file.path(TABLE_DIR, "23_key_results_summary.csv"))

guide <- read_csv("00_output_table_guide.csv")

guide_extra <- data.frame(file = c("22_formal_QA_checks.csv", "23_key_results_summary.csv"), purpose_cn = c("正式运行、筛选、ID、稳健判定、bootstrap、总体阈值和四格式图件的逐项QA",
  "最重要的ASV数量、reads覆盖率、稳健ASV数及正负总体阈值汇总"), stringsAsFactors = FALSE)

guide <- rbind(guide[!guide$file %in% guide_extra$file, , drop = FALSE], guide_extra)

write_csv(guide, file.path(TABLE_DIR, "00_output_table_guide.csv"))

if (any(!qa$pass)) {
  stop("正式结果QA未全部通过；请查看22_formal_QA_checks.csv。")
}

message("正式结果QA全部通过：", nrow(qa), "项；稳健ASV ", nrow(robust), "（z- ", nrow(negative),
  "，z+ ", nrow(positive), "）。")
