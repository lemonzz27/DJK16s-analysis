suppressPackageStartupMessages({
  library(vegan)
  library(permute)
})

set.seed(20260814)

options(stringsAsFactors = FALSE, width = 140)

project_dir <- normalizePath(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

input_asv <- file.path(project_dir, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv")

input_size <- file.path(project_dir, "00_输入文件", "02_样本分组与站点", "孔径分组_groupsize.csv")

input_month <- file.path(project_dir, "00_输入文件", "02_样本分组与站点", "月份分组_groupmonth.csv")

input_env <- file.path(project_dir, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv")

analysis_dir <- file.path(project_dir, "01_分析项目", "33_VPA时间环境空间")

table_dir <- file.path(analysis_dir, "结果表格", "基础_注释无关", "整体_时间环境空间_VPA")

figure_dir <- file.path(analysis_dir, "结果图片", "基础_注释无关", "整体_时间环境空间_VPA")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("[1/7] 读取并核对输入数据")

asv_raw <- read.csv(input_asv, check.names = FALSE, stringsAsFactors = FALSE)

size_raw <- read.csv(input_size, check.names = FALSE, stringsAsFactors = FALSE)

month_raw <- read.csv(input_month, check.names = FALSE, stringsAsFactors = FALSE)

env_raw <- read.csv(input_env, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(names(asv_raw)[1] == "mags")

sample_ids <- names(asv_raw)[-1]

stopifnot(!anyDuplicated(sample_ids), !anyDuplicated(size_raw$sample), !anyDuplicated(month_raw$sample),
  !anyDuplicated(env_raw$sample))

match_required <- function(ids, table_ids, table_name) {
  idx <- match(ids, table_ids)
  if (anyNA(idx)) {
    stop(table_name, " 缺少样本：", paste(ids[is.na(idx)], collapse = ", "))
  }
  idx
}

idx_size <- match_required(sample_ids, size_raw$sample, "孔径分组表")

idx_month <- match_required(sample_ids, month_raw$sample, "月份分组表")

idx_env <- match_required(sample_ids, env_raw$sample, "环境表")

site_from_id <- function(x) {
  ifelse(startsWith(x, "LHKX"), "LHKX", sub("[0-9].*$", "", x))
}

meta <- data.frame(sample = sample_ids, size = tolower(size_raw$group[idx_size]), month = as.character(month_raw$group[idx_month]),
  site = site_from_id(sample_ids), stringsAsFactors = FALSE)

meta$size <- factor(meta$size, levels = c("small", "large"))

month_levels <- unique(meta$month[order(match(meta$month, unique(meta$month)))])

month_levels <- sort(unique(meta$month))

meta$month <- factor(meta$month, levels = month_levels)

meta$site <- factor(meta$site, levels = c("KX", "TC", "BS", "LHKX"))

if (anyNA(meta$size) || anyNA(meta$month) || anyNA(meta$site)) {
  stop("样本编号解析或分组表中存在无法识别的值。")
}

asv_num <- as.matrix(asv_raw[, -1, drop = FALSE])

storage.mode(asv_num) <- "double"

rownames(asv_num) <- asv_raw$mags

taxa_keep <- rowSums(asv_num, na.rm = TRUE) > 0.1

asv_num <- asv_num[taxa_keep, , drop = FALSE]

comm_all <- t(asv_num)

rownames(comm_all) <- sample_ids

if (anyNA(comm_all) || any(comm_all < 0)) stop("ASV 表含 NA 或负值。")

lib_sizes <- rowSums(comm_all)

if (any(lib_sizes <= 0)) stop("至少一个样本的总丰度为 0。")

env_aligned <- env_raw[idx_env, , drop = FALSE]

rownames(env_aligned) <- sample_ids

env_candidates <- setdiff(names(env_aligned), c("sample", "Longitude", "Latitude", "Q"))

env_numeric <- as.data.frame(lapply(env_aligned[, env_candidates, drop = FALSE], as.numeric), check.names = FALSE)

rownames(env_numeric) <- sample_ids

if (anyNA(env_numeric)) stop("候选环境变量转换为数值后出现 NA。")

nonconstant <- vapply(env_numeric, function(x) is.finite(sd(x)) && sd(x) > 0, logical(1))

constant_removed <- names(env_numeric)[!nonconstant]

env_numeric <- env_numeric[, nonconstant, drop = FALSE]

env_scaled_all <- as.data.frame(scale(env_numeric), check.names = FALSE)

message("[2/7] 仅按环境变量自身的 VIF 筛查共线性")

calc_vif <- function(dat) {
  if (ncol(dat) == 1L)
    return(setNames(1, names(dat)))
  ans <- vapply(seq_along(dat), function(j) {
    fit <- lm(dat[[j]] ~ ., data = dat[, -j, drop = FALSE])
    r2 <- summary(fit)$r.squared
    if (!is.finite(r2) || r2 >= 1)
      Inf
    else 1/(1 - r2)
  }, numeric(1))
  names(ans) <- names(dat)
  ans
}

vif_threshold <- 10

vif_work <- env_scaled_all

vif_history <- list()

removed_vif <- character()

iteration <- 1L

repeat {
  vifs <- calc_vif(vif_work)
  remove_now <- if (max(vifs) > vif_threshold)
    names(which.max(vifs))
  else NA_character_
  vif_history[[iteration]] <- data.frame(iteration = iteration, variable = names(vifs), VIF = as.numeric(vifs),
    action = ifelse(names(vifs) == remove_now, "removed", ifelse(is.na(remove_now), "retained_final",
      "retained")), stringsAsFactors = FALSE)
  if (is.na(remove_now) || ncol(vif_work) <= 1L)
    break
  removed_vif <- c(removed_vif, remove_now)
  vif_work[[remove_now]] <- NULL
  iteration <- iteration + 1L
}

vif_history <- do.call(rbind, vif_history)

env_selected <- names(vif_work)

env_scaled_selected <- env_scaled_all[, env_selected, drop = FALSE]

final_vif <- data.frame(variable = names(calc_vif(env_scaled_selected)), VIF = as.numeric(calc_vif(env_scaled_selected)),
  stringsAsFactors = FALSE)

write.csv(vif_history, file.path(table_dir, "环境变量_VIF筛查过程.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(final_vif, file.path(table_dir, "环境变量_最终VIF.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cor_mat <- cor(env_numeric, use = "pairwise.complete.obs")

cor_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)

cor_long <- data.frame(variable_1 = rownames(cor_mat)[cor_idx[, 1]], variable_2 = colnames(cor_mat)[cor_idx[,
  2]], correlation = cor_mat[cor_idx], abs_correlation = abs(cor_mat[cor_idx]), stringsAsFactors = FALSE)

cor_long <- cor_long[order(-cor_long$abs_correlation), ]

write.csv(cor_long, file.path(table_dir, "环境变量_两两相关.csv"), row.names = FALSE, fileEncoding = "UTF-8")

make_design <- function(meta_sub, env_sub) {
  time_mm <- model.matrix(~month, data = meta_sub)[, -1, drop = FALSE]
  space_mm <- model.matrix(~site, data = meta_sub)[, -1, drop = FALSE]
  colnames(time_mm) <- paste0("T_", seq_len(ncol(time_mm)))
  colnames(space_mm) <- paste0("S_", seq_len(ncol(space_mm)))
  env_mm <- as.matrix(env_sub)
  colnames(env_mm) <- paste0("E_", make.names(colnames(env_mm), unique = TRUE))
  list(Time = as.data.frame(time_mm), Environment = as.data.frame(env_mm), Space = as.data.frame(space_mm))
}

make_formula <- function(response_name, constraints, conditions = character()) {
  rhs <- paste(constraints, collapse = " + ")
  if (length(conditions) > 0) {
    rhs <- paste0(rhs, " + Condition(", paste(conditions, collapse = " + "), ")")
  }
  as.formula(paste(response_name, "~", rhs), env = parent.frame())
}

fit_ordination <- function(method, response, model_data, constraints, conditions = character()) {
  response_obj <- response
  form <- make_formula("response_obj", constraints, conditions)
  if (method == "RDA") {
    rda(form, data = model_data)
  }
  else {
    dbrda(form, data = model_data, sqrt.dist = TRUE)
  }
}

extract_r2 <- function(model) {
  rr <- RsquareAdj(model)
  c(raw_R2 = unname(rr$r.squared), adjusted_R2 = unname(rr$adj.r.squared))
}

fraction_labels <- c(`[a]` = "Pure Time", `[b]` = "Pure Environment", `[c]` = "Pure Space", `[d]` = "Time ∩ Environment",
  `[e]` = "Time ∩ Space", `[f]` = "Environment ∩ Space", `[g]` = "Time ∩ Environment ∩ Space",
  `[h]` = "Residual")

extract_fractions <- function(vp, group_name, method_name) {
  z <- as.data.frame(vp$part$indfract)
  if (nrow(z) != 8L)
    stop("varpart 独立分量表不是预期的 8 行。")
  code <- sub(" .*", "", rownames(z))
  adj_col <- grep("Adj", names(z), value = TRUE)[1]
  data.frame(method = method_name, group = group_name, fraction_code = code, fraction = unname(fraction_labels[code]),
    adjusted_R2 = z[[adj_col]], percent = 100 * z[[adj_col]], directly_testable = code %in% c("[a]",
      "[b]", "[c]"), stringsAsFactors = FALSE)
}

make_permutation_control <- function(meta_sub, group_name, target, nperm = 999) {
  if (target %in% c("Time", "Environment")) {
    blocks <- if (group_name == "All")
      interaction(meta_sub$site, meta_sub$size, drop = TRUE)
    else meta_sub$site
    scheme <- if (group_name == "All")
      "within Site × Size"
    else "within Site"
  }
  else {
    blocks <- if (group_name == "All")
      interaction(meta_sub$month, meta_sub$size, drop = TRUE)
    else meta_sub$month
    scheme <- if (group_name == "All")
      "within Month × Size"
    else "within Month"
  }
  list(control = how(nperm = nperm, blocks = blocks), scheme = scheme)
}

test_pure_fraction <- function(method, response, X, meta_sub, group_name, target, fraction_adj_r2, seed,
  nperm = 999, parallel_n = 4) {
  model_data <- as.data.frame(do.call(cbind, unname(X)), check.names = FALSE)
  target_cols <- colnames(X[[target]])
  condition_cols <- unlist(lapply(X[setdiff(names(X), target)], colnames), use.names = FALSE)
  model <- fit_ordination(method, response, model_data, target_cols, condition_cols)
  perm <- make_permutation_control(meta_sub, group_name, target, nperm)
  set.seed(seed)
  tab <- anova.cca(model, permutations = perm$control, parallel = parallel_n)
  r2 <- extract_r2(model)
  data.frame(method = method, group = group_name, fraction = paste("Pure", target), VPA_adjusted_R2 = fraction_adj_r2,
    model_adjusted_R2 = r2[["adjusted_R2"]], df = tab$Df[1], F = tab$F[1], p_value = tab$`Pr(>F)`[1],
    permutations_requested = nperm, permutation_scheme = perm$scheme, significant_0.05 = tab$`Pr(>F)`[1] <
      0.05, stringsAsFactors = FALSE)
}

fit_totals <- function(method, response, X, group_name) {
  model_data <- as.data.frame(do.call(cbind, unname(X)), check.names = FALSE)
  set_names <- c("Time", "Environment", "Space")
  rows <- lapply(set_names, function(nm) {
    mod <- fit_ordination(method, response, model_data, colnames(X[[nm]]))
    rr <- extract_r2(mod)
    data.frame(method = method, group = group_name, model = paste0(nm, " total"), raw_R2 = rr[["raw_R2"]],
      adjusted_R2 = rr[["adjusted_R2"]], stringsAsFactors = FALSE)
  })
  all_cols <- unlist(lapply(X, colnames), use.names = FALSE)
  mod_full <- fit_ordination(method, response, model_data, all_cols)
  rr_full <- extract_r2(mod_full)
  rows[[length(rows) + 1L]] <- data.frame(method = method, group = group_name, model = "Full Time + Environment + Space",
    raw_R2 = rr_full[["raw_R2"]], adjusted_R2 = rr_full[["adjusted_R2"]], stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

message("[3/7] 运行 All、Large、Small 的 RDA 与 dbRDA VPA")

group_indices <- list(All = rep(TRUE, nrow(meta)), Large = meta$size == "large", Small = meta$size ==
  "small")

vp_objects <- list(RDA = list(), dbRDA = list())

fraction_results <- list()

test_results <- list()

total_results <- list()

seed_counter <- 20260814L

for (group_name in names(group_indices)) {
  message("  - 分组：", group_name)
  idx <- group_indices[[group_name]]
  meta_sub <- droplevels(meta[idx, , drop = FALSE])
  counts_sub <- comm_all[idx, , drop = FALSE]
  env_sub <- env_scaled_selected[idx, , drop = FALSE]
  X <- make_design(meta_sub, env_sub)
  comm_hel <- decostand(counts_sub, method = "hellinger")
  bray <- vegdist(counts_sub, method = "bray")
  vp_rda <- varpart(comm_hel, X$Time, X$Environment, X$Space)
  vp_dbrda <- varpart(bray, X$Time, X$Environment, X$Space, sqrt.dist = TRUE)
  vp_objects$RDA[[group_name]] <- vp_rda
  vp_objects$dbRDA[[group_name]] <- vp_dbrda
  frac_rda <- extract_fractions(vp_rda, group_name, "RDA-Hellinger")
  frac_dbrda <- extract_fractions(vp_dbrda, group_name, "dbRDA-Bray")
  fraction_results[[paste(group_name, "RDA")]] <- frac_rda
  fraction_results[[paste(group_name, "dbRDA")]] <- frac_dbrda
  total_results[[paste(group_name, "RDA")]] <- fit_totals("RDA", comm_hel, X, group_name)
  total_results[[paste(group_name, "dbRDA")]] <- fit_totals("dbRDA", bray, X, group_name)
  pure_map <- c(Time = "Pure Time", Environment = "Pure Environment", Space = "Pure Space")
  for (target in names(pure_map)) {
    seed_counter <- seed_counter + 1L
    adj <- frac_rda$adjusted_R2[match(pure_map[[target]], frac_rda$fraction)]
    test_results[[paste(group_name, "RDA", target)]] <- test_pure_fraction("RDA", comm_hel, X, meta_sub,
      group_name, target, adj, seed_counter)
    seed_counter <- seed_counter + 1L
    adj <- frac_dbrda$adjusted_R2[match(pure_map[[target]], frac_dbrda$fraction)]
    test_results[[paste(group_name, "dbRDA", target)]] <- test_pure_fraction("dbRDA", bray, X, meta_sub,
      group_name, target, adj, seed_counter)
  }
}

fractions <- do.call(rbind, fraction_results)

tests <- do.call(rbind, test_results)

totals <- do.call(rbind, total_results)

rownames(fractions) <- rownames(tests) <- rownames(totals) <- NULL

write.csv(fractions, file.path(table_dir, "VPA_全部调整R2分量.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(subset(fractions, method == "RDA-Hellinger"), file.path(table_dir, "VPA_RDA_Hellinger_分量.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

write.csv(subset(fractions, method == "dbRDA-Bray"), file.path(table_dir, "VPA_dbRDA_Bray_分量.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

write.csv(tests, file.path(table_dir, "VPA_纯效应受限置换检验.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(totals, file.path(table_dir, "VPA_各变量组总解释量与全模型.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

message("[4/7] 输出输入审计和设计信息")

input_audit <- data.frame(item = c("raw_taxa", "retained_nonzero_taxa", "samples", "large_samples", "small_samples",
  "sites", "sampling_months", "minimum_library_size", "maximum_library_size", "candidate_environment_variables",
  "selected_environment_variables", "VIF_threshold"), value = c(nrow(asv_raw), nrow(asv_num), nrow(meta),
  sum(meta$size == "large"), sum(meta$size == "small"), nlevels(meta$site), nlevels(meta$month), min(lib_sizes),
  max(lib_sizes), length(env_candidates), length(env_selected), vif_threshold), stringsAsFactors = FALSE)

write.csv(input_audit, file.path(table_dir, "输入与设计审计.csv"), row.names = FALSE, fileEncoding = "UTF-8")

env_decisions <- data.frame(category = c(rep("selected", length(env_selected)), rep("removed_by_VIF",
  length(removed_vif)), rep("removed_constant", length(constant_removed)), rep("excluded_from_environment",
  3)), variable = c(env_selected, removed_vif, constant_removed, "Longitude", "Latitude", "Q"), reason = c(rep(paste0("最终 VIF ≤ ",
  vif_threshold), length(env_selected)), rep(paste0("迭代 VIF > ", vif_threshold), length(removed_vif)),
  rep("无方差", length(constant_removed)), "由站点因子表示空间，避免重复计入", "由站点因子表示空间，避免重复计入",
  "按研究设定预先排除，并与距离衰减分析保持一致"), stringsAsFactors = FALSE)

write.csv(env_decisions, file.path(table_dir, "环境变量纳入与剔除清单.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

design_info <- data.frame(set = c("Time", "Environment", "Space"), representation = c(paste0("15-level sampling-month factor (",
  nlevels(meta$month) - 1, " dummy variables)"), paste(env_selected, collapse = ", "), paste0("4-level site factor (",
  nlevels(meta$site) - 1, " dummy variables)")), rationale = c("Non-parametric sampling-occasion effect; accommodates irregular month sequence",
  "All measured physicochemical variables retained after predictor-only VIF screening", "Only four repeatedly sampled stations; site factor is more defensible than dbMEM on duplicated coordinates"),
  stringsAsFactors = FALSE)

write.csv(design_info, file.path(table_dir, "VPA变量组定义.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("[5/7] 绘制 VPA 图和纯效应对比图")

plot_vpa_panels <- function(vp_list, method_label, file_stub) {
  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(1, 3), mar = c(2, 2, 4, 2), oma = c(0, 0, 2, 0))
    for (g in c("All", "Large", "Small")) {
      vp <- vp_list[[g]]
      frac <- fractions[fractions$method == method_label & fractions$group == g, ]
      explained <- 100 * (1 - frac$adjusted_R2[frac$fraction_code == "[h]"])
      plot(vp, Xnames = c("Time", "Env.", "Space"), bg = c("#E69F0088", "#009E7388", "#0072B288"),
        digits = 3, cutoff = -Inf, id.size = 0.9, cex = 1.05)
      title(main = sprintf("%s\nTotal adjusted R² = %.1f%%", g, explained), cex.main = 1.05)
    }
    mtext(paste0(method_label, " variation partitioning"), outer = TRUE, cex = 1.25, font = 2)
  }
  png(file.path(figure_dir, paste0(file_stub, ".png")), width = 3300, height = 1200, res = 260, bg = "white")
  draw()
  dev.off()
  pdf(file.path(figure_dir, paste0(file_stub, ".pdf")), width = 13, height = 5, useDingbats = FALSE)
  draw()
  dev.off()
}

plot_vpa_panels(vp_objects$RDA, "RDA-Hellinger", "VPA_RDA_Hellinger_All_Large_Small")

plot_vpa_panels(vp_objects$dbRDA, "dbRDA-Bray", "VPA_dbRDA_Bray_All_Large_Small")

plot_pure_bars <- function(method_label, file_stub) {
  z <- tests[tests$method == ifelse(method_label == "RDA-Hellinger", "RDA", "dbRDA"), ]
  z$group <- factor(z$group, levels = c("All", "Large", "Small"))
  z$fraction <- factor(z$fraction, levels = c("Pure Time", "Pure Environment", "Pure Space"))
  mat <- xtabs(100 * VPA_adjusted_R2 ~ fraction + group, data = z)
  pmat <- xtabs(p_value ~ fraction + group, data = z)
  cols <- c("#E69F00", "#009E73", "#0072B2")
  ymin <- min(0, mat, na.rm = TRUE)
  ymax <- max(0, mat, na.rm = TRUE)
  pad <- max(1, 0.18 * (ymax - ymin + 1))
  png(file.path(figure_dir, paste0(file_stub, ".png")), width = 2200, height = 1500, res = 250, bg = "white")
  old <- par(mar = c(5, 5, 4, 1) + 0.1)
  bp <- barplot(mat, beside = TRUE, col = cols, border = NA, ylim = c(ymin - pad, ymax + 2 * pad),
    ylab = "Pure adjusted R² (%)", main = paste0(method_label, ": independent fractions"), las = 1)
  abline(h = 0, col = "grey35", lwd = 1)
  text(bp, mat, labels = sprintf("%.1f%%\n%s", mat, ifelse(pmat < 0.05, "*", "ns")), pos = ifelse(mat >=
    0, 3, 1), cex = 0.78)
  legend("topright", legend = rownames(mat), fill = cols, border = NA, bty = "n")
  mtext("* restricted-permutation p < 0.05; ns = not significant", side = 1, line = 4, cex = 0.8)
  par(old)
  dev.off()
}

plot_pure_bars("RDA-Hellinger", "VPA_RDA_Hellinger_纯效应对比")

plot_pure_bars("dbRDA-Bray", "VPA_dbRDA_Bray_纯效应对比")

message("[6/7] 生成可直接使用的结果摘要")

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

fmt_p <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

summary_lines <- c("# 时间—环境—空间 VPA 结果摘要", "", "## 分析设定", "", paste0("- 样本：",
  nrow(meta), " 个（Large = ", sum(meta$size == "large"), "，Small = ", sum(meta$size == "small"),
  "），4 个站点 × 15 个采样月份 × 2 个孔径。"), paste0("- 群落矩阵：", nrow(asv_num),
  " 个非零 ASV；每样本 reads 为 ", min(lib_sizes), "–", max(lib_sizes), "。"), paste0("- 环境变量：预先排除 Q 后，",
  length(env_candidates), " 个候选变量经仅基于预测变量的 VIF ≤ ", vif_threshold, " 筛查后保留 ",
  length(env_selected), " 个：", paste(env_selected, collapse = "、"), "。"), paste0("- VIF 剔除：",
  ifelse(length(removed_vif), paste(removed_vif, collapse = "、"), "无"), "。"), "- 时间用 15 水平采样月份因子表示；空间用 4 水平站点因子表示；经纬度不再作为环境变量重复纳入。",
  "- 主分析为 Hellinger–RDA；Bray–Curtis dbRDA（sqrt.dist = TRUE）用于稳健性核对。",
  "- 纯时间和纯环境效应在站点内置换（All 另按孔径分层）；纯空间效应在月份内置换（All 另按孔径分层）；均请求 999 次置换。",
  "", "## 数值结果")

for (method_label in c("RDA-Hellinger", "dbRDA-Bray")) {
  summary_lines <- c(summary_lines, "", paste0("### ", method_label), "")
  test_method <- if (method_label == "RDA-Hellinger")
    "RDA"
  else "dbRDA"
  for (g in c("All", "Large", "Small")) {
    fz <- fractions[fractions$method == method_label & fractions$group == g, ]
    tz <- tests[tests$method == test_method & tests$group == g, ]
    explained <- 1 - fz$adjusted_R2[fz$fraction_code == "[h]"]
    get_test <- function(nm) tz[tz$fraction == nm, ][1, ]
    tt <- get_test("Pure Time")
    ee <- get_test("Pure Environment")
    ss <- get_test("Pure Space")
    shared <- fz[fz$fraction_code %in% c("[d]", "[e]", "[f]", "[g]"), ]
    largest_shared <- shared[which.max(shared$adjusted_R2), ]
    summary_lines <- c(summary_lines, paste0("- **", g, "**：总调整解释率 ", fmt_pct(explained),
      "；纯时间 ", fmt_pct(tt$VPA_adjusted_R2), "（p = ", fmt_p(tt$p_value), "）", "；纯环境 ",
      fmt_pct(ee$VPA_adjusted_R2), "（p = ", fmt_p(ee$p_value), "）", "；纯空间 ", fmt_pct(ss$VPA_adjusted_R2),
      "（p = ", fmt_p(ss$p_value), "）", "；最大共享分量为 ", largest_shared$fraction,
      " ", fmt_pct(largest_shared$adjusted_R2), "。"))
  }
}

negative_rows <- fractions[fractions$adjusted_R2 < 0, ]

summary_lines <- c(summary_lines, "", "## 解释边界", "", "- VPA 报告 adjusted R²；共享分量不能直接进行显著性检验，也不能单独证明因果。",
  "- 时间因子会吸收所有采样月份共同变化；与环境的共享分量可表述为季节环境梯度与时间演替的解释重叠。",
  "- 空间只有 4 个站点，纯空间效应代表稳定站点差异，不能仅凭该结果断言扩散限制。",
  "- All 结果包含未建模的孔径差异，Large/Small 分层结果应作为生态解释的主要依据。",
  paste0("- 小于 0 的 adjusted R² 分量共 ", nrow(negative_rows), " 个，已原样保留；它们应解释为接近 0 或不稳定，而不是负的生态作用。"),
  "", "详细数值见同目录 CSV；圆图中的圆位置和距离没有排序或效应大小含义。")

writeLines(summary_lines, file.path(analysis_dir, "结果解读.md"), useBytes = TRUE)

analysis_notes <- c("分析目的：将整体群落变异分解为时间、环境、空间的独立与共享解释量。",
  "", "时间：15 个实际采样月份的分类因子（14 个虚拟变量）。", "环境：预先排除 Q、保留 WT；其余实测水环境指标先在全体 120 样本中标准化，再按 VIF≤10 筛查；筛查不使用群落响应。",
  "空间：4 个固定重复采样站点的分类因子（3 个虚拟变量）。只有 4 个坐标，不进行不稳定的 dbMEM 展开。",
  "", "主结果：Hellinger 变换后的 ASV 矩阵进入 RDA 型 VPA。", "核对结果：原始抽平计数计算 Bray-Curtis 距离，进入 dbRDA 型 VPA，并使用 sqrt.dist=TRUE。",
  "解释量统一报告 adjusted R²。纯分量用 partial RDA/dbRDA 和受限置换检验；共享分量不做 p 值。",
  "负分量不删除、不归零，在解释时视为接近 0 或估计不稳定。", "", "置换结构：纯时间/纯环境在 Site 内置换（All 再按 Size 分层）；纯空间在 Month 内置换（All 再按 Size 分层）。",
  "All 分析未把孔径作为第四个解释组，因此用于全局描述；Large 与 Small 分层结果优先用于机制解释。")

writeLines(analysis_notes, file.path(analysis_dir, "分析说明.txt"), useBytes = TRUE)

message("[7/7] 保存会话信息并完成")

capture.output(sessionInfo(), file = file.path(table_dir, "R_sessionInfo.txt"))

message("完成。结果目录：", table_dir)

message("图片目录：", figure_dir)
