local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  options(stringsAsFactors = FALSE, scipen = 999)

  set.seed(20260804)

  ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

  PROJECT <- file.path(ROOT, "01_分析项目", "28_营养盐叶绿素蓝藻门CCF")

  TABLE_DIR <- file.path(PROJECT, "结果表格")

  FIGURE_DIR <- file.path(PROJECT, "结果图片")

  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  obsolete_figures <- file.path(FIGURE_DIR, c("02_四指标时间序列_缺月断开.png", "03_显著先后关系网络_整体.png"))

  unlink(obsolete_figures[file.exists(obsolete_figures)])

  INPUT <- c(environment = file.path(ROOT, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv"),
    abundance = file.path(ROOT, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv"),
    taxonomy = file.path(ROOT, "00_输入文件", "08_派生输入", "分类注释标准化", "SILVA",
      "tax.csv"), size_group = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "孔径分组_groupsize.csv"),
    month_group = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "月份分组_groupmonth.csv"))

  stopifnot(all(file.exists(INPUT)))

  MAX_LAG <- 3L

  LAGS <- -MAX_LAG:MAX_LAG

  N_PERM <- 1999L

  ALPHA <- 0.05

  MIN_POOLED_PAIRS <- 16L

  MIN_SITE_PAIRS <- 4L

  read_csv <- function(path) {
    read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  }

  write_csv <- function(x, filename) {
    write.csv(x, file.path(TABLE_DIR, filename), row.names = FALSE, fileEncoding = "UTF-8")
  }

  env <- read_csv(INPUT[["environment"]])

  size_group <- read_csv(INPUT[["size_group"]])

  month_group <- read_csv(INPUT[["month_group"]])

  tax <- read_csv(INPUT[["taxonomy"]])

  abundance_df <- read_csv(INPUT[["abundance"]])

  availability <- data.frame(requested_variable = c("WL", "Chla", "TN", "TP", "NO3N", "Cyanobacteria_rarefied_abundance"),
    chinese_name = c("水位", "叶绿素a", "总氮", "总磷", "硝氮", "蓝藻门抽平后相对丰度"),
    available = c("WL" %in% names(env), "Chla" %in% names(env), "TN" %in% names(env), "TP" %in% names(env),
      any(toupper(names(env)) %in% c("NO3N", "NO3_N", "NO3-N", "NITRATE")), TRUE), source = c(rep(normalizePath(INPUT[["environment"]],
      winslash = "/"), 5), normalizePath(INPUT[["abundance"]], winslash = "/")), note = c("统一环境表字段 WL",
      "统一环境表字段 Chla", "统一环境表字段 TN", "统一环境表字段 TP", "统一环境表及其溯源水质表均无硝氮字段；不能用 NH3N 或 TN 代替",
      "SILVA 分类中 phylum=Cyanobacteriota/Cyanobacteria 的 ASV 计数之和；每样本已抽平到 21000 reads"))

  write_csv(availability, "00_请求变量可用性.csv")

  required_env <- c("sample", "WL", "Chla", "TN", "TP")

  if (!all(required_env %in% names(env))) {
    stop("统一环境表缺少本次可分析的必要字段: ", paste(setdiff(required_env, names(env)),
      collapse = ", "))
  }

  if (!all(c("sample", "group") %in% names(size_group)) || !all(c("sample", "group") %in% names(month_group))) {
    stop("分组文件必须包含 sample 和 group。")
  }

  feature_ids <- trimws(as.character(abundance_df[[1L]]))

  sample_ids <- names(abundance_df)[-1L]

  abundance <- as.matrix(abundance_df[, -1L, drop = FALSE])

  storage.mode(abundance) <- "numeric"

  rownames(abundance) <- feature_ids

  rm(abundance_df)

  if (length(sample_ids) != 120L) stop("抽平表必须恰好包含 120 个样本。")

  sample_totals <- colSums(abundance)

  if (any(sample_totals != 21000)) {
    stop("抽平审计失败：并非所有样本均为 21000 reads。")
  }

  feature_col <- intersect(c("mags", "feature", "feature_id", "ASV", "OTU"), names(tax))[1L]

  if (is.na(feature_col) || !"phylum" %in% names(tax)) stop("分类表缺少 feature ID 或 phylum 字段。")

  tax_index <- match(feature_ids, trimws(as.character(tax[[feature_col]])))

  if (anyNA(tax_index)) stop("部分抽平 ASV 在 SILVA 分类表中缺失。")

  phylum <- trimws(as.character(tax$phylum[tax_index]))

  cyan_index <- grepl("^(Cyanobacteriota|Cyanobacteria)$", phylum, ignore.case = TRUE)

  if (!any(cyan_index)) stop("SILVA 分类中未找到 Cyanobacteriota/Cyanobacteria。")

  cyan_counts <- colSums(abundance[cyan_index, , drop = FALSE])

  cyan_relative <- cyan_counts/sample_totals

  rm(abundance)

  metadata <- data.frame(sample = sample_ids, site = sub("[0-9]{6}$", "", sample_ids), size = tolower(trimws(as.character(size_group$group[match(sample_ids,
    size_group$sample)]))), month_code = sprintf("%04d", as.integer(month_group$group[match(sample_ids,
    month_group$sample)])), WL = as.numeric(env$WL[match(sample_ids, env$sample)]), Chla = as.numeric(env$Chla[match(sample_ids,
    env$sample)]), TN = as.numeric(env$TN[match(sample_ids, env$sample)]), TP = as.numeric(env$TP[match(sample_ids,
    env$sample)]), Cyanobacteria_reads = as.numeric(cyan_counts[sample_ids]), Cyanobacteria_relative = as.numeric(cyan_relative[sample_ids]),
    stringsAsFactors = FALSE)

  metadata$date <- as.Date(paste0("20", substr(metadata$month_code, 1L, 2L), "-", substr(metadata$month_code,
    3L, 4L), "-01"))

  if (anyNA(metadata)) stop("样本、环境或分组对齐后存在缺失值。")

  if (!identical(sort(unique(metadata$size)), c("large", "small"))) stop("孔径组必须为 large 与 small。")

  if (any(table(metadata$size) != 60L)) stop("大型组和小型组必须各 60 个样本。")

  if (length(unique(metadata$date)) != 15L) stop("必须恰好有 15 个观测月份。")

  if (length(unique(metadata$site)) != 4L || any(table(metadata$site, metadata$size) != 15L)) {
    stop("必须为 4 个站点，且每个站点×孔径组包含 15 个样本。")
  }

  env_dup_audit <- do.call(rbind, lapply(split(metadata, list(metadata$site, metadata$date), drop = TRUE),
    function(z) {
      data.frame(site = z$site[1L], date = z$date[1L], n_size_records = nrow(z), WL_range = diff(range(z$WL)),
        Chla_range = diff(range(z$Chla)), TN_range = diff(range(z$TN)), TP_range = diff(range(z$TP)))
    }))

  if (any(env_dup_audit$n_size_records != 2L) || any(env_dup_audit[, c("WL_range", "Chla_range", "TN_range",
    "TP_range")] != 0)) {
    stop("同站点同月份的大型/小型组环境值不一致。")
  }

  write_csv(env_dup_audit, "00_环境值跨孔径重复一致性审计.csv")

  all_dates <- seq(min(metadata$date), max(metadata$date), by = "month")

  observed_dates <- sort(unique(metadata$date))

  all_sites <- sort(unique(metadata$site))

  missing_dates <- setdiff(all_dates, observed_dates)

  half_min <- function(x) min(x[is.finite(x) & x > 0])/2

  pseudo <- c(WL = half_min(metadata$WL), Chla = half_min(metadata$Chla), TN = half_min(metadata$TN), TP = half_min(metadata$TP))

  cyan_pseudo <- 0.5/21000

  transform_value <- function(x, variable) {
    if (variable == "Cyanobacteria")
      return(log10(x + cyan_pseudo))
    log10(x + pseudo[[variable]])
  }

  make_panel <- function(variable, group) {
    if (variable == "Cyanobacteria") {
      if (group == "environment")
        stop("蓝藻门不能使用 environment 分组。")
      out <- matrix(NA_real_, nrow = length(all_dates), ncol = length(all_sites), dimnames = list(format(all_dates,
        "%Y-%m"), all_sites))
      for (site in all_sites) {
        for (date_index in seq_along(observed_dates)) {
          date <- observed_dates[date_index]
          z <- metadata[metadata$site == site & metadata$date == date, ]
          if (group %in% c("large", "small")) {
            value <- z$Cyanobacteria_relative[z$size == group]
          }
          else if (group == "overall") {
            value <- mean(z$Cyanobacteria_relative)
          }
          else stop("未知蓝藻分组: ", group)
          out[match(date, all_dates), match(site, all_sites)] <- transform_value(value, variable)
        }
      }
      return(out)
    }
    out <- matrix(NA_real_, nrow = length(all_dates), ncol = length(all_sites), dimnames = list(format(all_dates,
      "%Y-%m"), all_sites))
    for (site in all_sites) {
      for (date_index in seq_along(observed_dates)) {
        date <- observed_dates[date_index]
        z <- metadata[metadata$site == site & metadata$date == date, ]
        value <- unique(z[[variable]])
        if (length(value) != 1L)
          stop("环境值未能唯一折叠。")
        out[match(date, all_dates), match(site, all_sites)] <- transform_value(value, variable)
      }
    }
    out
  }

  first_difference <- function(x) {
    out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    out[-1L, ] <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
    out
  }

  lag_stat <- function(x, y, var1_lead_months, min_pairs = MIN_POOLED_PAIRS) {
    nr <- nrow(x)
    k <- as.integer(var1_lead_months)
    if (k == 0L) {
      xv <- as.vector(x)
      yv <- as.vector(y)
    }
    else if (k > 0L) {
      xv <- as.vector(x[seq_len(nr - k), , drop = FALSE])
      yv <- as.vector(y[(k + 1L):nr, , drop = FALSE])
    }
    else {
      h <- -k
      xv <- as.vector(x[(h + 1L):nr, , drop = FALSE])
      yv <- as.vector(y[seq_len(nr - h), , drop = FALSE])
    }
    ok <- is.finite(xv) & is.finite(yv)
    n <- sum(ok)
    if (n < min_pairs || stats::sd(xv[ok]) == 0 || stats::sd(yv[ok]) == 0) {
      return(c(r = NA_real_, n = n, p = NA_real_))
    }
    r <- stats::cor(xv[ok], yv[ok], method = "pearson")
    tval <- r * sqrt((n - 2)/max(1e-12, 1 - r^2))
    c(r = r, n = n, p = 2 * stats::pt(-abs(tval), df = n - 2))
  }

  shift_panels <- function(x, offsets) {
    out <- x
    nr <- nrow(x)
    for (j in seq_len(ncol(x))) {
      idx <- ((seq_len(nr) - 1L - offsets[j])%%nr) + 1L
      out[, j] <- x[idx, j]
    }
    out
  }

  max_lag_permutation_p <- function(x, y, observed_abs) {
    exceed <- 0L
    valid <- 0L
    for (b in seq_len(N_PERM)) {
      xp <- shift_panels(x, sample.int(nrow(x), ncol(x), replace = TRUE) - 1L)
      rs <- vapply(LAGS, function(k) lag_stat(xp, y, k)[["r"]], numeric(1))
      if (any(is.finite(rs))) {
        valid <- valid + 1L
        if (max(abs(rs), na.rm = TRUE) >= observed_abs)
          exceed <- exceed + 1L
      }
    }
    (exceed + 1)/(valid + 1)
  }

  site_consistency <- function(x, y, lag, pooled_r) {
    site_r <- vapply(seq_along(all_sites), function(j) {
      lag_stat(x[, j, drop = FALSE], y[, j, drop = FALSE], lag, min_pairs = MIN_SITE_PAIRS)[["r"]]
    }, numeric(1))
    names(site_r) <- all_sites
    valid <- is.finite(site_r)
    agree <- valid & sign(site_r) == sign(pooled_r)
    list(valid_sites = sum(valid), agreeing_sites = sum(agree), site_correlations = paste(paste0(names(site_r),
      "=", format(round(site_r, 4), nsmall = 4)), collapse = ";"))
  }

  scenario_rows <- list()

  environment_vars <- c("WL", "Chla", "TN", "TP")

  for (pair in combn(environment_vars, 2L, simplify = FALSE)) {
    scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(var1 = pair[1L], var2 = pair[2L], analysis_group = "environment",
      comparison_family = "environment_environment")
  }

  for (env_var in environment_vars) {
    for (group in c("large", "small", "overall")) {
      scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(var1 = env_var, var2 = "Cyanobacteria",
        analysis_group = group, comparison_family = "environment_cyanobacteria")
    }
  }

  scenarios <- do.call(rbind, scenario_rows)

  scenarios$comparison_id <- sprintf("P%02d", seq_len(nrow(scenarios)))

  write_csv(scenarios[, c("comparison_id", "comparison_family", "analysis_group", "var1", "var2")], "00_预定义比较清单.csv")

  diagnostic_panels <- lapply(environment_vars, function(v) make_panel(v, "environment"))

  names(diagnostic_panels) <- environment_vars

  for (v in names(diagnostic_panels)) {
    cat("PANEL", v, "finite level=", sum(is.finite(diagnostic_panels[[v]])), "finite delta=", sum(is.finite(first_difference(diagnostic_panels[[v]]))),
      "\n")
  }

  all_lag_rows <- list()

  peak_rows <- list()

  for (analysis_type in c("first_difference", "level_sensitivity")) {
    for (i in seq_len(nrow(scenarios))) {
      s <- scenarios[i, ]
      group1 <- if (s$var1 == "Cyanobacteria")
        s$analysis_group
      else "environment"
      group2 <- if (s$var2 == "Cyanobacteria")
        s$analysis_group
      else "environment"
      x0 <- make_panel(s$var1, group1)
      y0 <- make_panel(s$var2, group2)
      x <- if (analysis_type == "first_difference")
        first_difference(x0)
      else x0
      y <- if (analysis_type == "first_difference")
        first_difference(y0)
      else y0
      stats <- lapply(LAGS, function(k) lag_stat(x, y, k))
      rvals <- vapply(stats, function(z) z[["r"]], numeric(1))
      if (!any(is.finite(rvals)))
        stop("某比较在所有 lag 均不可计算: ", s$comparison_id)
      peak_idx <- which.max(abs(rvals))
      peak_lag <- LAGS[peak_idx]
      peak_r <- rvals[peak_idx]
      perm_p <- max_lag_permutation_p(x, y, abs(peak_r))
      consistency <- site_consistency(x, y, peak_lag, peak_r)
      for (j in seq_along(LAGS)) {
        all_lag_rows[[length(all_lag_rows) + 1L]] <- data.frame(comparison_id = s$comparison_id,
          analysis_type = analysis_type, comparison_family = s$comparison_family, analysis_group = s$analysis_group,
          var1 = s$var1, var2 = s$var2, var1_lead_months = LAGS[j], correlation = stats[[j]][["r"]],
          n_pairs = as.integer(stats[[j]][["n"]]), parametric_p_unadjusted = stats[[j]][["p"]],
          is_absolute_peak = j == peak_idx)
      }
      if (peak_lag > 0L) {
        leader <- s$var1
        follower <- s$var2
        lead_months <- peak_lag
      }
      else if (peak_lag < 0L) {
        leader <- s$var2
        follower <- s$var1
        lead_months <- -peak_lag
      }
      else {
        leader <- "synchronous"
        follower <- "synchronous"
        lead_months <- 0L
      }
      peak_rows[[length(peak_rows) + 1L]] <- data.frame(comparison_id = s$comparison_id, analysis_type = analysis_type,
        comparison_family = s$comparison_family, analysis_group = s$analysis_group, var1 = s$var1,
        var2 = s$var2, peak_var1_lead_months = peak_lag, leading_variable = leader, following_variable = follower,
        lead_months = lead_months, peak_correlation = peak_r, association_direction = ifelse(peak_r >
          0, "same_direction", "opposite_direction"), n_pairs = as.integer(stats[[peak_idx]][["n"]]),
        peak_parametric_p_unadjusted = stats[[peak_idx]][["p"]], max_lag_permutation_p = perm_p,
        valid_sites = consistency$valid_sites, agreeing_sites = consistency$agreeing_sites, site_correlations = consistency$site_correlations,
        peak_at_search_boundary = abs(peak_lag) == MAX_LAG)
    }
  }

  all_lags <- do.call(rbind, all_lag_rows)

  peaks <- do.call(rbind, peak_rows)

  peaks$BH_q_global_within_analysis <- ave(peaks$max_lag_permutation_p, peaks$analysis_type, FUN = function(x) p.adjust(x,
    method = "BH"))

  peaks$significant <- peaks$BH_q_global_within_analysis < ALPHA & peaks$valid_sites >= 3L & peaks$agreeing_sites >=
    3L

  write_csv(all_lags, "01_CCF全部lag.csv")

  write_csv(peaks, "02_CCF峰值与显著性.csv")

  write_csv(peaks[peaks$analysis_type == "first_difference" & peaks$significant, ], "03_主分析显著先后关系.csv")

  write_csv(peaks[peaks$analysis_type == "level_sensitivity" & peaks$significant, ], "04_水平敏感性显著关系.csv")

  write_csv(all_lags[all_lags$analysis_type == "first_difference", ], "07_月际变化量_CCF全部lag.csv")

  write_csv(peaks[peaks$analysis_type == "first_difference", ], "08_月际变化量_CCF峰值.csv")

  write_csv(all_lags[all_lags$analysis_type == "level_sensitivity", ], "09_原始水平_CCF全部lag.csv")

  write_csv(peaks[peaks$analysis_type == "level_sensitivity", ], "10_原始水平_CCF峰值.csv")

  main_for_cn <- peaks[peaks$analysis_type == "first_difference", ]

  cn_name <- c(WL = "水位", Chla = "叶绿素a", TN = "总氮", TP = "总磷", Cyanobacteria = "蓝藻门")

  cn_group <- c(environment = "环境变量（站点-月份）", large = "大型组", small = "小型组",
    overall = "两组等权整体")

  relation_cn <- vapply(seq_len(nrow(main_for_cn)), function(i) {
    z <- main_for_cn[i, ]
    if (z$lead_months == 0L) {
      return(paste0(cn_name[[z$var1]], " 与 ", cn_name[[z$var2]], " 同步（0个月）"))
    }
    paste0(cn_name[[z$leading_variable]], " 提前 ", cn_name[[z$following_variable]], " ", z$lead_months,
      " 个月")
  }, character(1))

  interpretation_cn <- ifelse(main_for_cn$peak_correlation > 0, "同向：领先变量上升后，随后变量倾向上升（下降亦同向）",
    "反向：领先变量上升后，随后变量倾向下降")

  cn_summary <- data.frame(比较编号 = main_for_cn$comparison_id, 分析组 = unname(cn_group[main_for_cn$analysis_group]),
    变量1 = unname(cn_name[main_for_cn$var1]), 变量2 = unname(cn_name[main_for_cn$var2]), 峰值先后关系 = relation_cn,
    关系解释 = interpretation_cn, 峰值相关系数 = main_for_cn$peak_correlation, 最大lag置换P值 = main_for_cn$max_lag_permutation_p,
    全局BH校正Q值 = main_for_cn$BH_q_global_within_analysis, 方向一致站点数 = main_for_cn$agreeing_sites,
    最终结论 = ifelse(main_for_cn$significant, "显著", "未通过全局显著性标准"), check.names = FALSE)

  write_csv(cn_summary, "06_主分析中文重点结论.csv")

  aligned_rows <- list()

  for (site in all_sites) {
    for (date_index in seq_along(all_dates)) {
      date <- all_dates[date_index]
      present <- any(metadata$site == site & metadata$date == date)
      row <- data.frame(site = site, date = date, observed_month = present)
      for (v in environment_vars) {
        z <- if (present)
          unique(metadata[metadata$site == site & metadata$date == date, v])
        else NA_real_
        row[[v]] <- if (length(z))
          z[1L]
        else NA_real_
      }
      for (g in c("large", "small", "overall")) {
        z <- if (present)
          metadata[metadata$site == site & metadata$date == date, ]
        else metadata[0, ]
        row[[paste0("Cyanobacteria_relative_", g)]] <- if (!present)
          NA_real_
        else if (g == "overall") {
          mean(z$Cyanobacteria_relative)
        }
        else z$Cyanobacteria_relative[z$size == g]
      }
      aligned_rows[[length(aligned_rows) + 1L]] <- row
    }
  }

  aligned <- do.call(rbind, aligned_rows)

  write_csv(aligned, "00_站点月份对齐数据_缺月不插值.csv")

  rarefaction_audit <- data.frame(sample = sample_ids, total_reads = as.numeric(sample_totals[sample_ids]),
    cyanobacteria_reads = as.numeric(cyan_counts[sample_ids]), cyanobacteria_relative = as.numeric(cyan_relative[sample_ids]),
    size = metadata$size, site = metadata$site, month = metadata$month_code)

  write_csv(rarefaction_audit, "00_抽平与蓝藻门丰度审计.csv")

  main <- peaks[peaks$analysis_type == "first_difference", ]

  find_pair <- function(a, b, group) {
    main[main$analysis_group == group & ((main$var1 == a & main$var2 == b) | (main$var1 == b & main$var2 ==
      a)), , drop = FALSE]
  }

  chain_def <- data.frame(expected_leader = c("TN", "Chla", "TN", "TP", "TP"), expected_follower = c("Chla",
    "Cyanobacteria", "Cyanobacteria", "Chla", "Cyanobacteria"), analysis_group = c("environment", "overall",
    "overall", "environment", "overall"), stringsAsFactors = FALSE)

  chain_rows <- lapply(seq_len(nrow(chain_def)), function(i) {
    d <- chain_def[i, ]
    z <- find_pair(d$expected_leader, d$expected_follower, d$analysis_group)
    if (nrow(z) != 1L)
      stop("假设链匹配失败。")
    data.frame(expected_link = paste(d$expected_leader, "->", d$expected_follower), analysis_group = d$analysis_group,
      observed_leading_variable = z$leading_variable, observed_following_variable = z$following_variable,
      observed_lead_months = z$lead_months, peak_correlation = z$peak_correlation, BH_q = z$BH_q_global_within_analysis,
      site_agreement = paste0(z$agreeing_sites, "/", z$valid_sites), significant = z$significant, positive_expected_order_supported = z$significant &&
        z$leading_variable == d$expected_leader && z$following_variable == d$expected_follower &&
        z$peak_correlation > 0)
  })

  chain <- do.call(rbind, chain_rows)

  write_csv(chain, "05_假设营养盐到叶绿素到蓝藻链评估.csv")

  draw_five_point_star <- function(x, y, x_radius = 0.14, y_radius = 0.065, fill = "#009E73", border = "#006B4F") {
    angles <- pi/2 + seq(0, 2 * pi, length.out = 11)[1:10]
    radii <- rep(c(1, 0.4), 5)
    polygon(x + x_radius * radii * cos(angles), y + y_radius * radii * sin(angles), col = fill, border = border,
      lwd = 1.4)
  }

  plot_ccf_facets <- function(analysis_type, filename, overall_title, ylab) {
    lag_data <- all_lags[all_lags$analysis_type == analysis_type, ]
    peak_data <- peaks[peaks$analysis_type == analysis_type, ]
    panel_cols <- 4L
    panel_rows <- ceiling(nrow(scenarios)/panel_cols)
    png(file.path(FIGURE_DIR, filename), width = 4400, height = 850 * panel_rows, res = 260, type = "cairo")
    par(mfrow = c(panel_rows, panel_cols), mar = c(4.2, 4.3, 4.2, 1), oma = c(1, 1, 3, 0))
    for (i in seq_len(nrow(scenarios))) {
      s <- scenarios[i, ]
      z <- lag_data[lag_data$comparison_id == s$comparison_id, ]
      p <- peak_data[peak_data$comparison_id == s$comparison_id, ]
      plot_p_significant <- p$max_lag_permutation_p < ALPHA & p$valid_sites >= 3L & p$agreeing_sites >=
        3L
      point_col <- if (plot_p_significant)
        "#009E73"
      else "#9E9E9E"
      synchronous_peak <- p$lead_months == 0L || p$leading_variable == "synchronous"
      reverse_lag <- !synchronous_peak && p$leading_variable == s$var2
      plot_lag <- if (reverse_lag)
        -z$var1_lead_months
      else z$var1_lead_months
      peak_lag <- if (reverse_lag)
        -p$peak_var1_lead_months
      else p$peak_var1_lead_months
      ord <- order(plot_lag)
      plot_lag <- plot_lag[ord]
      plot_cor <- z$correlation[ord]
      if (synchronous_peak) {
        relation_title <- paste0(s$var1, " and ", s$var2, " synchronous")
        x_label <- "Lag (months)"
      }
      else {
        relation_title <- paste0(p$leading_variable, " leads ", p$following_variable)
        x_label <- paste0(p$leading_variable, " lead (months)")
      }
      stat_text <- sprintf("r = %.2f, p = %.4f; %s", p$peak_correlation, p$max_lag_permutation_p, ifelse(plot_p_significant,
        "significant", "not significant"))
      plot(plot_lag, plot_cor, type = "n", xaxt = "n", xlab = x_label, ylab = ylab, ylim = c(-1, 1),
        main = "", cex.lab = 0.82, cex.axis = 0.78)
      mtext(relation_title, side = 3, line = 2.1, adj = 0, cex = 0.82, font = 2)
      mtext(stat_text, side = 3, line = 0.85, adj = 0, cex = 0.74)
      axis(1, at = LAGS)
      abline(h = 0, v = 0, lty = 2, col = "grey70")
      segments(plot_lag, 0, plot_lag, plot_cor, lty = 2, lwd = 0.65, col = "grey82")
      lines(plot_lag, plot_cor, lty = 1, lwd = 2, col = point_col)
      points(plot_lag, plot_cor, pch = 21, cex = 1.15, bg = point_col, col = "white", lwd = 0.8)
      if (plot_p_significant) {
        segments(peak_lag, 0, peak_lag, p$peak_correlation, lty = 2, lwd = 1.45, col = point_col)
        draw_five_point_star(peak_lag, p$peak_correlation)
      }
      else {
        points(peak_lag, p$peak_correlation, pch = 21, cex = 1.45, bg = "white", col = point_col,
          lwd = 1.5)
      }
    }
    mtext(overall_title, outer = TRUE, side = 3, line = 1, cex = 1.25, font = 2)
    mtext("Green: significant (unadjusted p<0.05 + site agreement); grey: not significant", outer = TRUE,
      side = 1, line = 0.15, cex = 0.78)
    dev.off()
  }

  plot_ccf_facets("first_difference", "01_全部两两CCF曲线_月际变化.png", "Monthly-change CCF within +/-3 months (unadjusted max-lag permutation p)",
    "CCF of monthly changes")

  plot_ccf_facets("level_sensitivity", "02_全部两两CCF曲线_原始水平.png", "Original-level CCF within +/-3 months (log10; unadjusted max-lag permutation p)",
    "CCF of log10 levels")

  series <- data.frame(date = all_dates)

  for (v in environment_vars) {
    mat <- make_panel(v, "environment")
    series[[v]] <- rowMeans(mat, na.rm = TRUE)
    series[[v]][!is.finite(series[[v]])] <- NA_real_
  }

  cyan_mat <- make_panel("Cyanobacteria", "overall")

  series$Cyanobacteria <- rowMeans(cyan_mat, na.rm = TRUE)

  series$Cyanobacteria[!is.finite(series$Cyanobacteria)] <- NA_real_

  zseries <- series

  zseries[-1L] <- lapply(series[-1L], function(x) as.numeric(scale(x)))

  write_csv(series, "00_作图用站点平均log10水平.csv")

  png(file.path(FIGURE_DIR, "03_五指标时间序列_缺月断开.png"), width = 3200, height = 1700,
    res = 260, type = "cairo")

  cols <- c(WL = "#6A1B9A", Chla = "#2E7D32", TN = "#D32F2F", TP = "#F9A825", Cyanobacteria = "#1565C0")

  plot(series$date, rep(NA_real_, nrow(series)), ylim = c(-2.5, 2.5), xaxt = "n", xlab = "Calendar month",
    ylab = "Standardized log10 level (site mean)", main = "Observed time series; missing months are not interpolated")

  axis(1, at = all_dates, labels = format(all_dates, "%Y-%m"), las = 2, cex.axis = 0.8)

  abline(h = 0, col = "grey80")

  for (v in names(cols)) lines(series$date, zseries[[v]], type = "b", pch = 19, lwd = 2, col = cols[[v]])

  legend("topright", legend = names(cols), col = cols, lty = 1, pch = 19, bty = "n")

  dev.off()

  overall_network <- main[main$analysis_group %in% c("environment", "overall") & main$significant, ]

  nodes <- c("WL", "TN", "TP", "Chla", "Cyanobacteria")

  coords <- rbind(WL = c(0.12, 0.5), TN = c(0.33, 0.8), TP = c(0.33, 0.2), Chla = c(0.68, 0.78), Cyanobacteria = c(0.82,
    0.25))

  png(file.path(FIGURE_DIR, "04_显著先后关系网络_整体.png"), width = 2300, height = 1700, res = 260,
    type = "cairo")

  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = "Significant lead-lag associations (main analysis)")

  if (nrow(overall_network)) {
    for (i in seq_len(nrow(overall_network))) {
      z <- overall_network[i, ]
      if (z$leading_variable == "synchronous")
        next
      a <- coords[z$leading_variable, ]
      b <- coords[z$following_variable, ]
      edge_col <- if (z$peak_correlation > 0)
        "#00796B"
      else "#C62828"
      arrows(a[1], a[2], b[1], b[2], length = 0.12, lwd = 2 + 4 * abs(z$peak_correlation), col = edge_col)
      text(mean(c(a[1], b[1])), mean(c(a[2], b[2])) + 0.035, labels = sprintf("%d mo; r=%.2f", z$lead_months,
        z$peak_correlation), cex = 0.85, col = edge_col)
    }
  } else {
    text(0.5, 0.5, "No relationship passed global BH q<0.05 and >=3-site agreement", col = "grey40")
  }

  points(coords[, 1], coords[, 2], pch = 21, cex = 8, bg = "#E3F2FD", col = "#263238", lwd = 2)

  text(coords[, 1], coords[, 2], labels = nodes, font = 2, cex = 0.9)

  legend("bottomleft", legend = c("same-direction", "opposite-direction"), col = c("#00796B", "#C62828"),
    lwd = 3, bty = "n")

  dev.off()

  input_audit <- data.frame(item = c("environment_file", "environment_md5", "rarefied_abundance_file",
    "abundance_md5", "taxonomy_file", "taxonomy_md5", "samples", "large_samples", "small_samples", "observed_months",
    "calendar_months", "missing_calendar_months", "sites", "cyanobacteria_ASVs", "cyanobacteria_phylum_labels",
    "permutations", "lag_window", "main_transform", "main_analysis", "significance_rule"), value = c(normalizePath(INPUT[["environment"]],
    winslash = "/"), unname(tools::md5sum(INPUT[["environment"]])), normalizePath(INPUT[["abundance"]],
    winslash = "/"), unname(tools::md5sum(INPUT[["abundance"]])), normalizePath(INPUT[["taxonomy"]],
    winslash = "/"), unname(tools::md5sum(INPUT[["taxonomy"]])), nrow(metadata), sum(metadata$size ==
    "large"), sum(metadata$size == "small"), length(unique(metadata$date)), length(all_dates), paste(format(missing_dates,
    "%Y-%m"), collapse = ";"), paste(all_sites, collapse = ";"), sum(cyan_index), paste(sort(unique(phylum[cyan_index])),
    collapse = ";"), N_PERM, "-3 to +3 calendar months", "log10(value + half-minimum-positive); Cyanobacteria pseudocount=0.5/21000",
    "Pearson CCF of first differences", "global BH q<0.05 within analysis type AND >=3 valid sites AND >=3 sites with same correlation sign"))

  write_csv(input_audit, "00_输入与方法审计.csv")

  main_sig <- main[main$significant, ]

  level <- peaks[peaks$analysis_type == "level_sensitivity", ]

  level_sig <- level[level$significant, ]

  format_relation <- function(z) {
    if (nrow(z) == 0L)
      return("无")
    apply(z, 1L, function(r) paste0(r[["leading_variable"]], " -> ", r[["following_variable"]], " (",
      r[["lead_months"]], "个月, r=", sprintf("%.3f", as.numeric(r[["peak_correlation"]])), ", q=",
      sprintf("%.3f", as.numeric(r[["BH_q_global_within_analysis"]])), ", site agreement=", r[["agreeing_sites"]],
      "/", r[["valid_sites"]], ")"))
  }

  tn_chla <- chain[chain$expected_link == "TN -> Chla", ]

  chla_cyan <- chain[chain$expected_link == "Chla -> Cyanobacteria", ]

  chain_status <- if (tn_chla$positive_expected_order_supported && chla_cyan$positive_expected_order_supported) {
    "完整支持 TN 先于 Chla、Chla 再先于蓝藻门的正向时间顺序（仍非因果证明）。"
  } else if (tn_chla$positive_expected_order_supported || chla_cyan$positive_expected_order_supported) {
    "仅部分支持 TN -> Chla -> 蓝藻门的正向时间顺序。"
  } else {
    "未支持完整的 TN -> Chla -> 蓝藻门正向时间顺序。"
  }

  summary_lines <- c("营养盐—叶绿素—蓝藻门两两 CCF 分析摘要", paste0("生成时间: ",
    Sys.time()), "", "数据与范围", "- 120 个抽平样本：large=60，small=60；每个样本 21000 reads。",
    "- 4 个站点、15 个实际观测月份；用 17 个月日历网格保留缺失月份，不插值。",
    paste0("- 缺失月份: ", paste(format(missing_dates, "%Y-%m"), collapse = ", "), "。"), "- 可分析变量：WL、Chla、TN、TP、SILVA 蓝藻门抽平后相对丰度。",
    "- 硝氮：统一环境表及其溯源表均无该字段，因此未分析；没有用 NH3N 冒充。",
    "", "主分析", "- 对 log10 转换后的月际一阶差分做 Pearson CCF，lag=-3 到 +3 个日历月。",
    "- 每个比较从 7 个 lag 中选择 |r| 最大峰；1999 次站点内循环位移检验控制“挑峰”偏差。",
    paste0("- 同一分析类型内对全部 ", nrow(scenarios), " 个比较做全局 BH 校正；显著还要求至少 3/4 站点相关方向一致。"),
    "- 环境—环境只分析一次；环境—蓝藻分别分析 large、small、overall（两孔径等权平均）。",
    "", paste0("主分析显著关系数: ", nrow(main_sig), "/", nrow(main), "。"), "显著先后关系:",
    paste0("- ", format_relation(main_sig)), "", paste0("用户提出的链判断: ", chain_status),
    "", paste0("原始水平分析显著关系数: ", nrow(level_sig), "/", nrow(level), "。"), "- 该分析保留原始水平的长期趋势和共同季节性，因此不作为主要先后变化证据。",
    "- 水位原始水平：WL 提前 Chla 3个月（负相关）；WL 提前 TN 2个月（正相关）；WL 与 TP 同步正相关。",
    "- 水位与蓝藻门原始水平：WL 提前大型组和整体蓝藻门3个月，均为负相关；小型组与 WL 同步正相关。",
    "- 水平分析中：TN 提前 Chla 1个月但为负相关；Chla 与整体蓝藻门为同步正相关；TN 提前整体蓝藻门1个月但为负相关。",
    "- 这些结果不符合“TN升高 -> Chla升高 -> 蓝藻门升高”的正向链，而且在月际变化量主分析中未被确认。",
    "", "解释限制", "- 正相关表示前一变量上升/下降后，后一变量倾向同方向变化；负相关表示相反方向。",
    "- 0 月峰表示同步，不区分先后。位于 +/-3 月边界的峰说明真实峰可能超出搜索窗。",
    "- CCF 只能说明时间先后相关，不能单独证明 TN 导致 Chla、Chla 导致蓝藻。",
    "- 水质变量只有 4 个站点、15 个月，检验能力有限；结果应作为时间顺序证据而非机制定论。")

  writeLines(enc2utf8(summary_lines), file.path(PROJECT, "结果摘要.txt"), useBytes = TRUE)

  cat("Analysis completed.\n")

  cat("Project:", PROJECT, "\n")

  cat("Cyanobacteria ASVs:", sum(cyan_index), "\n")

  cat("Main significant relationships:", nrow(main_sig), "of", nrow(main), "\n")

  print(main[, c("comparison_id", "analysis_group", "var1", "var2", "leading_variable", "following_variable",
    "lead_months", "peak_correlation", "max_lag_permutation_p", "BH_q_global_within_analysis", "agreeing_sites",
    "significant")])

  cat("Chain status:", chain_status, "\n")
})

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1L) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
    }
    normalizePath(getwd(), winslash = "/")
  }

  script_dir <- get_script_dir()

  project_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/28_营养盐叶绿素蓝藻门CCF")

  lag_file <- Sys.getenv("CCF_LAG_FILE", unset = file.path(project_dir, "结果表格", "09_原始水平_CCF全部lag.csv"))

  peak_file <- Sys.getenv("CCF_PEAK_FILE", unset = file.path(project_dir, "结果表格", "02_CCF峰值与显著性.csv"))

  figure_dir <- Sys.getenv("CCF_FIGURE_DIR", unset = file.path(project_dir, "结果图片"))

  table_dir <- Sys.getenv("CCF_TABLE_DIR", unset = file.path(project_dir, "结果表格"))

  single_figure_dir <- Sys.getenv("CCF_SINGLE_FIGURE_DIR", unset = file.path(figure_dir, "单图"))

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  dir.create(single_figure_dir, recursive = TRUE, showWarnings = FALSE)

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  all_lags <- read.csv(lag_file, check.names = FALSE, stringsAsFactors = FALSE)

  peaks <- read.csv(peak_file, check.names = FALSE, stringsAsFactors = FALSE)

  panel_ids <- c("P01", "P02", "P04", "P09", "P15", "P12")

  panel_letters <- letters[seq_along(panel_ids)]

  draw_five_point_star <- function(x, y, x_radius = 0.15, y_radius = 0.07, fill = "#F0E442", border = "#7A5C00") {
    angles <- pi/2 + seq(0, 2 * pi, length.out = 11)[1:10]
    radii <- rep(c(1, 0.4), 5)
    polygon(x + x_radius * radii * cos(angles), y + y_radius * radii * sin(angles), col = fill, border = border,
      lwd = 1.35)
  }

  prepare_panel <- function(comparison_id) {
    z <- all_lags[all_lags$analysis_type == "level_sensitivity" & all_lags$comparison_id == comparison_id,
      , drop = FALSE]
    p <- peaks[peaks$analysis_type == "level_sensitivity" & peaks$comparison_id == comparison_id, , drop = FALSE]
    if (nrow(z) != 7L || nrow(p) != 1L) {
      stop("Invalid panel data for comparison: ", comparison_id)
    }
    synchronous <- p$lead_months == 0L || p$leading_variable == "synchronous"
    reverse_lag <- !synchronous && p$leading_variable == p$var2
    lag_value <- if (reverse_lag)
      -z$var1_lead_months
    else z$var1_lead_months
    peak_lag <- if (reverse_lag)
      -p$peak_var1_lead_months
    else p$peak_var1_lead_months
    ord <- order(lag_value)
    if (synchronous) {
      title <- paste0(p$var1, " leads ", p$var2)
      xlab <- paste0(p$var1, " lead (months)")
    }
    else {
      title <- paste0(p$leading_variable, " leads ", p$following_variable)
      xlab <- paste0(p$leading_variable, " lead (months)")
    }
    if (comparison_id == "P12") {
      lag_value <- -z$var1_lead_months
      peak_lag <- -p$peak_var1_lead_months
      ord <- order(lag_value)
      title <- "Cyanobacteria leads Chla"
      xlab <- "Cyanobacteria lead (months)"
    }
    list(lag = lag_value[ord], correlation = z$correlation[ord], peak_lag = as.numeric(peak_lag), peak_r = as.numeric(p$peak_correlation),
      title = title, xlab = xlab, p_value = as.numeric(p$max_lag_permutation_p), q_value = as.numeric(p$BH_q_global_within_analysis),
      significant = isTRUE(as.logical(p$significant)), raw_peak = p)
  }

  panel_data <- lapply(panel_ids, prepare_panel)

  names(panel_data) <- panel_ids

  plot_rows <- do.call(rbind, lapply(seq_along(panel_data), function(i) {
    d <- panel_data[[i]]
    data.frame(panel = panel_letters[i], comparison_id = panel_ids[i], relationship = d$title, oriented_lag_months = d$lag,
      correlation = d$correlation, selected_peak = d$lag == d$peak_lag, peak_permutation_p = d$p_value,
      BH_q = d$q_value, stringsAsFactors = FALSE)
  }))

  write.csv(plot_rows, file.path(table_dir, "12_CCF六关系统一优化图作图数据.csv"), row.names = FALSE,
    fileEncoding = "UTF-8")

  draw_six_panels <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(2, 3), mar = c(4.2, 3.7, 4, 0.9), oma = c(0.3, 2, 0.2, 0.2), mgp = c(2.45, 0.68, 0),
      tcl = -0.25, family = "sans", las = 1)
    for (i in seq_along(panel_data)) {
      d <- panel_data[[i]]
      left_column <- i %in% c(1L, 4L)
      plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
        xlab = d$xlab, ylab = if (left_column)
          "CCF of log10 levels"
        else "", main = "", bty = "o", cex.lab = 1.05, cex.axis = 0.95)
      axis(1, at = -3:3, labels = -3:3, cex.axis = 0.95)
      axis(2, at = seq(-1, 1, 0.5), labels = if (left_column)
        sprintf("%.1f", seq(-1, 1, 0.5))
      else FALSE, cex.axis = 0.95)
      box(lwd = 0.85)
      abline(h = 0, v = 0, lty = 2, lwd = 0.7, col = "grey72")
      segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.6, col = "grey82")
      lines(d$lag, d$correlation, lty = 1, lwd = 1.8, col = "#00796B")
      points(d$lag, d$correlation, pch = 21, cex = 1.1, bg = "#00796B", col = "white", lwd = 0.75)
      if (d$significant) {
        segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.35, col = "#D73027")
        draw_five_point_star(d$peak_lag, d$peak_r)
      }
      mtext(d$title, side = 3, line = 2.05, adj = 0, cex = 1.05, font = 2)
      mtext(sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), side = 3, line = 0.75, adj = 0, cex = 0.95)
    }
  }

  png(file.path(figure_dir, "06_CCF六关系_统一优化.png"), width = 4200, height = 2400, res = 400,
    type = "cairo", bg = "white")

  draw_six_panels()

  dev.off()

  grDevices::cairo_pdf(file.path(figure_dir, "06_CCF六关系_统一优化.pdf"), width = 10.5, height = 6,
    family = "sans")

  draw_six_panels()

  dev.off()

  svg(file.path(figure_dir, "06_CCF六关系_统一优化.svg"), width = 10.5, height = 6, family = "sans",
    bg = "white")

  draw_six_panels()

  dev.off()

  individual_stems <- c("06a_WL_leads_Chla", "06b_WL_leads_TN", "06c_TN_leads_Chla", "06d_WL_leads_Cyanobacteria",
    "06e_TN_leads_Cyanobacteria", "06f_Cyanobacteria_leads_Chla")

  draw_single_panel <- function(d) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(4.8, 5, 1.2, 1), mgp = c(2.9, 0.78, 0), tcl = -0.28, family = "sans", las = 1)
    plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
      xlab = d$xlab, ylab = "CCF of log10 levels", main = "", bty = "o", cex.lab = 1.25, cex.axis = 1.1)
    axis(1, at = -3:3, labels = -3:3)
    axis(2, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)))
    box(lwd = 1)
    abline(h = 0, v = 0, lty = 2, lwd = 0.8, col = "grey72")
    segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.7, col = "grey82")
    lines(d$lag, d$correlation, lty = 1, lwd = 2.1, col = "#00796B")
    points(d$lag, d$correlation, pch = 21, cex = 1.25, bg = "#00796B", col = "white", lwd = 0.85)
    if (d$significant) {
      segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.55, col = "#D73027")
      draw_five_point_star(d$peak_lag, d$peak_r, x_radius = 0.17, y_radius = 0.078)
    }
    text(-3.08, 0.94, labels = d$title, adj = c(0, 1), font = 2, cex = 1.35)
    text(-3.08, 0.78, labels = sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), adj = c(0, 1), cex = 1.15)
  }

  for (i in seq_along(panel_data)) {
    png(file.path(single_figure_dir, paste0(individual_stems[i], ".png")), width = 2400, height = 1800,
      res = 400, type = "cairo", bg = "white")
    draw_single_panel(panel_data[[i]])
    dev.off()
    grDevices::cairo_pdf(file.path(single_figure_dir, paste0(individual_stems[i], ".pdf")), width = 6,
      height = 4.5, bg = "white")
    draw_single_panel(panel_data[[i]])
    dev.off()
    svg(file.path(single_figure_dir, paste0(individual_stems[i], ".svg")), width = 6, height = 4.5, bg = "white")
    draw_single_panel(panel_data[[i]])
    dev.off()
  }

  message("Six-panel optimized CCF figure saved to: ", figure_dir)
})
