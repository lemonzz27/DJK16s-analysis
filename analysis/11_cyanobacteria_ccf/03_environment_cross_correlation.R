calculate_environment_ccf <- function(mode) {
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  options(stringsAsFactors = FALSE, scipen = 999)

  set.seed(20260817)

  ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

  PROJECT <- file.path(ROOT, "01_分析项目", "28_营养盐叶绿素蓝藻门CCF")

  SUPPLEMENT_MODE <- mode

  if (SUPPLEMENT_MODE == "AIR_TEMP") {
    OUTPUT_TAG <- "补充_气温AIRT_原始水平"
    REQUESTED_VARIABLES <- c("AIRT")
    REQUESTED_CN <- c(AIRT = "气温")
    REQUESTED_ALIASES <- list(AIRT = c("AIRT", "AIR_TEMP", "AIR-TEMP", "AIR_TEMPERATURE"))
  } else if (SUPPLEMENT_MODE == "N_P_RATIO") {
    OUTPUT_TAG <- "补充_氮磷比TNTP_原始水平"
    REQUESTED_VARIABLES <- c("TNTP")
    REQUESTED_CN <- c(TNTP = "氮磷比")
    REQUESTED_ALIASES <- list(TNTP = c("TNTP", "TN_TP", "TN:TP", "N_P_RATIO"))
  } else {
    OUTPUT_TAG <- "补充_TOC_WT_原始水平"
    REQUESTED_VARIABLES <- c("WT", "TOC")
    REQUESTED_CN <- c(WT = "水温", TOC = "总有机碳")
    REQUESTED_ALIASES <- list(WT = c("WT", "WATER_TEMP", "WATERTEMP", "TEMPERATURE", "TEMP"), TOC = c("TOC",
      "TOTAL_ORGANIC_CARBON", "TOTALORGANICCARBON"))
  }

  TABLE_DIR <- file.path(PROJECT, "结果表格", OUTPUT_TAG)

  FIGURE_DIR <- file.path(PROJECT, "结果图片", OUTPUT_TAG)

  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  INPUT <- c(environment = file.path(ROOT, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv"),
    air_temperature_derived = file.path(ROOT, "00_输入文件", "08_派生输入", "气温_AIRT_按现有水质匹配日期.csv"),
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

  find_column <- function(x, aliases) {
    idx <- match(toupper(aliases), toupper(names(x)), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx))
      names(x)[idx[1L]]
    else NA_character_
  }

  env <- read_csv(INPUT[["environment"]])

  if (SUPPLEMENT_MODE == "AIR_TEMP") {
    air_temperature <- read_csv(INPUT[["air_temperature_derived"]])
    if (!all(c("sample", "AIRT", "AIRT_day_gap_from_water_match") %in% names(air_temperature))) {
      stop("气温派生输入缺少sample、AIRT或匹配日期审计字段。")
    }
    air_index <- match(env$sample, air_temperature$sample)
    if (anyNA(air_index) || anyNA(air_temperature$AIRT[air_index])) {
      stop("气温派生输入未能与统一环境表样本一一完整匹配。")
    }
    if (any(air_temperature$AIRT_day_gap_from_water_match[air_index] != 0)) {
      stop("气温并非全部精确匹配既有水质匹配日期。")
    }
    env$AIRT <- as.numeric(air_temperature$AIRT[air_index])
  }

  if (SUPPLEMENT_MODE == "N_P_RATIO") {
    if (!all(c("sample", "TN", "TP") %in% names(env))) {
      stop("统一环境表缺少sample、TN或TP，不能计算氮磷比。")
    }
    env$TN <- as.numeric(env$TN)
    env$TP <- as.numeric(env$TP)
    if (anyNA(env[, c("TN", "TP")]) || any(env$TN <= 0) || any(env$TP <= 0)) {
      stop("TN或TP含缺失/非正值，不能直接计算并取log10氮磷比。")
    }
    env$TNTP <- env$TN/env$TP
    ratio_audit <- data.frame(sample = env$sample, TN = env$TN, TP = env$TP, TNTP_mass_ratio = env$TNTP,
      TNTP_molar_ratio = env$TNTP * 30.973761998/14.0067, definition = "mass ratio=TN/TP; molar ratio=(TN/14.0067)/(TP/30.973761998)")
    write_csv(ratio_audit, "00_氮磷比质量比与摩尔比_样本级.csv")
  }

  size_group <- read_csv(INPUT[["size_group"]])

  month_group <- read_csv(INPUT[["month_group"]])

  tax <- read_csv(INPUT[["taxonomy"]])

  abundance_df <- read_csv(INPUT[["abundance"]])

  matched_columns <- vapply(REQUESTED_VARIABLES, function(v) find_column(env, REQUESTED_ALIASES[[v]]),
    character(1))

  availability <- data.frame(requested_variable = REQUESTED_VARIABLES, chinese_name = unname(REQUESTED_CN[REQUESTED_VARIABLES]),
    available = !is.na(matched_columns), matched_column = ifelse(is.na(matched_columns), "", matched_columns),
    source = ifelse(REQUESTED_VARIABLES == "AIRT", normalizePath(INPUT[["air_temperature_derived"]],
      winslash = "/"), normalizePath(INPUT[["environment"]], winslash = "/")), note = ifelse(is.na(matched_columns),
      paste0("统一环境表没有", REQUESTED_VARIABLES, "字段；不使用其他指标替代"),
      ifelse(REQUESTED_VARIABLES == "AIRT", paste0("气温派生表字段已匹配: ", matched_columns,
        "；按既有水质matched_water_date精确连接，未插补"), ifelse(REQUESTED_VARIABLES ==
        "TNTP", paste0("由统一环境表TN/TP计算质量比: ", matched_columns, "；摩尔比仅相差常数30.973761998/14.0067，CCF结果相同"),
        paste0("统一环境表字段已匹配: ", matched_columns)))))

  write_csv(availability, "00_新增变量可用性审计.csv")

  if (!any(!is.na(matched_columns))) stop("统一环境表中未找到任何请求的新变量。")

  required_env <- c("sample", "WL", "Chla", "TN", "TP")

  if (!all(required_env %in% names(env))) {
    stop("统一环境表缺少必要字段: ", paste(setdiff(required_env, names(env)), collapse = ", "))
  }

  new_variables <- REQUESTED_VARIABLES[!is.na(matched_columns)]

  unavailable_variables <- REQUESTED_VARIABLES[is.na(matched_columns)]

  for (v in new_variables) env[[v]] <- as.numeric(env[[matched_columns[[v]]]])

  environment_variables <- c("WL", "Chla", "TN", "TP", new_variables)

  feature_ids <- trimws(as.character(abundance_df[[1L]]))

  sample_ids <- names(abundance_df)[-1L]

  abundance <- as.matrix(abundance_df[, -1L, drop = FALSE])

  storage.mode(abundance) <- "numeric"

  rownames(abundance) <- feature_ids

  rm(abundance_df)

  if (length(sample_ids) != 120L) stop("抽平表必须恰好包含120个样本。")

  sample_totals <- colSums(abundance)

  if (any(sample_totals != 21000)) stop("抽平审计失败：并非所有样本均为21000 reads。")

  feature_col <- intersect(c("mags", "feature", "feature_id", "ASV", "OTU"), names(tax))[1L]

  if (is.na(feature_col) || !"phylum" %in% names(tax)) stop("分类表缺少feature ID或phylum字段。")

  tax_index <- match(feature_ids, trimws(as.character(tax[[feature_col]])))

  if (anyNA(tax_index)) stop("部分抽平ASV在SILVA分类表中缺失。")

  phylum <- trimws(as.character(tax$phylum[tax_index]))

  cyan_index <- grepl("^(Cyanobacteriota|Cyanobacteria)$", phylum, ignore.case = TRUE)

  if (!any(cyan_index)) stop("SILVA分类中未找到Cyanobacteriota/Cyanobacteria。")

  cyan_counts <- colSums(abundance[cyan_index, , drop = FALSE])

  cyan_relative <- cyan_counts/sample_totals

  rm(abundance)

  metadata <- data.frame(sample = sample_ids, site = sub("[0-9]{6}$", "", sample_ids), size = tolower(trimws(as.character(size_group$group[match(sample_ids,
    size_group$sample)]))), month_code = sprintf("%04d", as.integer(month_group$group[match(sample_ids,
    month_group$sample)])), Cyanobacteria_reads = as.numeric(cyan_counts[sample_ids]), Cyanobacteria_relative = as.numeric(cyan_relative[sample_ids]),
    stringsAsFactors = FALSE)

  for (v in environment_variables) {
    metadata[[v]] <- as.numeric(env[[v]][match(sample_ids, env$sample)])
  }

  metadata$date <- as.Date(paste0("20", substr(metadata$month_code, 1L, 2L), "-", substr(metadata$month_code,
    3L, 4L), "-01"))

  required_metadata <- c("sample", "site", "size", "month_code", "date", "Cyanobacteria_relative", environment_variables)

  if (anyNA(metadata[, required_metadata, drop = FALSE])) {
    stop("样本、环境或分组对齐后存在缺失值。")
  }

  if (!identical(sort(unique(metadata$size)), c("large", "small")) || any(table(metadata$size) != 60L)) {
    stop("大型组和小型组必须各60个样本。")
  }

  if (length(unique(metadata$date)) != 15L) stop("必须恰好有15个观测月份。")

  if (length(unique(metadata$site)) != 4L || any(table(metadata$site, metadata$size) != 15L)) {
    stop("必须为4个站点，且每个站点×孔径组包含15个样本。")
  }

  env_dup_audit <- do.call(rbind, lapply(split(metadata, list(metadata$site, metadata$date), drop = TRUE),
    function(z) {
      row <- data.frame(site = z$site[1L], date = z$date[1L], n_size_records = nrow(z))
      for (v in environment_variables) row[[paste0(v, "_range")]] <- diff(range(z[[v]]))
      row
    }))

  range_cols <- paste0(environment_variables, "_range")

  if (any(env_dup_audit$n_size_records != 2L) || any(env_dup_audit[, range_cols] != 0)) {
    stop("同站点同月份的大型/小型组环境值不一致。")
  }

  write_csv(env_dup_audit, "00_环境值跨孔径重复一致性审计.csv")

  all_dates <- seq(min(metadata$date), max(metadata$date), by = "month")

  observed_dates <- sort(unique(metadata$date))

  all_sites <- sort(unique(metadata$site))

  missing_dates <- setdiff(all_dates, observed_dates)

  half_min <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (!length(positive))
      stop("变量没有正值，不能建立log10伪计数。")
    min(positive)/2
  }

  pseudo <- vapply(environment_variables, function(v) half_min(metadata[[v]]), numeric(1))

  cyan_pseudo <- 0.5/21000

  transform_value <- function(x, variable) {
    if (variable == "Cyanobacteria")
      return(log10(x + cyan_pseudo))
    if (variable == "AIRT")
      return(x)
    log10(x + pseudo[[variable]])
  }

  make_panel <- function(variable, group) {
    out <- matrix(NA_real_, nrow = length(all_dates), ncol = length(all_sites), dimnames = list(format(all_dates,
      "%Y-%m"), all_sites))
    for (site in all_sites) {
      for (date_index in seq_along(observed_dates)) {
        date <- observed_dates[date_index]
        z <- metadata[metadata$site == site & metadata$date == date, , drop = FALSE]
        if (variable == "Cyanobacteria") {
          if (group %in% c("large", "small")) {
            value <- z$Cyanobacteria_relative[z$size == group]
          }
          else if (group == "overall") {
            value <- mean(z$Cyanobacteria_relative)
          }
          else stop("未知蓝藻门分组: ", group)
        }
        else {
          value <- unique(z[[variable]])
          if (length(value) != 1L)
            stop("环境值未能唯一折叠: ", variable)
        }
        out[match(date, all_dates), match(site, all_sites)] <- transform_value(value, variable)
      }
    }
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

  base_variables <- c("WL", "Chla", "TN", "TP")

  for (new_var in new_variables) {
    for (base_var in base_variables) {
      scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(var1 = new_var, var2 = base_var, analysis_group = "environment",
        comparison_family = "new_environment_environment")
    }
    for (group in c("large", "small", "overall")) {
      scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(var1 = new_var, var2 = "Cyanobacteria",
        analysis_group = group, comparison_family = "new_environment_cyanobacteria")
    }
  }

  if (length(new_variables) >= 2L) {
    for (pair in combn(new_variables, 2L, simplify = FALSE)) {
      scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(var1 = pair[1L], var2 = pair[2L], analysis_group = "environment",
        comparison_family = "new_environment_environment")
    }
  }

  scenarios <- unique(do.call(rbind, scenario_rows))

  scenarios$structural_coupling <- SUPPLEMENT_MODE == "N_P_RATIO" & ((scenarios$var1 == "TNTP" & scenarios$var2 %in%
    c("TN", "TP")) | (scenarios$var2 == "TNTP" & scenarios$var1 %in% c("TN", "TP")))

  scenarios$interpretation_note <- ifelse(scenarios$structural_coupling, "TNTP is calculated from TN and TP; this comparison is mathematically coupled and is not independent ecological evidence",
    "")

  scenarios$comparison_id <- sprintf("SUP%02d", seq_len(nrow(scenarios)))

  scenarios <- scenarios[, c("comparison_id", "comparison_family", "analysis_group", "var1", "var2", "structural_coupling",
    "interpretation_note")]

  write_csv(scenarios, "00_新增比较清单.csv")

  all_lag_rows <- list()

  peak_rows <- list()

  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    group1 <- if (s$var1 == "Cyanobacteria")
      s$analysis_group
    else "environment"
    group2 <- if (s$var2 == "Cyanobacteria")
      s$analysis_group
    else "environment"
    x <- make_panel(s$var1, group1)
    y <- make_panel(s$var2, group2)
    stats <- lapply(LAGS, function(k) lag_stat(x, y, k))
    rvals <- vapply(stats, function(z) z[["r"]], numeric(1))
    if (!any(is.finite(rvals)))
      stop("某比较在所有lag均不可计算: ", s$comparison_id)
    peak_idx <- which.max(abs(rvals))
    peak_lag <- LAGS[peak_idx]
    peak_r <- rvals[peak_idx]
    perm_p <- max_lag_permutation_p(x, y, abs(peak_r))
    consistency <- site_consistency(x, y, peak_lag, peak_r)
    for (j in seq_along(LAGS)) {
      all_lag_rows[[length(all_lag_rows) + 1L]] <- data.frame(comparison_id = s$comparison_id, analysis_type = "original_level",
        comparison_family = s$comparison_family, analysis_group = s$analysis_group, var1 = s$var1,
        var2 = s$var2, var1_lead_months = LAGS[j], correlation = stats[[j]][["r"]], n_pairs = as.integer(stats[[j]][["n"]]),
        parametric_p_unadjusted = stats[[j]][["p"]], is_absolute_peak = j == peak_idx)
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
    significant <- perm_p < ALPHA && consistency$valid_sites >= 3L && consistency$agreeing_sites >= 3L
    peak_rows[[length(peak_rows) + 1L]] <- data.frame(comparison_id = s$comparison_id, analysis_type = "original_level",
      comparison_family = s$comparison_family, analysis_group = s$analysis_group, var1 = s$var1, var2 = s$var2,
      peak_var1_lead_months = peak_lag, leading_variable = leader, following_variable = follower, lead_months = lead_months,
      peak_correlation = peak_r, association_direction = ifelse(peak_r > 0, "same_direction", "opposite_direction"),
      structural_coupling = s$structural_coupling, interpretation_note = s$interpretation_note, n_pairs = as.integer(stats[[peak_idx]][["n"]]),
      peak_parametric_p_unadjusted = stats[[peak_idx]][["p"]], max_lag_permutation_p_unadjusted = perm_p,
      valid_sites = consistency$valid_sites, agreeing_sites = consistency$agreeing_sites, site_correlations = consistency$site_correlations,
      peak_at_search_boundary = abs(peak_lag) == MAX_LAG, significant_uncorrected = significant)
    cat(sprintf("%s: %s vs %s | %s | lag=%+d r=%.4f p=%.4f sites=%d/%d\n", s$comparison_id, s$var1, s$var2,
      s$analysis_group, peak_lag, peak_r, perm_p, consistency$agreeing_sites, consistency$valid_sites))
  }

  all_lags <- do.call(rbind, all_lag_rows)

  peaks <- do.call(rbind, peak_rows)

  significant_peaks <- peaks[peaks$significant_uncorrected, , drop = FALSE]

  cyan_peaks <- peaks[peaks$comparison_family == "new_environment_cyanobacteria", , drop = FALSE]

  write_csv(all_lags, "01_新增变量原始水平_CCF全部lag.csv")

  write_csv(peaks, "02_新增变量原始水平_CCF峰值_未校正p.csv")

  write_csv(significant_peaks, "03_显著新增关系_未校正p小于0.05.csv")

  write_csv(cyan_peaks, "04_新增变量与蓝藻门_CCF峰值.csv")

  aligned_rows <- list()

  for (site in all_sites) {
    for (date_index in seq_along(all_dates)) {
      date <- all_dates[date_index]
      z <- metadata[metadata$site == site & metadata$date == date, , drop = FALSE]
      row <- data.frame(site = site, date = date, observed = nrow(z) > 0L)
      for (v in environment_variables) {
        row[[v]] <- if (!nrow(z))
          NA_real_
        else unique(z[[v]])
      }
      for (g in c("large", "small", "overall")) {
        row[[paste0("Cyanobacteria_relative_", g)]] <- if (!nrow(z)) {
          NA_real_
        }
        else if (g == "overall") {
          mean(z$Cyanobacteria_relative)
        }
        else {
          z$Cyanobacteria_relative[z$size == g]
        }
      }
      aligned_rows[[length(aligned_rows) + 1L]] <- row
    }
  }

  write_csv(do.call(rbind, aligned_rows), "00_站点月份对齐数据_缺月不插值.csv")

  audit <- data.frame(item = c("environment_source", "abundance_source", "taxonomy_source", "samples",
    "large_samples", "small_samples", "sites", "observed_months", "missing_calendar_months", "cyanobacteria_ASVs",
    "rarefaction_depth", "new_variables_requested", "new_variables_available", "new_variables_unavailable",
    "analysis_level", "transform", "lag_window", "permutations", "significance_rule", "multiple_testing_correction"),
    value = c(normalizePath(INPUT[["environment"]], winslash = "/"), normalizePath(INPUT[["abundance"]],
      winslash = "/"), normalizePath(INPUT[["taxonomy"]], winslash = "/"), nrow(metadata), sum(metadata$size ==
      "large"), sum(metadata$size == "small"), length(all_sites), length(observed_dates), length(missing_dates),
      sum(cyan_index), 21000, paste(REQUESTED_VARIABLES, collapse = ";"), paste(new_variables, collapse = ";"),
      paste(REQUESTED_VARIABLES[is.na(matched_columns)], collapse = ";"), "original level only", paste0("positive variables: log10(value + half-minimum-positive); ",
        "Cyanobacteria: log10(value + 0.5/21000); AIRT: original degrees Celsius (signed scale)"),
      "-3 to +3 calendar months", N_PERM, "max-lag permutation p<0.05 AND >=3 valid sites AND >=3 agreeing sites",
      "none (user requested unadjusted p)"))

  write_csv(audit, "00_输入与方法审计.csv")

  variable_label <- c(AIRT = "Air temperature", TNTP = "TN:TP mass ratio", WT = "Water temperature", TOC = "TOC",
    WL = "Water level", Chla = "Chla", CODmn = "CODmn", BOD5 = "BOD5", TN = "TN", TP = "TP", Cyanobacteria = "Cyanobacteria")

  group_label <- c(environment = "environment", large = "large", small = "small", overall = "overall")

  draw_peak_symbol <- function(x, y, status_col) {
    points(x, y, pch = 21, cex = 1.5, bg = status_col, col = "#424242", lwd = 1.2)
  }

  plot_new_variable <- function(new_var) {
    ids <- scenarios$comparison_id[scenarios$var1 == new_var | scenarios$var2 == new_var]
    ids <- ids[!vapply(ids, function(id) {
      s <- scenarios[scenarios$comparison_id == id, ]
      all(c(s$var1, s$var2) %in% new_variables) && new_var != s$var1
    }, logical(1))]
    if (!length(ids))
      return(invisible(NULL))
    ncol_plot <- 2L
    nrow_plot <- ceiling(length(ids)/ncol_plot)
    png(file.path(FIGURE_DIR, paste0("01_", new_var, "参与的原始水平CCF曲线.png")), width = 2000,
      height = max(1300, 620 * nrow_plot), res = 180)
    oldpar <- par(mfrow = c(nrow_plot, ncol_plot), mar = c(4.3, 4.4, 4.2, 1.2), oma = c(1.5, 1.5, 3.2,
      0.5))
    on.exit({
      par(oldpar)
      dev.off()
    }, add = TRUE)
    for (id in ids) {
      s <- scenarios[scenarios$comparison_id == id, ]
      z <- all_lags[all_lags$comparison_id == id, ]
      p <- peaks[peaks$comparison_id == id, ]
      if (s$var1 == new_var) {
        oriented_lag <- z$var1_lead_months
        peak_lag <- p$peak_var1_lead_months
        other <- s$var2
      }
      else {
        oriented_lag <- -z$var1_lead_months
        peak_lag <- -p$peak_var1_lead_months
        other <- s$var1
      }
      ord <- order(oriented_lag)
      col <- if (isTRUE(p$significant_uncorrected))
        "#009E73"
      else "#9E9E9E"
      plot(oriented_lag[ord], z$correlation[ord], type = "o", pch = 16, lwd = 2, col = col, ylim = c(-1,
        1), xaxt = "n", xlab = paste0(new_var, " lead (months)"), ylab = "Pearson CCF (original levels)",
        main = paste0(variable_label[[new_var]], " vs ", variable_label[[other]], " | ", group_label[[s$analysis_group]]))
      axis(1, at = LAGS)
      abline(h = 0, v = 0, lty = 3, col = "grey65")
      draw_peak_symbol(peak_lag, p$peak_correlation, col)
      timing <- if (peak_lag > 0) {
        paste0(new_var, " leads by ", peak_lag, " mo")
      }
      else if (peak_lag < 0) {
        paste0(new_var, " lags by ", -peak_lag, " mo")
      }
      else "synchronous"
      mtext(sprintf("peak %s; r=%+.2f; p=%.4f; sites=%d/%d; %s", timing, p$peak_correlation, p$max_lag_permutation_p_unadjusted,
        p$agreeing_sites, p$valid_sites, ifelse(p$significant_uncorrected, "significant", "not significant")),
        side = 3, line = 0.25, cex = 0.72)
    }
    if (length(ids)%%ncol_plot)
      plot.new()
    mtext(paste0(new_var, " supplementary original-level CCF (unadjusted permutation p)"), outer = TRUE,
      side = 3, line = 1.5, cex = 1.25, font = 2)
    mtext("Green: significant (p<0.05 + site agreement); grey: not significant", outer = TRUE, side = 1,
      line = 0.15, cex = 0.78)
  }

  for (v in new_variables) plot_new_variable(v)

  if (nrow(cyan_peaks)) {
    cyan_plot <- cyan_peaks
    cyan_plot$new_variable <- cyan_plot$var1
    cyan_plot$new_variable_lead_months <- cyan_plot$peak_var1_lead_months
    cyan_plot$label <- paste(cyan_plot$new_variable, cyan_plot$analysis_group, sep = " | ")
    cyan_plot <- cyan_plot[order(match(cyan_plot$new_variable, new_variables), match(cyan_plot$analysis_group,
      c("large", "small", "overall"))), ]
    y <- rev(seq_len(nrow(cyan_plot)))
    cols <- ifelse(cyan_plot$peak_correlation > 0, "#00897B", "#D32F2F")
    png(file.path(FIGURE_DIR, "02_新增变量与蓝藻门_提前滞后汇总.png"), width = 2200, height = max(1200,
      180 * nrow(cyan_plot) + 600), res = 200)
    par(mar = c(5.5, 9.5, 4.5, 6.5))
    plot(c(-3.6, 5.3), c(0.5, nrow(cyan_plot) + 0.5), type = "n", axes = FALSE, xlab = "Signed timing (months): negative = new variable lags; positive = new variable leads",
      ylab = "", main = paste0(paste(unname(variable_label[new_variables]), collapse = " / "), " vs Cyanobacteria: original-level CCF"))
    axis(1, at = LAGS)
    axis(2, at = y, labels = cyan_plot$label, las = 1)
    abline(v = 0, lty = 3, col = "grey55")
    segments(0, y, cyan_plot$new_variable_lead_months, y, col = cols, lwd = 3)
    points(cyan_plot$new_variable_lead_months, y, pch = 21, bg = cols, col = cols, cex = 1.3)
    text(3.25, y, labels = sprintf("r=%+.2f; p=%.4f; sites=%d/%d", cyan_plot$peak_correlation, cyan_plot$max_lag_permutation_p_unadjusted,
      cyan_plot$agreeing_sites, cyan_plot$valid_sites), adj = 0, col = cols, cex = 0.84)
    mtext("Green: r>0 same-direction; red: r<0 opposite-direction", side = 3, line = 0.35)
    dev.off()
  }

  summary_lines <- c(paste0(paste(unname(REQUESTED_CN[REQUESTED_VARIABLES]), collapse = "与"), "补充原始水平CCF摘要"),
    paste0("生成时间: ", Sys.time()), "", paste0("样本: ", nrow(metadata), "（large=", sum(metadata$size ==
      "large"), ", small=", sum(metadata$size == "small"), "）；4站点×15个观测月份；缺月不插值。"),
    paste0("本次实际可分析的新变量: ", paste(new_variables, collapse = "、"), "。"), if (length(unavailable_variables)) {
      paste0("未获得变量: ", paste(unavailable_variables, collapse = "、"), "；未使用其他环境指标替代。")
    } else "全部请求变量均已从统一环境表读取并纳入。", paste0("新增比较数: ",
      nrow(peaks), "；未校正p<0.05且站点一致的显著关系: ", sum(peaks$significant_uncorrected),
      "。"), paste0("方法: 原始水平Pearson CCF；正值指标做log10转换，蓝藻门加0.5/21000后log10；",
      "气温保留原始摄氏度（含负值，不能log10）；lag=-3到+3个月；", "1999次站点内循环位移的最大-lag置换检验。"),
    "显著标准: 未校正置换p<0.05，且至少3个有效站点、至少3站点相关方向一致；不做多重检验校正。",
    "正lag表示var1先变化，负lag表示var2先变化；r正负表示同向/反向，不代表时间方向。",
    if (SUPPLEMENT_MODE == "N_P_RATIO") {
      paste0("氮磷比定义为TN/TP质量比；摩尔比=质量比×30.973761998/14.0067，", "只相差正常数，所以Pearson CCF的r、p和峰值lag完全相同。")
    } else "", if (SUPPLEMENT_MODE == "N_P_RATIO") {
      "TNTP与TN、TP的比较存在组成型数学耦合，仅作完整展示，不能视为独立生态证据。"
    } else "", "原始水平容易受共同季节性、趋势和自相关影响，CCF只能表示时间先后相关，不能证明因果。",
    "", "峰值结果:", apply(peaks, 1L, function(z) sprintf("- %s vs %s | %s: %s -> %s, %s个月, r=%+.3f, p=%.4f, sites=%s/%s, significant=%s",
      z[["var1"]], z[["var2"]], z[["analysis_group"]], z[["leading_variable"]], z[["following_variable"]],
      z[["lead_months"]], as.numeric(z[["peak_correlation"]]), as.numeric(z[["max_lag_permutation_p_unadjusted"]]),
      z[["agreeing_sites"]], z[["valid_sites"]], z[["significant_uncorrected"]])))

  writeLines(summary_lines, file.path(PROJECT, paste0(OUTPUT_TAG, "摘要.txt")), useBytes = TRUE)

  cat("Completed supplementary original-level CCF.\n")

  cat("Mode:", SUPPLEMENT_MODE, "\n")

  cat("Available new variables:", paste(new_variables, collapse = ", "), "\n")

  cat("Unavailable requested variables:", paste(unavailable_variables, collapse = ", "), "\n")

  cat("Comparisons:", nrow(peaks), "\n")

  cat("Significant unadjusted:", sum(peaks$significant_uncorrected), "\n")

  cat("Output:", normalizePath(TABLE_DIR, winslash = "/"), "\n")
}

for (mode in c("TOC_WT", "AIR_TEMP", "N_P_RATIO")) {
  calculate_environment_ccf(mode)
}

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

  PROJECT <- file.path(ROOT, "01_分析项目", "28_营养盐叶绿素蓝藻门CCF")

  TABLE_ROOT <- file.path(PROJECT, "结果表格")

  FIGURE_DIR <- Sys.getenv("CCF_SUPPLEMENTARY_SINGLE_DIR", unset = file.path(PROJECT, "结果图片", "单图"))

  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  read_csv <- function(path) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  }

  draw_five_point_star <- function(x, y, x_radius = 0.17, y_radius = 0.078, fill = "#F0E442", border = "#7A5C00") {
    angles <- pi/2 + seq(0, 2 * pi, length.out = 11)[1:10]
    radii <- rep(c(1, 0.4), 5)
    polygon(x + x_radius * radii * cos(angles), y + y_radius * radii * sin(angles), col = fill, border = border,
      lwd = 1.35)
  }

  load_panel <- function(subdir, variable, title, xlab, ylab, filename, title_cex = 1.1, other_variable = "Cyanobacteria") {
    table_dir <- file.path(TABLE_ROOT, subdir)
    lag_data <- read_csv(file.path(table_dir, "01_新增变量原始水平_CCF全部lag.csv"))
    peak_data <- read_csv(file.path(table_dir, "02_新增变量原始水平_CCF峰值_未校正p.csv"))
    z <- lag_data[lag_data$var1 == variable & lag_data$var2 == other_variable & lag_data$analysis_group %in%
      c("overall", "environment"), , drop = FALSE]
    p <- peak_data[peak_data$var1 == variable & peak_data$var2 == other_variable & peak_data$analysis_group %in%
      c("overall", "environment"), , drop = FALSE]
    if (nrow(z) != 7L || nrow(p) != 1L) {
      stop("Expected seven lag rows and one peak row for ", variable, " vs ", other_variable)
    }
    ord <- order(z$var1_lead_months)
    list(lag = as.numeric(z$var1_lead_months[ord]), correlation = as.numeric(z$correlation[ord]), peak_lag = as.numeric(p$peak_var1_lead_months),
      peak_r = as.numeric(p$peak_correlation), p_value = as.numeric(p$max_lag_permutation_p_unadjusted),
      significant = isTRUE(as.logical(p$significant_uncorrected)), title = title, xlab = xlab, ylab = ylab,
      filename = filename, title_cex = title_cex)
  }

  panels <- list(load_panel("补充_TOC_WT_原始水平", "WT", "Water temperature leads Cyanobacteria",
    "Water temperature lead (months)", "CCF of log10 levels", "06g_Water_temperature_leads_Cyanobacteria.png",
    title_cex = 1.22), load_panel("补充_气温AIRT_原始水平", "AIRT", "Air temperature leads Cyanobacteria",
    "Air temperature lead (months)", "CCF of original levels", "06h_Air_temperature_leads_Cyanobacteria.png",
    title_cex = 1.1), load_panel("补充_氮磷比TNTP_原始水平", "TNTP", "TN:TP mass ratio leads Chla",
    "TN:TP mass ratio lead (months)", "CCF of log10 levels", "06i_TNTP_mass_ratio_leads_Chla.png", title_cex = 1.22,
    other_variable = "Chla"))

  draw_single_panel <- function(d) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(4.8, 5, 1.2, 1), mgp = c(2.9, 0.78, 0), tcl = -0.28, family = "sans", las = 1)
    plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
      xlab = d$xlab, ylab = d$ylab, main = "", bty = "o", cex.lab = 1.25, cex.axis = 1.1)
    axis(1, at = -3:3, labels = -3:3)
    axis(2, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)))
    box(lwd = 1)
    abline(h = 0, v = 0, lty = 2, lwd = 0.8, col = "grey72")
    segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.7, col = "grey82")
    lines(d$lag, d$correlation, lty = 1, lwd = 2.1, col = "#00796B")
    points(d$lag, d$correlation, pch = 21, cex = 1.25, bg = "#00796B", col = "white", lwd = 0.85)
    if (d$significant) {
      segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.55, col = "#D73027")
      draw_five_point_star(d$peak_lag, d$peak_r)
    }
    text(-3.08, 0.94, labels = d$title, adj = c(0, 1), font = 2, cex = d$title_cex)
    text(-3.08, 0.78, labels = sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), adj = c(0, 1), cex = 1.15)
  }

  for (d in panels) {
    png(file.path(FIGURE_DIR, d$filename), width = 2400, height = 1800, res = 400, type = "cairo", bg = "white")
    draw_single_panel(d)
    dev.off()
    grDevices::cairo_pdf(file.path(FIGURE_DIR, sub("\\.png$", ".pdf", d$filename)), width = 6, height = 4.5,
      bg = "white")
    draw_single_panel(d)
    dev.off()
    svg(file.path(FIGURE_DIR, sub("\\.png$", ".svg", d$filename)), width = 6, height = 4.5, bg = "white")
    draw_single_panel(d)
    dev.off()
  }

  message("Styled supplementary CCF single panels saved to: ", FIGURE_DIR)
})
