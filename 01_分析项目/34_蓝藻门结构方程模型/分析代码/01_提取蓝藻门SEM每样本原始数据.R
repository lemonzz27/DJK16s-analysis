options(stringsAsFactors = FALSE, scipen = 999)

locate_script_dir <- function() {
  source_files <- vapply(sys.frames(), function(x) {
    if (is.null(x$ofile))
      ""
    else as.character(x$ofile)
  }, character(1))
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files) > 0) {
    return(dirname(normalizePath(tail(source_files, 1), winslash = "/", mustWork = TRUE)))
  }
  script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_argument) > 0) {
    script_file <- sub("^--file=", "", script_argument[1])
    return(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- locate_script_dir()

config_file <- file.path(dirname(script_dir), "配置", "SEM_config.R")

if (!file.exists(config_file)) stop("找不到配置文件：", config_file)

source(config_file, encoding = "UTF-8")

OUTPUT_DIR <- file.path(PROJECT, "原始数据")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

INPUT <- c(abundance = file.path(ROOT, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv"),
  taxonomy = file.path(ROOT, "00_输入文件", "08_派生输入", "分类注释标准化", "SILVA",
    "tax.csv"), environment = file.path(ROOT, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv"),
  location = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "地点分组_grouplocation.csv"),
  month = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "月份分组_groupmonth.csv"),
  season = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "季节分组_groupseason.csv"),
  rain = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "降雨分组_grouprain.csv"),
  size = file.path(ROOT, "00_输入文件", "02_样本分组与站点", "孔径分组_groupsize.csv"))

if (!all(file.exists(INPUT))) {
  stop("以下输入文件不存在：", paste(names(INPUT)[!file.exists(INPUT)], collapse = ", "))
}

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(x, filename) {
  write.csv(x, file.path(OUTPUT_DIR, filename), row.names = FALSE, fileEncoding = "UTF-8")
}

read_group <- function(path, sample_ids, group_name) {
  group_data <- read_csv(path)
  if (!all(c("sample", "group") %in% names(group_data))) {
    stop(group_name, " 分组文件必须包含 sample 和 group 两列。")
  }
  values <- as.character(group_data$group[match(sample_ids, group_data$sample)])
  if (anyNA(values)) {
    stop(group_name, " 分组存在未匹配样本。")
  }
  values
}

abundance_table <- read_csv(INPUT[["abundance"]])

taxonomy <- read_csv(INPUT[["taxonomy"]])

environment <- read_csv(INPUT[["environment"]])

if (!"mags" %in% names(abundance_table) || !"mags" %in% names(taxonomy)) {
  stop("丰度表与分类表必须都包含 mags 列。")
}

if (!"sample" %in% names(environment)) {
  stop("环境表必须包含 sample 列。")
}

sample_ids <- names(abundance_table)[-1]

if (length(sample_ids) != 120L) {
  warning("丰度表样本数为 ", length(sample_ids), "，项目预期为 120。")
}

if (anyDuplicated(sample_ids)) stop("丰度表存在重复样本名。")

abundance <- as.matrix(abundance_table[, -1, drop = FALSE])

storage.mode(abundance) <- "numeric"

row.names(abundance) <- abundance_table$mags

total_reads <- colSums(abundance)

if (any(total_reads != 21000)) {
  stop("并非所有样本都抽平到 21000 reads。")
}

taxonomy <- taxonomy[match(row.names(abundance), taxonomy$mags), , drop = FALSE]

if (anyNA(taxonomy$mags)) {
  stop("部分 ASV 在 SILVA 标准化分类表中未匹配。")
}

phylum <- trimws(as.character(taxonomy$phylum))

cyanobacteria_index <- grepl("^(Cyanobacteriota|Cyanobacteria)$", phylum, ignore.case = TRUE)

if (!any(cyanobacteria_index)) {
  stop("SILVA 分类中未找到 Cyanobacteriota/Cyanobacteria。")
}

cyanobacteria_reads <- colSums(abundance[cyanobacteria_index, , drop = FALSE])

cyanobacteria_asv_detected <- colSums(abundance[cyanobacteria_index, , drop = FALSE] > 0)

cyanobacteria_relative <- cyanobacteria_reads/total_reads

site <- read_group(INPUT[["location"]], sample_ids, "地点")

month_code <- sprintf("%04d", as.integer(read_group(INPUT[["month"]], sample_ids, "月份")))

season <- read_group(INPUT[["season"]], sample_ids, "季节")

rain_group <- read_group(INPUT[["rain"]], sample_ids, "降雨")

size_fraction <- tolower(read_group(INPUT[["size"]], sample_ids, "孔径"))

if (!all(size_fraction %in% c("large", "small"))) {
  stop("孔径分组只能是 large 或 small。")
}

sample_date <- as.Date(paste0("20", substr(month_code, 1, 2), "-", substr(month_code, 3, 4), "-01"))

site_month_id <- paste(site, format(sample_date, "%Y-%m"), sep = "_")

environment_matched <- environment[match(sample_ids, environment$sample), setdiff(names(environment),
  "sample"), drop = FALSE]

if (anyNA(match(sample_ids, environment$sample))) {
  stop("环境表存在未匹配样本。")
}

sample_data <- data.frame(sample = sample_ids, site = site, month_code = month_code, date = sample_date,
  site_month_id = site_month_id, season = season, rain_group = rain_group, size_fraction = size_fraction,
  environment_matched, total_reads = as.numeric(total_reads[sample_ids]), cyanobacteria_reads = as.numeric(cyanobacteria_reads[sample_ids]),
  cyanobacteria_ASVs_detected = as.numeric(cyanobacteria_asv_detected[sample_ids]), cyanobacteria_relative = as.numeric(cyanobacteria_relative[sample_ids]),
  cyanobacteria_percent = 100 * as.numeric(cyanobacteria_relative[sample_ids]), taxonomy_source = "SILVA standardized phylum=Cyanobacteriota/Cyanobacteria",
  stringsAsFactors = FALSE, check.names = FALSE)

sample_data$TN_TP_mass_ratio <- with(sample_data, ifelse(TP > 0, TN/TP, NA_real_))

sample_data$TN_TP_molar_ratio <- sample_data$TN_TP_mass_ratio * 30.973761998/14.0067

preferred_order <- c("sample", "site", "month_code", "date", "site_month_id", "season", "rain_group",
  "size_fraction", "Longitude", "Latitude", "WT", "pH", "DO", "CODmn", "BOD5", "NH3N", "TN", "TP",
  "TN_TP_mass_ratio", "TN_TP_molar_ratio", "Cu", "Zn", "F", "Se", "As", "Pb", "Chla", "Q", "WL", "total_reads",
  "cyanobacteria_reads", "cyanobacteria_ASVs_detected", "cyanobacteria_relative", "cyanobacteria_percent",
  "taxonomy_source")

preferred_order <- preferred_order[preferred_order %in% names(sample_data)]

sample_data <- sample_data[, c(preferred_order, setdiff(names(sample_data), preferred_order)), drop = FALSE]

sample_data <- sample_data[order(sample_data$site, sample_data$date, sample_data$size_fraction), , drop = FALSE]

row.names(sample_data) <- NULL

environment_columns <- setdiff(names(environment), "sample")

environment_range_audit <- do.call(rbind, lapply(split(sample_data, sample_data$site_month_id), function(z) {
  ranges <- vapply(environment_columns, function(v) {
    x <- as.numeric(z[[v]])
    diff(range(x, na.rm = TRUE))
  }, numeric(1))
  data.frame(site_month_id = z$site_month_id[1], n_samples = nrow(z), max_environment_range = max(ranges,
    na.rm = TRUE), stringsAsFactors = FALSE)
}))

if (any(environment_range_audit$n_samples != 2L)) {
  stop("部分站点-月份不是恰好两个孔径样本。")
}

if (any(environment_range_audit$max_environment_range != 0)) {
  stop("同一站点-月份的 large/small 环境值并非完全一致。")
}

write_csv(sample_data, "蓝藻门SEM_每样本原始数据_120样本.csv")

write_csv(environment_range_audit, "蓝藻门SEM_环境值跨孔径一致性审计.csv")

site_month_data <- do.call(rbind, lapply(split(sample_data, sample_data$site_month_id), function(z) {
  large <- z[z$size_fraction == "large", , drop = FALSE]
  small <- z[z$size_fraction == "small", , drop = FALSE]
  if (nrow(large) != 1L || nrow(small) != 1L) {
    stop("站点-月份缺少 large 或 small 样本：", z$site_month_id[1])
  }
  out <- z[1, c("site", "month_code", "date", "site_month_id", "season", "rain_group", environment_columns,
    "TN_TP_mass_ratio", "TN_TP_molar_ratio"), drop = FALSE]
  out$cyanobacteria_reads_large <- large$cyanobacteria_reads
  out$cyanobacteria_reads_small <- small$cyanobacteria_reads
  out$cyanobacteria_relative_large <- large$cyanobacteria_relative
  out$cyanobacteria_relative_small <- small$cyanobacteria_relative
  out$cyanobacteria_relative_overall <- mean(z$cyanobacteria_relative)
  out$cyanobacteria_percent_large <- large$cyanobacteria_percent
  out$cyanobacteria_percent_small <- small$cyanobacteria_percent
  out$cyanobacteria_percent_overall <- mean(z$cyanobacteria_percent)
  out
}))

site_month_data <- site_month_data[order(site_month_data$site, as.Date(site_month_data$date)), , drop = FALSE]

row.names(site_month_data) <- NULL

write_csv(site_month_data, "蓝藻门SEM_站点月份数据_60观测.csv")

source_audit <- data.frame(Item = c("abundance_source", "taxonomy_source", "environment_source", "sample_rows",
  "site_month_rows", "sites", "months_per_site", "rarefaction_depth", "cyanobacteria_ASVs_in_SILVA",
  "large_samples", "small_samples", "environment_duplicate_check", "TN_TP_mass_ratio", "TN_TP_molar_ratio"),
  Value = c(INPUT[["abundance"]], INPUT[["taxonomy"]], INPUT[["environment"]], nrow(sample_data), nrow(site_month_data),
    length(unique(sample_data$site)), paste(table(sample_data$site)/2, collapse = ";"), unique(sample_data$total_reads),
    sum(cyanobacteria_index), sum(sample_data$size_fraction == "large"), sum(sample_data$size_fraction ==
      "small"), "passed: all paired environment ranges are zero", "TN / TP (mass ratio)", "(TN / TP) * 30.973761998 / 14.0067"))

write_csv(source_audit, "蓝藻门SEM_数据来源与计算审计.csv")

cat("数据提取完成：\n")

cat("- 每样本表：", nrow(sample_data), " 行\n")

cat("- 站点-月份表：", nrow(site_month_data), " 行\n")

cat("- SILVA 蓝藻门 ASV：", sum(cyanobacteria_index), " 个\n")

cat("输出目录：", normalizePath(OUTPUT_DIR, winslash = "/"), "\n")
