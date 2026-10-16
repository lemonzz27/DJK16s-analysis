Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")

.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "Rlib"), .libPaths()))

suppressPackageStartupMessages(library(vegan))

project_dir <- normalizePath(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

analysis_dir <- file.path(project_dir, "01_分析项目", "33_VPA时间环境空间")

table_dir <- file.path(analysis_dir, "结果表格", "基础_注释无关", "整体_时间环境空间_VPA")

figure_dir <- file.path(project_dir, "outputs", "20261002_manuscript_revision", "S9")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

asv_raw <- read.csv(file.path(project_dir, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv"),
  check.names = FALSE)

size_raw <- read.csv(file.path(project_dir, "00_输入文件", "02_样本分组与站点", "孔径分组_groupsize.csv"),
  check.names = FALSE)

month_raw <- read.csv(file.path(project_dir, "00_输入文件", "02_样本分组与站点", "月份分组_groupmonth.csv"),
  check.names = FALSE)

env_raw <- read.csv(file.path(project_dir, "00_输入文件", "03_环境数据", "环境因子_分析用_env.csv"),
  check.names = FALSE)

sample_ids <- names(asv_raw)[-1]

site_from_id <- function(x) ifelse(startsWith(x, "LHKX"), "LHKX", sub("[0-9].*$", "", x))

meta <- data.frame(sample = sample_ids, size = factor(tolower(size_raw$group[match(sample_ids, size_raw$sample)]),
  levels = c("small", "large")), month = factor(as.character(month_raw$group[match(sample_ids, month_raw$sample)])),
  site = factor(site_from_id(sample_ids), levels = c("KX", "TC", "BS", "LHKX")))

asv_num <- as.matrix(asv_raw[, -1, drop = FALSE])

storage.mode(asv_num) <- "double"

asv_num <- asv_num[rowSums(asv_num) > 0.1, , drop = FALSE]

comm_all <- t(asv_num)

rownames(comm_all) <- sample_ids

selected <- read.csv(file.path(table_dir, "环境变量纳入与剔除清单.csv"), check.names = FALSE)

selected <- selected$variable[selected$category == "selected"]

env_aligned <- env_raw[match(sample_ids, env_raw$sample), , drop = FALSE]

env_scaled <- as.data.frame(scale(env_aligned[, selected, drop = FALSE]), check.names = FALSE)

make_design <- function(meta_sub, env_sub) {
  time_mm <- model.matrix(~month, data = meta_sub)[, -1, drop = FALSE]
  space_mm <- model.matrix(~site, data = meta_sub)[, -1, drop = FALSE]
  list(Time = as.data.frame(time_mm), Environment = as.data.frame(env_sub), Space = as.data.frame(space_mm))
}

fractions <- read.csv(file.path(table_dir, "VPA_全部调整R2分量.csv"), check.names = FALSE)

groups <- list(All = rep(TRUE, nrow(meta)), Large = meta$size == "large", Small = meta$size == "small")

vp_objects <- list(dbRDA = list())

for (g in names(groups)) {
  idx <- groups[[g]]
  X <- make_design(droplevels(meta[idx, , drop = FALSE]), env_scaled[idx, , drop = FALSE])
  counts <- comm_all[idx, , drop = FALSE]
  vp_objects$dbRDA[[g]] <- varpart(vegdist(counts, "bray"), X$Time, X$Environment, X$Space, sqrt.dist = TRUE)
}

plot_vpa_panels <- function(vp_list, method_label, file_stub) {
  display_label <- c(All = "All", Large = "ConF", Small = "USF")
  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(1, 3), mar = c(2, 2, 4, 2), oma = c(0, 0, 2, 0))
    for (g in c("All", "Large", "Small")) {
      frac <- fractions[fractions$method == method_label & fractions$group == g, ]
      explained <- 100 * (1 - frac$adjusted_R2[frac$fraction_code == "[h]"])
      plot(vp_list[[g]], Xnames = c("Time", "Env.", "Space"), bg = c("#E69F0088", "#009E7388",
        "#0072B288"), digits = 3, cutoff = -Inf, id.size = 0.9, cex = 1.05)
      title(main = sprintf("%s\nTotal adjusted R² = %.1f%%", display_label[[g]], explained), cex.main = 1.05)
    }
    mtext(paste0(method_label, " variation partitioning"), outer = TRUE, cex = 1.25, font = 2)
  }
  png(file.path(figure_dir, paste0(file_stub, ".png")), width = 3300, height = 1200, res = 260, bg = "white")
  draw()
  dev.off()
  pdf(file.path(figure_dir, paste0(file_stub, ".pdf")), width = 13, height = 5, useDingbats = FALSE)
  draw()
  dev.off()
  svg(file.path(figure_dir, paste0(file_stub, ".svg")), width = 13, height = 5, bg = "white")
  draw()
  dev.off()
}

plot_vpa_panels(vp_objects$dbRDA, "dbRDA-Bray", "VPA_dbRDA_Bray_All_ConF_USF")
