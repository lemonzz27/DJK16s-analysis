options(stringsAsFactors = FALSE)

Sys.setlocale("LC_ALL", "English_United States.utf8")

root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

out <- file.path(root, "work/supplement_revision/data")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

.libPaths(unique(c(file.path(root, "99_运行环境/R_library"), .libPaths())))

save_csv <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE, fileEncoding = "UTF-8")

src <- file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/06_community_composition/05_plot_monthly_composition.R")

menv <- new.env()

for (e in parse(src, encoding = "UTF-8")) {
  eval(e, menv)
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("prepare_monthly_level")))
    break
}

menv$top_n <- 1000000L

for (lev in c("phylum", "order", "genus")) {
  z <- menv$prepare_monthly_level(menv$dataset_definitions$SILVA_all, lev)
  save_csv(z$data, paste0("monthly_all_", lev, ".csv"))
}

gsrc <- file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/06_community_composition/01_mbpd_composition_functions.R")

genv <- new.env(parent = globalenv())

for (e in parse(gsrc, encoding = "UTF-8")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
    eval(e, genv)
}

genv$top_n <- 1000000L

genv$category_levels <- "Size"

genv$group_levels <- c("Large", "Small")

genv$group_contract <- list(Size = list(order = c("large", "small"), labels = c(large = "Large", small = "Small")))

genv$group_metadata <- list(Size = data.frame(SampleID = menv$size_meta$SampleID, GroupRaw = menv$size_meta$SizeRaw))

mat <- menv$abundance_all[menv$candidate_ids, , drop = FALSE]

ab <- list(ASV_ID = rownames(mat), matrix = mat, n_asv = nrow(mat), n_samples = ncol(mat))

raw_candidate <- read.csv(file.path(root, "00_输入文件/05_致病菌候选/候选致病菌原始来源表_04_pathogen_candidate_abundance.csv"),
  row.names = 1, check.names = FALSE)

genv$abundance_matrix <- as.matrix(raw_candidate)

genv$taxonomy <- read.csv(file.path(root, "01_分析项目/09_物种组成/结果表格/MBPD/致病菌_MBPD七层级分组组成/MBPD_ASV_taxonomy_LCA_2755.csv"),
  check.names = FALSE)

z <- genv$prepare_level_data("genus")

save_csv(z$data, "candidate_all_genera.csv")

renv <- new.env()

rsrc <- file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/06_community_composition/00_taxonomy_network_functions.R")

owd <- getwd()

setwd(file.path(root, "work/taxonomic_network_palette"))

for (e in parse(rsrc, encoding = "UTF-8")) {
  eval(e, renv)
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("phylum_abundance")))
    break
}

setwd(owd)

save_csv(renv$phylum_abundance, "radial_all_phyla.csv")

alpha <- read.csv(file.path(root, "01_分析项目/03_Alpha多样性/结果表格/基础_注释无关/整体_Alpha计算/alpha1.csv"),
  row.names = 1)

gr <- read.csv(file.path(root, "00_输入文件/02_样本分组与站点/孔径分组_groupsize.csv"))

alpha$group <- factor(gr$group[match(rownames(alpha), gr$sample)], levels = c("large", "small"))

res <- lapply(c("Richness", "Shannon", "Simpson", "Chao1"), function(metric) {
  d <- data.frame(Value = alpha[[metric]], group = alpha$group)
  a <- summary(aov(Value ~ group, data = d))[[1]]
  data.frame(Metric = metric, n_ConF = sum(d$group == "large"), Mean_ConF = mean(d$Value[d$group ==
    "large"]), SD_ConF = sd(d$Value[d$group == "large"]), n_USF = sum(d$group == "small"), Mean_USF = mean(d$Value[d$group ==
    "small"]), SD_USF = sd(d$Value[d$group == "small"]), F = a[1, "F value"], df1 = a[1, "Df"], df2 = a[2,
    "Df"], P = a[1, "Pr(>F)"])
})

save_csv(do.call(rbind, res), "alpha_summary_anova.csv")

cat("Complete exports finished\n")
