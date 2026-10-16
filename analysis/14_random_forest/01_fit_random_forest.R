Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")

.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(randomForest)
  library(rfPermute)
  library(dplyr)
  library(scales)
})

ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

OUT <- file.path(ROOT, "outputs/20261002_manuscript_revision/RF")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

input_root <- file.path(ROOT, "99_运行环境/串行运行工作区/SILVA")

ENV_FILE <- file.path(input_root, "env.csv")

PCOA_FILE <- file.path(input_root, "PCOA1年际变化/1.csv")

GROUP_FILE <- file.path(input_root, "随机森林/groupsize.csv")

SOURCE <- file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/14_random_forest/00_reference_random_forest.R")

for (f in c(ENV_FILE, PCOA_FILE, GROUP_FILE, SOURCE)) stopifnot(file.exists(f))

write.csv(data.frame(role = c("environment", "response", "group", "predecessor_script"), path = c(ENV_FILE,
  PCOA_FILE, GROUP_FILE, SOURCE), md5 = unname(tools::md5sum(c(ENV_FILE, PCOA_FILE, GROUP_FILE, SOURCE)))),
  file.path(OUT, "input_provenance.csv"), row.names = FALSE)

env_data <- read.csv(ENV_FILE, stringsAsFactors = FALSE, check.names = FALSE)

pcoa_data <- read.csv(PCOA_FILE, stringsAsFactors = FALSE, check.names = FALSE)[, c("sample", "PCoA1")]

group_data <- read.csv(GROUP_FILE, stringsAsFactors = FALSE, check.names = FALSE)[, c("sample", "group")]

stopifnot("Q" %in% names(env_data), !anyDuplicated(pcoa_data$sample), !anyDuplicated(group_data$sample))

env_data$Q <- NULL

write.csv(env_data, file.path(OUT, "environment_Q_removed.csv"), row.names = FALSE)

dat <- inner_join(env_data, pcoa_data, by = "sample") %>%
  left_join(group_data, by = "sample")

stopifnot(!anyNA(dat$group))

dat$group <- trimws(as.character(dat$group))

dat$RF_Group <- case_when(tolower(dat$group) %in% c("large", "大型") ~ "ConF", tolower(dat$group) %in%
  c("small", "小型") ~ "USF", TRUE ~ dat$group)

response_col <- "PCoA1"

predictor_cols <- setdiff(names(dat), c("sample", "group", "RF_Group", response_col, "Longitude", "Latitude"))

predictor_cols <- predictor_cols[sapply(dat[predictor_cols], is.numeric)]

stopifnot(length(predictor_cols) == 19L, !"Q" %in% predictor_cols, nrow(dat) == 120L, sum(dat$RF_Group ==
  "ConF") == 60L, sum(dat$RF_Group == "USF") == 60L)

dat[[response_col]] <- as.numeric(dat[[response_col]])

write.csv(dat[, c("sample", "RF_Group", response_col, predictor_cols)], file.path(OUT, "model_input_120_samples_19_predictors.csv"),
  row.names = FALSE)

write.csv(data.frame(predictor = predictor_cols), file.path(OUT, "predictors.csv"), row.names = FALSE)

sig_label <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p <
  0.05, "*", ""))))

extract_pvalue <- function(rf_obj, vars, measure = "%IncMSE") {
  pval_obj <- rf_obj$pval
  p_out <- setNames(rep(NA_real_, length(vars)), vars)
  if (is.null(pval_obj))
    return(p_out)
  pmat <- if (length(dim(pval_obj)) == 3)
    pval_obj[, , 2, drop = TRUE]
  else pval_obj
  pmat <- as.data.frame(pmat, check.names = FALSE)
  use_measure <- if (measure %in% colnames(pmat))
    measure
  else colnames(pmat)[1]
  common_vars <- intersect(vars, rownames(pmat))
  p_out[common_vars] <- pmat[common_vars, use_measure]
  p_out
}

N_TREE <- 500L

N_REP <- 999L

run_rf_one_group <- function(input_data, group_name) {
  model_data <- input_data %>%
    select(all_of(c(response_col, predictor_cols))) %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  stopifnot(nrow(model_data) == nrow(input_data))
  message("Starting ", group_name, ": ", nrow(model_data), " samples; ", Sys.time())
  rng_before <- .Random.seed
  rf_fit <- rfPermute(as.formula(paste(response_col, "~ .")), data = model_data, ntree = N_TREE, num.rep = N_REP,
    importance = TRUE)
  saveRDS(list(fit = rf_fit, rng_before = rng_before, rng_after = .Random.seed, model_data = model_data),
    file.path(OUT, paste0("RF_model_", group_name, ".rds")))
  imp <- data.frame(importance(rf_fit, scale = TRUE), check.names = FALSE)
  imp$Variable <- rownames(imp)
  pvals <- extract_pvalue(rf_fit, imp$Variable, "%IncMSE")
  imp$P_value <- pvals[imp$Variable]
  imp$P_adj_BH <- p.adjust(imp$P_value, "BH")
  imp$Sig <- sig_label(imp$P_adj_BH)
  imp$Group <- group_name
  imp$N_sample <- nrow(model_data)
  imp <- arrange(imp, desc(`%IncMSE`))
  rf_core <- if (!is.null(rf_fit$rf))
    rf_fit$rf
  else rf_fit
  summary <- data.frame(Group = group_name, N_sample = nrow(model_data), N_predictor = length(predictor_cols),
    MSE = tail(rf_core$mse, 1), Variance_explained_percent = tail(rf_core$rsq, 1) * 100)
  write.csv(imp, file.path(OUT, paste0("RF_importance_", group_name, ".csv")), row.names = FALSE)
  write.csv(summary, file.path(OUT, paste0("RF_summary_", group_name, ".csv")), row.names = FALSE)
  message("Finished ", group_name, ": ", Sys.time())
  list(importance = imp, summary = summary)
}

set.seed(123)

rf_all <- run_rf_one_group(dat, "All")

rf_conf <- run_rf_one_group(filter(dat, RF_Group == "ConF"), "ConF")

rf_usf <- run_rf_one_group(filter(dat, RF_Group == "USF"), "USF")

importance_all <- bind_rows(rf_all$importance, rf_conf$importance, rf_usf$importance) %>%
  mutate(Group = factor(Group, levels = c("All", "ConF", "USF")), Label = sprintf("%.2f %s", `%IncMSE`,
    Sig))

summary_all <- bind_rows(rf_all$summary, rf_conf$summary, rf_usf$summary)

write.csv(importance_all, file.path(OUT, "RF_importance_All_ConF_USF.csv"), row.names = FALSE)

write.csv(importance_all, file.path(OUT, "Table_S15_RF_importance_replacement.csv"), row.names = FALSE)

write.csv(summary_all, file.path(OUT, "RF_model_summary.csv"), row.names = FALSE)

write.csv(summary_all, file.path(OUT, "Table_S15_RF_model_summary_replacement.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

source(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/14_random_forest/02_plot_random_forest.R"), encoding = "UTF-8")
