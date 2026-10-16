local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
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

  stopifnot(all(file.exists(c(ENV_FILE, PCOA_FILE, GROUP_FILE))))

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
})

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")

  .libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/R_library"), .libPaths()))

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
    library(scales)
  })

  OUT <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "outputs/20261002_manuscript_revision/RF")

  importance_all <- read.csv(file.path(OUT, "RF_importance_All_ConF_USF.csv"), check.names = FALSE)

  importance_all$Sig[is.na(importance_all$Sig)] <- ""

  importance_all$Group <- factor(importance_all$Group, levels = c("All", "ConF", "USF"))

  importance_all$Label <- sprintf("%.2f %s", importance_all$`%IncMSE`, importance_all$Sig)

  make_plot <- function(groups, base_size = 9) {
    pd <- importance_all %>%
      filter(Group %in% groups) %>%
      group_by(Group) %>%
      arrange(`%IncMSE`, .by_group = TRUE) %>%
      mutate(Variable_facet = paste(Variable, Group, sep = "___")) %>%
      ungroup()
    pd$Variable_facet <- factor(pd$Variable_facet, levels = pd$Variable_facet)
    max_abs <- max(abs(pd$`%IncMSE`))
    pd$Label_y <- pd$`%IncMSE` + ifelse(pd$`%IncMSE` >= 0, 1, -1) * 0.025 * max_abs
    pd$Label_hjust <- ifelse(pd$`%IncMSE` >= 0, 0, 1)
    ggplot(pd, aes(x = Variable_facet, y = `%IncMSE`, fill = `%IncMSE`)) + geom_hline(yintercept = 0,
      linewidth = 0.35, color = "black") + geom_col(width = 0.78, color = "black", linewidth = 0.25) +
      geom_text(aes(y = Label_y, label = Label, hjust = Label_hjust), size = base_size/ggplot2::.pt,
        color = "black", family = "Arial") + scale_y_continuous(expand = expansion(mult = c(0.15,
      0.22))) + coord_flip(clip = "off") + facet_wrap(~Group, ncol = 1, scales = "free_y") + scale_x_discrete(labels = function(x) sub("___.*$",
      "", x)) + scale_fill_gradientn(colors = c(brewer.pal(5, "Blues")[2], brewer.pal(5, "YlGnBu")[1],
      brewer.pal(6, "OrRd")[4]), values = rescale(c(min(pd$`%IncMSE`), median(pd$`%IncMSE`), max(pd$`%IncMSE`))),
      guide = "none") + labs(x = NULL, y = "Permutation importance (scaled %IncMSE)") + theme_test(base_size = base_size,
      base_family = "Arial") + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white",
      color = "black", linewidth = 0.4), strip.text = element_text(face = "bold", size = base_size +
      1), axis.text = element_text(color = "black", size = base_size), axis.title.x = element_text(size = base_size +
      1), plot.margin = ggplot2::margin(5, 8, 5, 5))
  }

  export_plot <- function(p, name, width_mm, height_mm) {
    ggsave(file.path(OUT, paste0(name, ".png")), p, width = width_mm, height = height_mm, units = "mm",
      dpi = 600, device = ragg::agg_png, bg = "white")
    ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = width_mm, height = height_mm, units = "mm",
      device = cairo_pdf, bg = "white")
    ggsave(file.path(OUT, paste0(name, ".svg")), p, width = width_mm, height = height_mm, units = "mm",
      device = svglite::svglite, bg = "white")
  }

  export_plot(make_plot(c("ConF", "USF")), "Figure_S10_RF_ConF_USF_without_Q", 170, 204)

  export_plot(make_plot(c("All", "ConF", "USF")), "RF_importance_All_ConF_USF_without_Q", 183, 275)

  for (g in c("All", "ConF", "USF")) export_plot(make_plot(g), paste0("RF_importance_", g, "_without_Q"),
    170, 106)
})
