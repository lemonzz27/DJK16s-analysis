rm(list = ls())

pkgs <- c("ggplot2", "RColorBrewer", "randomForest", "rfPermute", "dplyr", "scales")

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("缺少 R 包：", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(randomForest)
  library(rfPermute)
  library(dplyr)
  library(scales)
})

set.seed(123)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- Sys.getenv("DJK16S_ACTIVE_ROOT", unset = get_script_dir())

project_dir <- dirname(script_dir)

ENV_FILE <- file.path(project_dir, "env.csv")

PCOA_FILE <- file.path(project_dir, "PCOA1年际变化", "1.csv")

GROUP_FILE <- file.path(script_dir, "groupsize.csv")

OUT_DIR <- file.path(script_dir, "outputs", "rf_size_group")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_TREE <- 500

N_REP <- 999

env_data <- read.csv(ENV_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

pcoa_data <- read.csv(PCOA_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

group_data <- read.csv(GROUP_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

env_data <- env_data %>%
  rename(sample = sample)

if (!all(c("sample", "PCoA1") %in% names(pcoa_data))) {
  stop("PCoA 文件必须包含 sample 和 PCoA1 两列。")
}

pcoa_data <- pcoa_data[, c("sample", "PCoA1"), drop = FALSE]

group_data <- group_data %>%
  rename(sample = sample, group = group) %>%
  select(sample, group)

dat <- env_data %>%
  inner_join(pcoa_data, by = "sample") %>%
  left_join(group_data, by = "sample")

if (any(is.na(dat$group))) {
  stop("有样本在 groupsize 文件中找不到 large/small 分组，请检查 sample 名是否一致。")
}

dat$group <- trimws(as.character(dat$group))

dat$RF_Group <- case_when(tolower(dat$group) %in% c("large", "大型") ~ "ConF", tolower(dat$group) %in%
  c("small", "小型") ~ "USF", TRUE ~ dat$group)

response_col <- "PCoA1"

predictor_cols <- setdiff(names(dat), c("sample", "group", "RF_Group", response_col, "Longitude", "Latitude"))

predictor_cols <- predictor_cols[sapply(dat[predictor_cols], is.numeric)]

dat[[response_col]] <- as.numeric(dat[[response_col]])

message("响应变量：", response_col)

message("环境因子：", paste(predictor_cols, collapse = ", "))

message("All 样本数：", nrow(dat))

message("ConF 样本数：", sum(dat$RF_Group == "ConF"))

message("USF 样本数：", sum(dat$RF_Group == "USF"))

sig_label <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}

extract_pvalue <- function(rf_obj, vars, measure = "%IncMSE") {
  pval_obj <- rf_obj$pval
  p_out <- rep(NA_real_, length(vars))
  names(p_out) <- vars
  if (is.null(pval_obj))
    return(p_out)
  if (length(dim(pval_obj)) == 3) {
    pmat <- pval_obj[, , 2, drop = TRUE]
  }
  else {
    pmat <- pval_obj
  }
  pmat <- as.data.frame(pmat, check.names = FALSE)
  use_measure <- if (measure %in% colnames(pmat))
    measure
  else colnames(pmat)[1]
  common_vars <- intersect(vars, rownames(pmat))
  p_out[common_vars] <- pmat[common_vars, use_measure]
  return(p_out)
}

run_rf_one_group <- function(input_data, group_name) {
  model_data <- input_data %>%
    select(all_of(c(response_col, predictor_cols))) %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  rf_formula <- as.formula(paste(response_col, "~ ."))
  rf_fit <- rfPermute(rf_formula, data = model_data, ntree = N_TREE, num.rep = N_REP, importance = TRUE)
  imp <- data.frame(importance(rf_fit, scale = TRUE), check.names = FALSE)
  imp$Variable <- rownames(imp)
  pvals <- extract_pvalue(rf_fit, vars = imp$Variable, measure = "%IncMSE")
  imp$P_value <- pvals[imp$Variable]
  imp$P_adj_BH <- p.adjust(imp$P_value, method = "BH")
  imp$Sig <- sig_label(imp$P_adj_BH)
  imp$Group <- group_name
  imp$N_sample <- nrow(model_data)
  imp <- imp %>%
    arrange(desc(`%IncMSE`))
  rf_core <- if (!is.null(rf_fit$rf))
    rf_fit$rf
  else rf_fit
  summary_df <- data.frame(Group = group_name, N_sample = nrow(model_data), N_predictor = length(predictor_cols),
    MSE = if (!is.null(rf_core$mse))
      tail(rf_core$mse, 1)
    else NA, Variance_explained_percent = if (!is.null(rf_core$rsq))
      tail(rf_core$rsq, 1) * 100
    else NA)
  return(list(importance = imp, summary = summary_df))
}

rf_all <- run_rf_one_group(dat, "All")

rf_large <- run_rf_one_group(dat %>%
  filter(RF_Group == "ConF"), "ConF")

rf_small <- run_rf_one_group(dat %>%
  filter(RF_Group == "USF"), "USF")

importance_all <- bind_rows(rf_all$importance, rf_large$importance, rf_small$importance) %>%
  mutate(Group = factor(Group, levels = c("All", "ConF", "USF")), Label = sprintf("%.2f %s", `%IncMSE`,
    Sig))

summary_all <- bind_rows(rf_all$summary, rf_large$summary, rf_small$summary)

write.csv(importance_all, file = file.path(OUT_DIR, "RF_importance_All_Large_Small.csv"), row.names = FALSE)

write.csv(summary_all, file = file.path(OUT_DIR, "RF_model_summary.csv"), row.names = FALSE)

plot_data <- importance_all %>%
  group_by(Group) %>%
  arrange(`%IncMSE`, .by_group = TRUE) %>%
  mutate(Variable_facet = paste(Variable, Group, sep = "___")) %>%
  ungroup()

plot_data$Variable_facet <- factor(plot_data$Variable_facet, levels = plot_data$Variable_facet)

label_fun <- function(x) gsub("___.*$", "", x)

max_abs <- max(abs(plot_data$`%IncMSE`), na.rm = TRUE)

plot_data <- plot_data %>%
  mutate(Label_y = ifelse(`%IncMSE` >= 0, `%IncMSE` + 0.03 * max_abs, `%IncMSE` - 0.03 * max_abs),
    Label_hjust = ifelse(`%IncMSE` >= 0, 0, 1))

p_all <- ggplot(plot_data, aes(x = Variable_facet, y = `%IncMSE`, fill = `%IncMSE`)) + geom_hline(yintercept = 0,
  linewidth = 0.4, color = "black") + geom_col(width = 0.78, color = "black", linewidth = 0.25) + geom_text(aes(y = Label_y,
  label = Label, hjust = Label_hjust), size = 3.2, color = "black") + scale_y_continuous(expand = expansion(mult = c(0.18,
  0.16))) + coord_flip(clip = "off") + facet_wrap(~Group, ncol = 1, scales = "free_y") + scale_x_discrete(labels = label_fun) +
  scale_fill_gradientn(colors = c(brewer.pal(5, "Blues")[2], brewer.pal(5, "YlGnBu")[1], brewer.pal(6,
    "OrRd")[4]), values = rescale(c(min(plot_data$`%IncMSE`, na.rm = TRUE), median(plot_data$`%IncMSE`,
    na.rm = TRUE), max(plot_data$`%IncMSE`, na.rm = TRUE))), guide = "none") + labs(x = NULL, y = "Increase in MSE (%)",
  subtitle = "Stars: BH-adjusted permutation P (*<0.05, **<0.01, ***<0.001)") + theme_test(base_size = 12) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",
    linewidth = 0.35), strip.background = element_rect(fill = "white", color = "black", linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 12), axis.text.y = element_text(color = "black",
      size = 10), axis.text.x = element_text(color = "black", size = 10), axis.title.x = element_text(color = "black",
      size = 12), plot.margin = ggplot2::margin(8, 35, 8, 35))

print(p_all)

ggsave(file.path(OUT_DIR, "RF_importance_All_Large_Small.png"), p_all, width = 8, height = 12, dpi = 300)

ggsave(file.path(OUT_DIR, "RF_importance_All_Large_Small.pdf"), p_all, width = 8, height = 12)

plot_one_group <- function(gname) {
  pd <- importance_all %>%
    filter(Group == gname) %>%
    arrange(`%IncMSE`) %>%
    mutate(Variable = factor(Variable, levels = Variable))
  max_abs_one <- max(abs(pd$`%IncMSE`), na.rm = TRUE)
  pd <- pd %>%
    mutate(Label_y = ifelse(`%IncMSE` >= 0, `%IncMSE` + 0.03 * max_abs_one, `%IncMSE` - 0.03 * max_abs_one),
      Label_hjust = ifelse(`%IncMSE` >= 0, 0, 1))
  p <- ggplot(pd, aes(x = Variable, y = `%IncMSE`, fill = `%IncMSE`)) + geom_hline(yintercept = 0,
    linewidth = 0.4, color = "black") + geom_col(width = 0.78, color = "black", linewidth = 0.25) +
    geom_text(aes(y = Label_y, label = Label, hjust = Label_hjust), size = 3.4, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.16))) + coord_flip(clip = "off") + scale_fill_gradientn(colors = c(brewer.pal(5,
    "Blues")[2], brewer.pal(5, "YlGnBu")[1], brewer.pal(6, "OrRd")[4]), values = rescale(c(min(pd$`%IncMSE`,
    na.rm = TRUE), median(pd$`%IncMSE`, na.rm = TRUE), max(pd$`%IncMSE`, na.rm = TRUE))), guide = "none") +
    labs(title = gname, x = NULL, y = "Increase in MSE (%)", subtitle = "Stars: BH-adjusted permutation P (*<0.05, **<0.01, ***<0.001)") +
    theme_test(base_size = 12) + theme(panel.grid = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35), plot.title = element_text(hjust = 0.5,
      face = "bold", size = 14), axis.text.y = element_text(color = "black", size = 10), axis.text.x = element_text(color = "black",
      size = 10), axis.title.x = element_text(color = "black", size = 12), plot.margin = ggplot2::margin(8,
      18, 8, 8))
  ggsave(file.path(OUT_DIR, paste0("RF_importance11_", gname, ".png")), p, width = 7, height = 6, dpi = 300)
  return(p)
}

p_All <- plot_one_group("All")

p_Large <- plot_one_group("ConF")

p_Small <- plot_one_group("USF")

message("分析完成！结果已输出到：", OUT_DIR)

plot_one_group <- function(gname) {
  pd <- importance_all %>%
    dplyr::filter(Group == gname) %>%
    dplyr::arrange(`%IncMSE`) %>%
    dplyr::mutate(Variable = factor(Variable, levels = Variable))
  if (nrow(pd) == 0) {
    warning("没有找到分组：", gname)
    return(NULL)
  }
  max_abs_one <- max(abs(pd$`%IncMSE`), na.rm = TRUE)
  if (!is.finite(max_abs_one) || max_abs_one == 0)
    max_abs_one <- 1
  pd <- pd %>%
    dplyr::mutate(Label_y = ifelse(`%IncMSE` >= 0, `%IncMSE` + 0.01 * max_abs_one, `%IncMSE` - 0.01 *
      max_abs_one), Label_hjust = ifelse(`%IncMSE` >= 0, 0, 1), Group = as.character(Group))
  p <- ggplot(pd, aes(x = Variable, y = `%IncMSE`, fill = `%IncMSE`)) + geom_hline(yintercept = 0,
    linewidth = 0.4, color = "black") + geom_col(width = 0.78, color = "black", linewidth = 0.25) +
    geom_text(aes(y = Label_y, label = Label, hjust = Label_hjust), size = 3.6, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.16))) + coord_flip(clip = "off") + scale_fill_gradientn(colors = c(brewer.pal(5,
    "Blues")[2], brewer.pal(5, "YlGnBu")[1], brewer.pal(6, "OrRd")[4]), values = scales::rescale(c(min(pd$`%IncMSE`,
    na.rm = TRUE), median(pd$`%IncMSE`, na.rm = TRUE), max(pd$`%IncMSE`, na.rm = TRUE))), guide = "none") +
    labs(title = gname, x = NULL, y = "Increase in MSE (%)", subtitle = "Stars: BH-adjusted permutation P (*<0.05, **<0.01, ***<0.001)") +
    theme_test(base_size = 12) + theme(panel.grid = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35), plot.title = element_text(hjust = 0.5,
      face = "bold", size = 15), axis.text.y = element_text(color = "black", size = 10), axis.text.x = element_text(color = "black",
      size = 10), axis.title.x = element_text(color = "black", size = 12), plot.margin = ggplot2::margin(8,
      35, 8, 8))
  print(p)
  ggsave(filename = file.path(OUT_DIR, paste0("RF_importance_", gname, ".png")), plot = p, width = 10,
    height = 6, dpi = 300)
  ggsave(filename = file.path(OUT_DIR, paste0("RF_importance11_", gname, ".pdf")), plot = p, width = 10,
    height = 6)
  write.csv(pd, file = file.path(OUT_DIR, paste0("RF_importance11_", gname, ".csv")), row.names = FALSE)
  return(p)
}

p_All <- plot_one_group("All")

p_Large <- plot_one_group("ConF")

p_Small <- plot_one_group("USF")

message("分析完成！已分别输出 All、Large、Small 三张随机森林图。结果目录：",
  OUT_DIR)
