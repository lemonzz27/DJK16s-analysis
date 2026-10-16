local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA/PCOA1年际变化"))
  library(ggplot2)

  library(dplyr)

  library(vegan)

  source(normalizePath(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/gamm_functions.R"),
    winslash = "/", mustWork = TRUE), encoding = "UTF-8")

  tpm_data <- read.csv("tpm.csv", row.names = 1, stringsAsFactors = FALSE)

  group_month <- read.csv("groupmonth.csv", stringsAsFactors = FALSE)

  group_size <- read.csv("groupsize.csv", stringsAsFactors = FALSE)

  group_location <- read.csv("grouplocation.csv", stringsAsFactors = FALSE)

  colnames(group_month) <- c("sample", "month_group")

  colnames(group_size) <- c("sample", "size_group")

  colnames(group_location) <- c("sample", "site_group")

  tpm_t <- t(tpm_data)

  dist_matrix <- vegdist(tpm_t, method = "bray")

  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

  pcoa1_var <- round(pcoa_res$eig[1]/sum(pcoa_res$eig) * 100, 2)

  pcoa_df <- data.frame(sample = rownames(tpm_t), PCoA1 = pcoa_res$points[, 1], stringsAsFactors = FALSE)

  write.csv(pcoa_df, "1.csv")

  merged_data <- pcoa_df %>%
    inner_join(group_month, by = "sample") %>%
    inner_join(group_size, by = "sample") %>%
    inner_join(group_location, by = "sample")

  merged_data$Date <- as.Date(paste0("20", merged_data$month_group, "01"), format = "%Y%m%d")

  merged_data$cluster <- interaction(merged_data$site_group, merged_data$size_group, drop = TRUE)

  data_all <- merged_data

  data_all$size_group <- "All"

  plot_data <- bind_rows(merged_data, data_all)

  plot_data$size_group <- factor(plot_data$size_group, levels = c("large", "small", "All"))

  my_colors <- c(large = "#D95F02", small = "#1B9E77", All = "#4DBBD5")

  pcoa1_gamm <- fit_repeated_time_gamm(data = merged_data, response = "PCoA1", cluster_col = "cluster",
    k = 6L, include_overall = TRUE, use_random_intercept = TRUE)

  gamm_statistics <- add_gamm_bh_adjustment(pcoa1_gamm$statistics, include_overall = TRUE)

  gamm_predictions <- pcoa1_gamm$predictions

  write.csv(gamm_statistics, "GAMM_PCoA1_Model_Statistics.csv", row.names = FALSE)

  write.csv(gamm_predictions, "GAMM_PCoA1_Predictions.csv", row.names = FALSE)

  p_pcoa1 <- ggplot(plot_data, aes(x = Date, y = PCoA1, color = size_group, fill = size_group)) + geom_ribbon(data = gamm_predictions,
    aes(x = Date, ymin = lower, ymax = upper, fill = size_group), inherit.aes = FALSE, alpha = 0.2, color = NA) +
    geom_line(data = gamm_predictions, aes(x = Date, y = fit, color = size_group), inherit.aes = FALSE,
      linewidth = 0.9) + geom_point(alpha = 0.4, position = position_jitter(width = 5)) + scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) + scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") +
    theme_classic() + labs(title = "Trend of Beta Diversity (PCoA1) over Time", subtitle = gamm_subtitle(gamm_statistics,
    include_overall = TRUE), x = "Month", y = paste0("PCoA1 (", pcoa1_var, "%)"), color = "Group", fill = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8.5,
      lineheight = 1.05))

  ggsave("PCoA1_Trend.png", plot = p_pcoa1, width = 8, height = 6, dpi = 300)

  if (interactive()) print(p_pcoa1)
})
