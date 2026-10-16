required_packages <- c("ggplot2", "dplyr", "lubridate")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}

library(ggplot2)

library(dplyr)

library(lubridate)

source(normalizePath(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/07_alpha_diversity/00_gamm_functions.R"),
  winslash = "/", mustWork = TRUE), encoding = "UTF-8")

alpha_data <- read.csv("alpha1.csv", stringsAsFactors = FALSE)

colnames(alpha_data)[1] <- "sample"

group_month <- read.csv("groupmonth.csv", stringsAsFactors = FALSE)

group_size <- read.csv("groupsize.csv", stringsAsFactors = FALSE)

group_location <- read.csv("grouplocation.csv", stringsAsFactors = FALSE)

colnames(group_month)[2] <- "month_group"

colnames(group_size)[2] <- "size_group"

colnames(group_location)[2] <- "site_group"

merged_data <- alpha_data %>%
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

alpha_metrics <- c("Richness", "Shannon", "Chao1")

gamm_fits <- lapply(alpha_metrics, function(metric_name) {
  fit_repeated_time_gamm(data = merged_data, response = metric_name, cluster_col = "cluster", k = 6L,
    include_overall = TRUE, use_random_intercept = TRUE)
})

names(gamm_fits) <- alpha_metrics

gamm_statistics <- bind_rows(lapply(gamm_fits, `[[`, "statistics"))

gamm_statistics <- add_gamm_bh_adjustment(gamm_statistics, include_overall = TRUE)

gamm_predictions <- bind_rows(lapply(gamm_fits, `[[`, "predictions"))

write.csv(gamm_statistics, "GAMM_Alpha_Model_Statistics.csv", row.names = FALSE)

write.csv(gamm_predictions, "GAMM_Alpha_Predictions.csv", row.names = FALSE)

plot_alpha_gamm <- function(metric_name, plot_title, y_label, file_name) {
  curve_data <- gamm_predictions %>%
    filter(metric == metric_name)
  statistics_row <- gamm_statistics %>%
    filter(metric == metric_name)
  p <- ggplot(plot_data, aes(x = Date, y = .data[[metric_name]], color = size_group, fill = size_group)) +
    geom_ribbon(data = curve_data, aes(x = Date, ymin = lower, ymax = upper, fill = size_group),
      inherit.aes = FALSE, alpha = 0.2, color = NA) + geom_line(data = curve_data, aes(x = Date,
    y = fit, color = size_group), inherit.aes = FALSE, linewidth = 0.9) + geom_point(alpha = 0.4,
    position = position_jitter(width = 5)) + scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors) +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") + theme_classic() + labs(title = plot_title,
    subtitle = gamm_subtitle(statistics_row, include_overall = TRUE), x = "Month", y = y_label, color = "Group",
    fill = "Group") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8.5,
    lineheight = 1.05))
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
  p
}

p_richness <- plot_alpha_gamm("Richness", "Trend of Alpha Diversity (Richness) over Time", "Richness",
  "Richness_Trend.png")

p_shannon <- plot_alpha_gamm("Shannon", "Trend of Alpha Diversity (Shannon Index) over Time", "Shannon Index",
  "Shannon_Trend.png")

p_chao1 <- plot_alpha_gamm("Chao1", "Trend of Alpha Diversity (Chao1 Index) over Time", "Chao1 Index",
  "chao1_Trend.png")

if (interactive()) {
  print(p_richness)
  print(p_shannon)
  print(p_chao1)
}
