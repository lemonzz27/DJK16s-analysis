local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA/a多样性"))
  library(vegan)

  alpha <- function(x, tree = NULL, base = exp(1)) {
    est <- estimateR(x)
    Richness <- est[1, ]
    Chao1 <- est[2, ]
    ACE <- est[4, ]
    Shannon <- diversity(x, index = "shannon", base = base)
    Simpson <- diversity(x, index = "simpson")
    Pielou <- Shannon/log(Richness, base)
    goods_coverage <- 1 - rowSums(x == 1)/rowSums(x)
    result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
    if (!is.null(tree)) {
      PD_whole_tree <- picante::pd(x, tree, include.root = FALSE)[1]
      names(PD_whole_tree) <- "PD_whole_tree"
      result <- cbind(result, PD_whole_tree)
    }
    return(result)
  }

  otu <- read.csv("tpm.csv", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

  otu <- t(otu)

  otu <- as.data.frame(otu)

  sample_depths <- rowSums(otu)

  min_depth <- min(sample_depths)

  set.seed(123)

  otu_rarefied <- rrarefy(otu, sample = min_depth)

  otu_rarefied <- as.data.frame(otu_rarefied)

  alpha_all <- alpha(otu_rarefied, base = 2)

  write.csv(alpha_all, "alpha1.csv", quote = FALSE)
})

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA/a多样性"))
  rm(list = ls())

  library(tidyverse)

  library(ggpubr)

  library(gghalves)

  library(multcomp)

  library(doBy)

  analyze_and_plot_alpha_cou <- function(alpha_file, group_file, metrics, output_prefix, group_levels,
    color_map, plot_width = 5.2, plot_height = 4.5) {
    alpha <- read.csv(alpha_file, row.name = 1)
    group_info <- read.csv(group_file)
    samples <- rownames(alpha)
    group_vec <- group_info$group[match(samples, group_info$sample)]
    alpha$group <- factor(group_vec, levels = group_levels)
    alpha <- alpha[!is.na(alpha$group), ]
    display_labels <- if (identical(group_levels, c("large", "small"))) {
      c(large = "ConF", small = "USF")
    }
    else {
      stats::setNames(group_levels, group_levels)
    }
    results <- list()
    for (metric in metrics) {
      if (!metric %in% colnames(alpha)) {
        warning(paste("指标", metric, "不存在于数据中，跳过"))
        next
      }
      data <- data.frame(Value = alpha[[metric]], group = alpha$group)
      fit <- aov(Value ~ group, data = data)
      tuk <- glht(fit, linfct = mcp(group = "Tukey"))
      Tukey_res <- cld(tuk, level = 0.05, decreasing = TRUE)
      Tukey_letters <- Tukey_res$mcletters$Letters
      y_range <- diff(range(data$Value))
      stat <- data %>%
        group_by(group) %>%
        summarise(Position = max(Value) + (y_range * 0.12), .groups = "drop") %>%
        mutate(sig = Tukey_letters[as.character(group)])
      plot <- ggplot(data, aes(x = group, y = Value, fill = group)) + geom_half_violin(side = "l",
        trim = TRUE, alpha = 0.35, color = NA, width = 0.8) + geom_boxplot(width = 0.12, outlier.shape = NA,
        alpha = 0.8, size = 0.5, color = "black") + geom_jitter(shape = 21, size = 2.2, width = 0.1,
        alpha = 0.8, stroke = 0.3, color = "white") + geom_text(data = stat, aes(x = group, y = Position,
        label = sig), size = 6, fontface = "bold", color = "black", inherit.aes = FALSE) + scale_fill_manual(values = color_map,
        breaks = group_levels, labels = display_labels) + scale_x_discrete(labels = display_labels) +
        labs(x = NULL, y = metric) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_rect(linewidth = 1.2,
        color = "black"), axis.text.x = element_text(color = "black", size = 12, face = "bold", angle = if (length(group_levels) >
        6)
        45
      else 0, hjust = if (length(group_levels) > 6)
        1
      else 0.5), axis.text.y = element_text(color = "black", size = 11), axis.title.y = element_text(size = 15,
        face = "bold"), axis.ticks = element_line(color = "black", linewidth = 0.6), legend.position = "none")
      ggsave(filename = paste0(output_prefix, "_", metric, ".png"), plot = plot, width = plot_width,
        height = plot_height, dpi = 600)
      results[[metric]] <- list(data = data, plot = plot)
    }
    return(results)
  }

  metrics_to_analyze <- c("Richness", "Shannon", "Simpson", "Chao1")

  levels_1 <- c("large", "small")

  levels_2 <- c("spring", "summer", "autumn", "winter")

  levels_3 <- c("dry", "wet")

  levels_4 <- c("BS", "KX", "LHKX", "TC")

  colors_1 <- c(large = "#D95F02", small = "#1B9E77")

  colors_2 <- c(spring = "#95D5B2", summer = "#FFBA08", autumn = "#F4A261", winter = "#A8DADC")

  colors_3 <- c(dry = "#DCA04E", wet = "#2E74B5")

  colors_4 <- c(BS = "#95D5B2", KX = "#FFBA08", LHKX = "#F4A261", TC = "#A8DADC")

  levels_5 <- sort(unique(read.csv("groupmonth.csv", stringsAsFactors = FALSE)$group))

  colors_5 <- setNames(grDevices::hcl.colors(length(levels_5), palette = "Dark 3"), levels_5)

  results_size_all <- analyze_and_plot_alpha_cou(alpha_file = "alpha1.csv", group_file = "groupsize.csv",
    metrics = metrics_to_analyze, output_prefix = "alpha_size", group_levels = levels_1, color_map = colors_1,
    plot_width = 4.2, plot_height = 5.4)

  results_season_all <- analyze_and_plot_alpha_cou(alpha_file = "alpha1.csv", group_file = "groupseason.csv",
    metrics = metrics_to_analyze, output_prefix = "alpha_season", group_levels = levels_2, color_map = colors_2)

  results_rain_all <- analyze_and_plot_alpha_cou(alpha_file = "alpha1.csv", group_file = "grouprain.csv",
    metrics = metrics_to_analyze, output_prefix = "alpha_rain", group_levels = levels_3, color_map = colors_3)

  results_month_all <- analyze_and_plot_alpha_cou(alpha_file = "alpha1.csv", group_file = "groupmonth.csv",
    metrics = metrics_to_analyze, output_prefix = "alpha_month", group_levels = levels_5, color_map = colors_5)

  results_area_all <- analyze_and_plot_alpha_cou(alpha_file = "alpha1.csv", group_file = "grouplocation.csv",
    metrics = metrics_to_analyze, output_prefix = "alpha_location", group_levels = levels_4, color_map = colors_4)
})
