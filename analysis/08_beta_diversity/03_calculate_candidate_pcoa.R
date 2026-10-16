library(vegan)

library(ggplot2)

otu_raw <- read.csv("tpm.csv", row.names = 1, stringsAsFactors = FALSE)

otu <- data.frame(t(otu_raw))

colors_ele <- c(large = "#D95F02", small = "#1B9E77")

colors_ew <- c(spring = "#95D5B2", summer = "#FFBA08", autumn = "#F4A261", winter = "#A8DADC")

colors_cou <- c(dry = "#DCA04E", wet = "#2E74B5")

colors_peo <- c(BS = "#95D5B2", KX = "#FFBA08", LHKX = "#F4A261", TC = "#A8DADC")

month_levels <- sort(unique(read.csv("groupmonth.csv", stringsAsFactors = FALSE)$group))

colors_month <- setNames(grDevices::hcl.colors(length(month_levels), palette = "Dark 3"), month_levels)

cat("正在计算全局 Bray-Curtis 距离矩阵与 PCoA 坐标...\n")

bray_dis <- vegdist(otu, method = "bray")

pcoa <- cmdscale(bray_dis, k = 2, eig = TRUE)

exp_var <- pcoa$eig/sum(pcoa$eig)

pcoa_global <- data.frame(pcoa$point)[1:2]

colnames(pcoa_global) <- c("x", "y")

pcoa_global$sample <- rownames(pcoa_global)

plot_pcoa_cns <- function(global_coords, otu_data, exp_var, group_file, color_map, title_text) {
  display_labels <- if (setequal(names(color_map), c("large", "small"))) {
    c(large = "ConF", small = "USF")
  }
  else {
    stats::setNames(names(color_map), names(color_map))
  }
  group_info <- read.csv(group_file, stringsAsFactors = FALSE)
  group_info$group <- trimws(group_info$group)
  group_info$group <- factor(group_info$group, levels = names(color_map))
  common_samples <- intersect(global_coords$sample, group_info$sample)
  plot_df <- global_coords[match(common_samples, global_coords$sample), ]
  group_sub <- group_info[match(common_samples, group_info$sample), ]
  plot_df$group <- group_sub$group
  valid_idx <- !is.na(plot_df$group)
  plot_df <- plot_df[valid_idx, ]
  group_sub <- group_sub[valid_idx, ]
  valid_samples <- plot_df$sample
  otu_sub <- otu_data[valid_samples, ]
  adonis_res <- adonis2(otu_sub ~ group, data = group_sub, distance = "bray", permutations = 999)
  r2_val <- sprintf("%.3f", adonis_res$R2[1])
  p_val <- sprintf("%.3f", adonis_res$`Pr(>F)`[1])
  stat_label <- paste0("atop(bold('PERMANOVA'), italic(R)^2 == ", r2_val, " ~ italic(P) == ", p_val,
    ")")
  p <- ggplot(plot_df, aes(x = x, y = y, color = group, fill = group)) + stat_ellipse(geom = "polygon",
    level = 0.95, alpha = 0.1, linetype = 2, linewidth = 0.3, show.legend = FALSE) + geom_point(size = 4,
    alpha = 0.85, shape = 21, color = "white", stroke = 0.3) + scale_color_manual(values = color_map,
    labels = display_labels) + scale_fill_manual(values = color_map, labels = display_labels) + scale_y_continuous(expand = expansion(mult = c(0.1,
    0.2))) + labs(x = paste0("PCoA1 (", round(100 * exp_var[1], 1), "%)"), y = paste0("PCoA2 (",
    round(100 * exp_var[2], 1), "%)"), title = title_text) + theme_bw() + theme(panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1.2, color = "black"), axis.title = element_text(size = 15,
      face = "bold"), axis.text = element_text(size = 12, color = "black"), legend.title = element_blank(),
    legend.text = element_text(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    text = element_text(family = "sans")) + annotate("text", x = Inf, y = Inf, label = stat_label,
    parse = TRUE, size = 5, hjust = 1.1, vjust = 1.5, color = "black")
  return(p)
}

cat("正在生成图表...\n")

p_ele <- plot_pcoa_cns(pcoa_global, otu, exp_var, "groupsize.csv", colors_ele, "size")

p_ew <- plot_pcoa_cns(pcoa_global, otu, exp_var, "groupseason.csv", colors_ew, "season")

p_rain <- plot_pcoa_cns(pcoa_global, otu, exp_var, "grouprain.csv", colors_cou, "rain")

p_month <- plot_pcoa_cns(pcoa_global, otu, exp_var, "groupmonth.csv", colors_month, "month")

p_peo <- plot_pcoa_cns(pcoa_global, otu, exp_var, "grouplocation.csv", colors_peo, "location")

ggsave("PCoA_size.png", p_ele, width = 7.5, height = 6, dpi = 600)

ggsave("PCoA_season.png", p_ew, width = 7.5, height = 6, dpi = 600)

ggsave("PCoA_rain.png", p_rain, width = 7.5, height = 6, dpi = 600)

ggsave("PCoA_month.png", p_month, width = 7.5, height = 6, dpi = 600)

ggsave("PCoA_location.png", p_peo, width = 7.5, height = 6, dpi = 600)

cat("所有 5 张 PCoA 图片已成功生成并保存！\n")
