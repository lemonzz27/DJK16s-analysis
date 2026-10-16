library(vegan)

library(ggplot2)

set.seed(20260809)

otu_raw <- read.csv("tpm.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

otu <- data.frame(t(otu_raw), check.names = FALSE)

month_info <- read.csv("groupmonth.csv", stringsAsFactors = FALSE, check.names = FALSE)

size_info <- read.csv("groupsize.csv", stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c("sample", "group")

if (!all(required_columns %in% names(month_info))) {
  stop("groupmonth.csv 必须包含 sample 和 group 两列。")
}

if (!all(required_columns %in% names(size_info))) {
  stop("groupsize.csv 必须包含 sample 和 group 两列。")
}

if (anyDuplicated(month_info$sample)) {
  stop("groupmonth.csv 中存在重复的 sample。")
}

if (anyDuplicated(size_info$sample)) {
  stop("groupsize.csv 中存在重复的 sample。")
}

sample_order <- rownames(otu)

metadata <- data.frame(sample = sample_order, month = month_info$group[match(sample_order, month_info$sample)],
  size = size_info$group[match(sample_order, size_info$sample)], stringsAsFactors = FALSE)

metadata$month <- trimws(as.character(metadata$month))

metadata$size <- tolower(trimws(as.character(metadata$size)))

valid_sample <- complete.cases(metadata[, c("month", "size")]) & metadata$size %in% c("small", "large")

if (any(!valid_sample)) {
  warning(sum(!valid_sample), " 个样本缺少有效的月份或粒级信息，已从分析中剔除。")
}

metadata <- metadata[valid_sample, , drop = FALSE]

otu <- otu[metadata$sample, , drop = FALSE]

if (nrow(metadata) < 4) {
  stop("有效样本数不足，无法进行 PCoA 和 PERMANOVA。")
}

month_levels <- sort(unique(metadata$month))

if (length(month_levels) != 15) {
  stop("当前有效月份数为 ", length(month_levels), "，本配色按 15 个月设计，请检查 groupmonth.csv。")
}

metadata$month <- factor(metadata$month, levels = month_levels)

metadata$size <- factor(metadata$size, levels = c("small", "large"))

month_colors <- c("#7F0000", "#B30000", "#D7301F", "#EF6548", "#FC8D59", "#00441B", "#006D2C", "#238B45",
  "#41AB5D", "#74C476", "#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6")

names(month_colors) <- month_levels

size_shapes <- c(small = 21, large = 24)

size_labels <- c(small = "USF", large = "ConF")

cat("正在计算 Bray-Curtis 距离和 PCoA 坐标……\n")

bray_dis <- vegdist(otu, method = "bray")

pcoa_result <- cmdscale(bray_dis, k = 2, eig = TRUE)

explained_variance <- pcoa_result$eig/sum(pcoa_result$eig)

plot_data <- data.frame(sample = rownames(pcoa_result$points), PCoA1 = pcoa_result$points[, 1], PCoA2 = pcoa_result$points[,
  2], stringsAsFactors = FALSE)

plot_data$month <- metadata$month[match(plot_data$sample, metadata$sample)]

plot_data$size <- metadata$size[match(plot_data$sample, metadata$sample)]

set.seed(20260809)

permanova_month <- adonis2(otu ~ month, data = metadata, method = "bray", permutations = 999)

set.seed(20260809)

permanova_size <- adonis2(otu ~ size, data = metadata, method = "bray", permutations = 999)

extract_permanova <- function(result, factor_name) {
  data.frame(factor = factor_name, Df = result$Df[1], SumOfSqs = result$SumOfSqs[1], R2 = result$R2[1],
    F = result$F[1], P = result$`Pr(>F)`[1], permutations = 999, distance = "Bray-Curtis", model = "single-factor PERMANOVA",
    stringsAsFactors = FALSE)
}

permanova_table <- rbind(extract_permanova(permanova_month, "month"), extract_permanova(permanova_size,
  "size"))

write.csv(permanova_table, "PCoA_month_size_PERMANOVA.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("\nPERMANOVA 结果：\n")

print(permanova_table, row.names = FALSE)

format_p <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return("< 0.001")
  }
  sprintf("= %.3f", p_value)
}

month_r2 <- permanova_table$R2[permanova_table$factor == "month"]

month_p <- permanova_table$P[permanova_table$factor == "month"]

size_r2 <- permanova_table$R2[permanova_table$factor == "size"]

size_p <- permanova_table$P[permanova_table$factor == "size"]

stat_label <- paste0("PERMANOVA (Bray-Curtis; 999 permutations)\n", "Month: R² = ", sprintf("%.3f",
  month_r2), ", P ", format_p(month_p), "\n", "Size: R² = ", sprintf("%.3f", size_r2), ", P ", format_p(size_p))

small_data <- plot_data[plot_data$size == "small", , drop = FALSE]

large_data <- plot_data[plot_data$size == "large", , drop = FALSE]

p_combined <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) + stat_ellipse(data = small_data, geom = "polygon",
  type = "t", level = 0.95, fill = "#009E73", colour = NA, alpha = 0.1, show.legend = FALSE) + stat_ellipse(data = small_data,
  geom = "path", type = "t", level = 0.95, colour = "#009E73", linewidth = 0.75, linetype = 2, show.legend = FALSE) +
  stat_ellipse(data = large_data, geom = "polygon", type = "t", level = 0.95, fill = "#D55E00", colour = NA,
    alpha = 0.1, show.legend = FALSE) + stat_ellipse(data = large_data, geom = "path", type = "t",
  level = 0.95, colour = "#D55E00", linewidth = 0.75, linetype = 2, show.legend = FALSE) + geom_point(aes(fill = month,
  shape = size), size = 4, colour = "#2B2B2B", stroke = 0.3, alpha = 0.9) + scale_fill_manual(values = month_colors,
  breaks = month_levels, name = "Month") + scale_shape_manual(values = size_shapes, breaks = names(size_shapes),
  labels = size_labels, name = "Size") + labs(x = paste0("PCoA1 (", round(100 * explained_variance[1],
  1), "%)"), y = paste0("PCoA2 (", round(100 * explained_variance[2], 1), "%)"), title = "PCoA: month and size",
  subtitle = "Color = month; shape and 95% ellipse = size") + annotate("label", x = Inf, y = Inf, label = stat_label,
  hjust = 1.04, vjust = 1.08, size = 4.1, lineheight = 1.05, linewidth = 0.25, fill = "white", colour = "black") +
  guides(fill = guide_legend(order = 1, ncol = 1, override.aes = list(shape = 21, size = 4, colour = "#2B2B2B",
    alpha = 1)), shape = guide_legend(order = 2, ncol = 1, override.aes = list(fill = "grey70", colour = "#2B2B2B",
    size = 4, alpha = 1))) + theme_bw(base_family = "sans") + theme(panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", linewidth = 1.1), axis.title = element_text(size = 16,
    face = "bold"), axis.text = element_text(size = 13, colour = "black"), plot.title = element_text(size = 17,
    face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5), legend.position = "right",
  legend.box = "vertical", legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 12),
  plot.margin = margin(8, 8, 8, 8))

print(p_combined)

ggsave("PCoA_month_size_combined.png", p_combined, width = 9.2, height = 7.2, dpi = 600, bg = "white")

ggsave("PCoA_month_size_combined.pdf", p_combined, width = 9.2, height = 7.2, device = cairo_pdf, bg = "white")

write.csv(plot_data, "PCoA_month_size_scores.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("\n已输出：\n", "1. PCoA_month_size_combined.png\n", "2. PCoA_month_size_combined.pdf\n", "3. PCoA_month_size_PERMANOVA.csv\n",
  "4. PCoA_month_size_scores.csv\n")
