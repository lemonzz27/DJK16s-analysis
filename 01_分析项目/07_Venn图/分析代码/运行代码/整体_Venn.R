library(tidyverse)

library(VennDiagram)

library(grid)

tpm_data <- read.csv(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "00_输入文件/01_ASV与序列/ASV原始reads_36097ASV_120样本.csv"),
  row.names = 1, check.names = FALSE)

group_data <- read.csv(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "00_输入文件/02_样本分组与站点/孔径分组_groupsize.csv"),
  check.names = FALSE)

colnames(tpm_data) <- gsub("\\.", "_", colnames(tpm_data))

group_data$sample <- gsub("\\.", "_", group_data$sample)

tpm_long <- tpm_data %>%
  rownames_to_column(var = "ASV") %>%
  pivot_longer(cols = -ASV, names_to = "sample", values_to = "abundance")

plot_data_summary <- tpm_long %>%
  left_join(group_data, by = "sample") %>%
  filter(abundance > 0)

asvs_small <- plot_data_summary %>%
  filter(group == "small") %>%
  pull(ASV) %>%
  unique()

asvs_large <- plot_data_summary %>%
  filter(group == "large") %>%
  pull(ASV) %>%
  unique()

n_small <- length(asvs_small)

n_large <- length(asvs_large)

n_overlap <- length(intersect(asvs_small, asvs_large))

cat("Small =", n_small, "\n")

cat("Large =", n_large, "\n")

cat("Overlap =", n_overlap, "\n")

need_invert <- n_large > n_small

cat("是否需要左右镜像：", need_invert, "\n")

venn_plot <- draw.pairwise.venn(area1 = n_small, area2 = n_large, cross.area = n_overlap, category = c("USF",
  "ConF"), fill = c("#1B9E77", "#D95F02"), alpha = c(0.75, 0.75), col = c("black", "black"), lwd = c(1.5,
  1.5), lty = c("solid", "solid"), euler.d = TRUE, scaled = TRUE, inverted = need_invert, rotation.degree = 0,
  cex = c(2.5, 2.3, 2.5), fontface = "bold", fontfamily = "sans", cat.cex = 1.95, cat.fontface = "bold",
  cat.fontfamily = "sans", cat.col = c("#1B9E77", "#D95F02"), cat.pos = c(-180, 0), cat.dist = c(0.055,
    0.055), print.mode = c("raw", "percent"), sigdigs = 2, ind = FALSE)

png(filename = "Figure_Venn_Small_Left_Large_Right.png", width = 2400, height = 2400, res = 300, type = "cairo-png",
  bg = "transparent")

grid.newpage()

grid.draw(venn_plot)

grid.text("Number of Detected ASVs", x = 0.5, y = 0.95, gp = gpar(fontsize = 22, fontface = "bold", fontfamily = "sans",
  col = "black"))

dev.off()

cairo_pdf(filename = "Figure_Venn_Small_Left_Large_Right.pdf", width = 8, height = 8, bg = "transparent")

grid.newpage()

grid.draw(venn_plot)

grid.text("Number of Detected ASVs", x = 0.5, y = 0.95, gp = gpar(fontsize = 22, fontface = "bold", fontfamily = "sans",
  col = "black"))

dev.off()

cat("\n最终位置：\n")

cat("左边 = Small =", n_small, "\n")

cat("右边 = Large =", n_large, "\n")

cat("共同 ASVs =", n_overlap, "\n")

cat("\n图片已保存：\n")

cat("Figure_Venn_Small_Left_Large_Right.png\n")

cat("Figure_Venn_Small_Left_Large_Right.pdf\n")
