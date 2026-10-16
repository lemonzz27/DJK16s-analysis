palette <- c(Pseudomonadota = "#DF89FF", Actinomycetota = "#73C000", Bacteroidota = "#FF8805", Cyanobacteriota = "#4C463E",
  unclassified = "#00C4FF", Acidobacteriota = "#D3B3B0", Chloroflexota = "#FF5584", Bdellovibrionota = "#A6761D",
  Thermoproteota = "#00BD94", Planctomycetota = "#6699CC", Deinococcota = "#E6AB02", `Other phyla` = "#C0C0C0")

draw_legend <- function() {
  par(mar = c(0.1, 0.1, 0.1, 0.1), family = "Arial", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 10.7), ylim = c(0, 5.2), asp = 1)
  text(0.28, 4.92, "Network legend", adj = c(0, 0.5), cex = 1.25, font = 2, col = "#202124")
  text(0.28, 4.3, "Node color (Phylum)", adj = c(0, 0.5), cex = 0.98, font = 2, col = "#202124")
  labels <- names(palette)
  first <- labels[1:6]
  second <- labels[7:12]
  rows <- seq(3.72, 0.72, length.out = 6)
  for (i in seq_along(first)) {
    points(0.55, rows[i], pch = 21, cex = 1.55, bg = palette[first[i]], col = "#4A4A4A", lwd = 0.65)
    text(0.88, rows[i], first[i], adj = c(0, 0.5), cex = 0.92, col = "#202124")
  }
  for (i in seq_along(second)) {
    points(3.1, rows[i], pch = 21, cex = 1.55, bg = palette[second[i]], col = "#4A4A4A", lwd = 0.65)
    text(3.43, rows[i], second[i], adj = c(0, 0.5), cex = 0.92, col = "#202124")
  }
  segments(5.82, 0.48, 5.82, 4.5, col = "#D9D9D9", lwd = 0.9)
  text(6.15, 4.3, "Node size (Mean)", adj = c(0, 0.5), cex = 0.98, font = 2, col = "#202124")
  size_values <- c(10, 20, 30, 40)
  size_x <- c(6.5, 7.35, 8.3, 9.35)
  size_cex <- c(0.8, 1.25, 1.7, 2.15)
  for (i in seq_along(size_values)) {
    points(size_x[i], 3.48, pch = 21, cex = size_cex[i], bg = "#C0C0C0", col = "#4A4A4A", lwd = 0.65)
    text(size_x[i], 2.88, as.character(size_values[i]), cex = 0.86, col = "#202124")
  }
  text(6.15, 2.25, "Edge color (Association)", adj = c(0, 0.5), cex = 0.98, font = 2, col = "#202124")
  segments(6.32, 1.58, 7.9, 1.58, col = "#EE7993", lwd = 4.1, lend = 1)
  text(8.2, 1.58, "Positive", adj = c(0, 0.5), cex = 0.92, col = "#202124")
  segments(6.32, 0.78, 7.9, 0.78, col = "#3FB17D", lwd = 4.1, lend = 1)
  text(8.2, 0.78, "Negative", adj = c(0, 0.5), cex = 0.92, col = "#202124")
}

output_dir <- paste0(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/11_静态共现网络/"),
  "结果图片/SILVA/整体_静态网络_出现率20pct_增强版")

png(file.path(output_dir, "Static_Network_Legend.png"), width = 6.65, height = 3.2, units = "in", res = 600,
  type = "cairo", bg = "white", pointsize = 10)

draw_legend()

dev.off()

svg(file.path(output_dir, "Static_Network_Legend.svg"), width = 6.65, height = 3.2, bg = "white", pointsize = 10,
  family = "Arial")

draw_legend()

dev.off()

cairo_pdf(file.path(output_dir, "Static_Network_Legend.pdf"), width = 6.65, height = 3.2, bg = "white",
  family = "Arial", pointsize = 10)

draw_legend()

dev.off()
