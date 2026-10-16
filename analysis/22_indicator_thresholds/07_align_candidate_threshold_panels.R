Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")

.libPaths(c(file.path(getwd(), "Rlib"), .libPaths()))

library(ggplot2)

library(dplyr)

library(grid)

src <- readLines(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/21_candidate_dynamics/03_plot_candidate_wide_layout.R"), encoding = "UTF-8")

eval(parse(text = src[1:(grep("^combined_plot <-", src)[1] - 1)]))

top <- p_main + annotation_custom(ggplotGrob(p_water_inset), xmin = as.numeric(as.Date("2025-01-01")),
  xmax = as.numeric(as.Date("2025-07-25")), ymin = 143, ymax = 246) + annotate("text", x = as.Date("2025-08-01"),
  y = 244, label = "Overall GAMM: richness P<0.001; relative abundance P<0.001 (logit scale)", hjust = 0,
  size = 3.3, colour = "#444444") + theme(legend.position = "inside", legend.position.inside = c(0.99,
  0.918), legend.justification = c(1, 1), legend.direction = "horizontal", legend.text = element_text(size = 10),
  legend.key.width = unit(9, "mm"), legend.background = element_rect(fill = "white", colour = NA),
  legend.margin = margin(0, 0, 0, 0))

root <- "01_分析项目/36_TITAN2致病菌ASV水平_水位"

cache <- file.path(root, "结果表格/WL/ASV/plotdata_figure1_community_sumz")

profile_curves <- read.csv(file.path(cache, "03_shape_preserving_curve_plotdata.csv"))

threshold_heights <- read.csv(file.path(cache, "05_plotted_robust_thresholds.csv"))

n_analyzed <- 476

n_robust <- 181

src2 <- readLines(file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/22_indicator_thresholds/02_plot_community_thresholds.R"),
  encoding = "UTF-8")

eval(parse(text = src2[grep("^direction_colors <-", src2)[1]:(grep("^output_prefix <-", src2)[1] - 1)]))

bottom <- sumz_plot + labs(title = NULL, subtitle = NULL) + scale_y_continuous(expand = expansion(mult = c(0.06,
  0.18))) + guides(colour = guide_legend(order = 2, nrow = 1, override.aes = list(linetype = "solid",
  alpha = 1, linewidth = 1.2)), linetype = guide_legend(order = 1, nrow = 1, override.aes = list(colour = "#4A4A4A",
  alpha = c(0.32, 1), linewidth = c(0.9, 1.15)))) + theme(legend.position = "inside", legend.position.inside = c(0.005,
  0.995), legend.justification = c(0, 1), legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0,
  "mm"), legend.margin = margin(2, 3, 2, 3), legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = NA),
  axis.text = element_text(size = 11), axis.title = element_text(size = 12, face = "bold"), panel.border = element_rect(fill = NA,
    colour = "black", linewidth = 0.8), plot.margin = margin(5, 8, 5, 8))

g1 <- ggplotGrob(top)

g2 <- ggplotGrob(bottom)

stopifnot(length(g1$widths) == length(g2$widths), length(g1$heights) == length(g2$heights))

g1$widths <- g2$widths <- unit.pmax(g1$widths, g2$widths)

g1$heights <- g2$heights <- unit.pmax(g1$heights, g2$heights)

xi <- which(g2$layout$name == "xlab-b")

g2$grobs[[xi]]$children[[1]]$y <- g2$grobs[[xi]]$children[[1]]$y + unit(22, "pt")

out <- "Pathogen_TITAN_aligned_20260921"

dir.create(out, showWarnings = FALSE)

save_plot <- function(name, grobs, preview = FALSE) {
  w <- 11.8
  h <- 6.3 * length(grobs)
  for (ext in if (preview)
    "png"
  else c("png", "pdf", "svg")) {
    file <- file.path(out, paste0(name, if (preview)
      "_preview"
    else "", ".", ext))
    if (ext == "png")
      ragg::agg_png(file, width = w, height = h, units = "in", res = if (preview)
        110
      else 600, background = "white")
    if (ext == "pdf")
      cairo_pdf(file, width = w, height = h, bg = "white")
    if (ext == "svg")
      svglite::svglite(file, width = w, height = h, bg = "white")
    grid.newpage()
    for (i in seq_along(grobs)) {
      pushViewport(viewport(x = 0.5, y = 1 - (i - 0.5)/length(grobs), width = 1, height = 1/length(grobs)))
      grid.draw(grobs[[i]])
      popViewport()
    }
    dev.off()
  }
}

save_plot("Pathogen_temporal_aligned", list(g1))

save_plot("TITAN_sumz_aligned", list(g2))

save_plot("Pathogen_TITAN_pair_aligned", list(g1, g2))

save_plot("Pathogen_TITAN_pair_aligned", list(g1, g2), TRUE)
