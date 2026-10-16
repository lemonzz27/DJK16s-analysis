options(stringsAsFactors = FALSE)

input_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/icamp_assembly_plots")

output_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "main_figure_panels_large_font_20260902/Figure_4/iCAMP")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_cols <- c("Heterogeneous Selection", "Homogeneous Selection", "Dispersal Limitation", "Homogenizing Dispersal",
  "Drift and Others")

process_colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")

draw_icamp_month <- function() {
  dat <- read.csv(file.path(input_dir, "plot_data_month.csv"), check.names = FALSE)
  values <- as.matrix(dat[process_cols])
  values <- values/rowSums(values)
  labels <- dat$Month
  par(mar = c(6.6, 5.8, 1, 14), family = "sans", mgp = c(3.1, 0.85, 0))
  mids <- barplot(t(values), col = process_colors, border = "black", lwd = 0.45, width = 0.7, space = 0.45,
    axes = FALSE, names.arg = rep("", nrow(values)), ylim = c(0, 1))
  axis(1, at = mids, labels = FALSE, tick = FALSE)
  usr <- par("usr")
  text(mids, usr[3L] - 0.035 * (usr[4L] - usr[3L]), labels, srt = 60, adj = 1, xpd = NA, cex = 1)
  axis(2, at = seq(0, 1, by = 0.2), labels = paste0(seq(0, 100, by = 20), "%"), las = 1, tcl = -0.25,
    cex.axis = 1.05)
  box(bty = "l", lwd = 0.8)
  mtext("Relative importance", side = 2, line = 3.4, cex = 1.12)
  mtext("Month", side = 1, line = 5, cex = 1.12)
  legend(x = usr[2L] + 0.025 * (usr[2L] - usr[1L]), y = usr[4L], xjust = 0, yjust = 1, xpd = NA, legend = process_cols,
    fill = process_colors, border = "black", bty = "n", title = "Assembly mechanisms", cex = 1, y.intersp = 1.1)
}

save_base <- function(stem, width, height, draw_fun) {
  grDevices::png(file.path(output_dir, paste0(stem, ".png")), width = width, height = height, units = "in",
    res = 600, type = "cairo", bg = "white")
  draw_fun()
  grDevices::dev.off()
  grDevices::cairo_pdf(file.path(output_dir, paste0(stem, ".pdf")), width = width, height = height,
    bg = "white")
  draw_fun()
  grDevices::dev.off()
}

save_base("iCAMP_month_large_font", 12, 5.5, draw_icamp_month)

draw_ratio_trend <- function() {
  obs <- read.csv(file.path(input_dir, "stochastic_deterministic_ratio_pairwise_data.csv"), check.names = FALSE)
  pred <- read.csv(file.path(input_dir, "stochastic_deterministic_ratio_GAMM_predictions.csv"), check.names = FALSE)
  obs$Date <- as.Date(obs$Date)
  pred$Date <- as.Date(pred$Date)
  groups <- c("ConF", "USF", "All")
  colors <- c(ConF = "#D55E00", USF = "#009E73", All = "#0072B2")
  xlim <- range(c(obs$Date, pred$Date)) + c(-7, 7)
  ylim <- c(0, 2)
  par(mar = c(6.2, 6.1, 1, 1.2), family = "sans", mgp = c(3.2, 0.85, 0), xaxs = "i", yaxs = "i")
  plot(xlim, ylim, type = "n", axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim)
  for (g in groups) {
    z <- pred[pred$Group == g, ]
    z <- z[order(z$Date), ]
    polygon(c(z$Date, rev(z$Date)), c(z$lower, rev(z$upper)), col = grDevices::adjustcolor(colors[g],
      alpha.f = 0.18), border = NA)
  }
  abline(h = 1, col = "#666666", lty = 2, lwd = 1)
  set.seed(20260824)
  for (g in groups) {
    z <- obs[obs$Group == g & obs$Ratio >= ylim[1L] & obs$Ratio <= ylim[2L], ]
    points(z$Date + stats::runif(nrow(z), -5, 5), z$Ratio, pch = 16, cex = 0.58, col = grDevices::adjustcolor(colors[g],
      alpha.f = 0.36))
  }
  for (g in groups) {
    z <- pred[pred$Group == g, ]
    z <- z[order(z$Date), ]
    lines(z$Date, z$fit, col = colors[g], lwd = 2.6)
  }
  month_ticks_all <- sort(unique(obs$Date))
  month_ticks <- month_ticks_all[seq(1L, length(month_ticks_all), by = 2L)]
  axis(1, at = month_ticks, labels = FALSE, tcl = -0.22)
  usr <- par("usr")
  text(month_ticks, usr[3L] - 0.025 * diff(usr[3:4]), format(month_ticks, "%Y-%m"), srt = 45, adj = 1,
    xpd = NA, cex = 1.05)
  axis(2, at = seq(0, 2, by = 0.5), labels = sprintf("%.1f", seq(0, 2, by = 0.5)), las = 1, tcl = -0.22,
    cex.axis = 1.05)
  box(bty = "o", lwd = 0.85)
  mtext("Month", side = 1, line = 4.9, cex = 1.15, font = 2)
  mtext("Stochastic / deterministic ratio", side = 2, line = 3.8, cex = 1.15, font = 2)
  text(usr[2L] - 0.01 * diff(usr[1:2]), 1.05, "Ratio = 1", adj = 1, col = "#555555", cex = 0.9)
  legend("topright", inset = c(0.025, 0.025), legend = groups, col = colors[groups], lwd = 2.6, pch = 16,
    pt.cex = 0.8, bty = "n", title = "Group", cex = 1.02, y.intersp = 1.15)
}

save_base("iCAMP_stochastic_deterministic_ratio_trend_large_font", 12, 6, draw_ratio_trend)

message("Saved enlarged iCAMP panels to: ", output_dir)
