options(stringsAsFactors = FALSE)

input_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work_artifact/icamp_data")

output_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work_artifact/fig4a_updated")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_cols <- c("Heterogeneous Selection", "Homogeneous Selection", "Dispersal Limitation", "Homogenizing Dispersal",
  "Drift and Others")

process_labels <- c("Heterogeneous\nSelection", "Homogeneous\nSelection", "Dispersal\nLimitation", "Homogenizing\nDispersal",
  "Drift and\nOthers")

process_colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")

draw_icamp_month <- function() {
  dat <- read.csv(file.path(input_dir, "plot_data_month.csv"), check.names = FALSE)
  values <- as.matrix(dat[process_cols])
  values <- values/rowSums(values)
  labels <- dat$Month
  par(mar = c(6.6, 5.8, 1, 14), family = "sans", mgp = c(3.1, 0.85, 0), xaxs = "r", yaxs = "i")
  mids <- barplot(t(values), col = process_colors, border = "black", lwd = 0.45, width = 0.7, space = 0.45,
    axes = FALSE, names.arg = rep("", nrow(values)), ylim = c(0, 1.025))
  axis(1, at = mids, labels = FALSE, tick = FALSE)
  usr <- par("usr")
  text(mids, usr[3L] - 0.035 * (usr[4L] - usr[3L]), labels, srt = 60, adj = 1, xpd = NA, cex = 1)
  axis(2, at = seq(0, 1, by = 0.2), labels = paste0(seq(0, 100, by = 20), "%"), las = 1, tcl = -0.25,
    cex.axis = 1.05)
  box(bty = "o", lwd = 0.8)
  mtext("Relative importance", side = 2, line = 3.4, cex = 1.12)
  mtext("Month", side = 1, line = 5, cex = 1.12)
  legend(x = usr[2L] + 0.025 * (usr[2L] - usr[1L]), y = usr[4L], xjust = 0, yjust = 1, xpd = NA, legend = process_labels,
    fill = process_colors, border = "black", bty = "n", title = "Assembly\nmechanisms", cex = 1,
    y.intersp = 1.35, x.intersp = 0.8)
}

png_path <- file.path(output_dir, "Fig4A_iCAMP_month.png")

pdf_path <- file.path(output_dir, "Fig4A_iCAMP_month.pdf")

backup_path <- file.path(output_dir, "Fig4A_iCAMP_month_before_closed_frame.png")

if (file.exists(png_path) && !file.exists(backup_path)) {
  file.copy(png_path, backup_path, overwrite = FALSE)
}

grDevices::png(png_path, width = 12, height = 5.5, units = "in", res = 600, type = "cairo", bg = "white")

draw_icamp_month()

grDevices::dev.off()

grDevices::cairo_pdf(pdf_path, width = 12, height = 5.5, bg = "white")

draw_icamp_month()

grDevices::dev.off()

message("Updated Figure 4A: ", png_path)
