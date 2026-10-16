rm(list = ls())

ROOT <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

PROJECT <- file.path(ROOT, "01_分析项目", "28_营养盐叶绿素蓝藻门CCF")

TABLE_ROOT <- file.path(PROJECT, "结果表格")

FIGURE_DIR <- Sys.getenv("CCF_SUPPLEMENTARY_SINGLE_DIR", unset = file.path(PROJECT, "结果图片", "单图"))

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

draw_five_point_star <- function(x, y, x_radius = 0.17, y_radius = 0.078, fill = "#F0E442", border = "#7A5C00") {
  angles <- pi/2 + seq(0, 2 * pi, length.out = 11)[1:10]
  radii <- rep(c(1, 0.4), 5)
  polygon(x + x_radius * radii * cos(angles), y + y_radius * radii * sin(angles), col = fill, border = border,
    lwd = 1.35)
}

load_panel <- function(subdir, variable, title, xlab, ylab, filename, title_cex = 1.1, other_variable = "Cyanobacteria") {
  table_dir <- file.path(TABLE_ROOT, subdir)
  lag_data <- read_csv(file.path(table_dir, "01_新增变量原始水平_CCF全部lag.csv"))
  peak_data <- read_csv(file.path(table_dir, "02_新增变量原始水平_CCF峰值_未校正p.csv"))
  z <- lag_data[lag_data$var1 == variable & lag_data$var2 == other_variable & lag_data$analysis_group %in%
    c("overall", "environment"), , drop = FALSE]
  p <- peak_data[peak_data$var1 == variable & peak_data$var2 == other_variable & peak_data$analysis_group %in%
    c("overall", "environment"), , drop = FALSE]
  if (nrow(z) != 7L || nrow(p) != 1L) {
    stop("Expected seven lag rows and one peak row for ", variable, " vs ", other_variable)
  }
  ord <- order(z$var1_lead_months)
  list(lag = as.numeric(z$var1_lead_months[ord]), correlation = as.numeric(z$correlation[ord]), peak_lag = as.numeric(p$peak_var1_lead_months),
    peak_r = as.numeric(p$peak_correlation), p_value = as.numeric(p$max_lag_permutation_p_unadjusted),
    significant = isTRUE(as.logical(p$significant_uncorrected)), title = title, xlab = xlab, ylab = ylab,
    filename = filename, title_cex = title_cex)
}

panels <- list(load_panel("补充_TOC_WT_原始水平", "WT", "Water temperature leads Cyanobacteria",
  "Water temperature lead (months)", "CCF of log10 levels", "06g_Water_temperature_leads_Cyanobacteria.png",
  title_cex = 1.22), load_panel("补充_气温AIRT_原始水平", "AIRT", "Air temperature leads Cyanobacteria",
  "Air temperature lead (months)", "CCF of original levels", "06h_Air_temperature_leads_Cyanobacteria.png",
  title_cex = 1.1), load_panel("补充_氮磷比TNTP_原始水平", "TNTP", "TN:TP mass ratio leads Chla",
  "TN:TP mass ratio lead (months)", "CCF of log10 levels", "06i_TNTP_mass_ratio_leads_Chla.png", title_cex = 1.22,
  other_variable = "Chla"))

draw_single_panel <- function(d) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.8, 5, 1.2, 1), mgp = c(2.9, 0.78, 0), tcl = -0.28, family = "sans", las = 1)
  plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
    xlab = d$xlab, ylab = d$ylab, main = "", bty = "o", cex.lab = 1.25, cex.axis = 1.1)
  axis(1, at = -3:3, labels = -3:3)
  axis(2, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)))
  box(lwd = 1)
  abline(h = 0, v = 0, lty = 2, lwd = 0.8, col = "grey72")
  segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.7, col = "grey82")
  lines(d$lag, d$correlation, lty = 1, lwd = 2.1, col = "#00796B")
  points(d$lag, d$correlation, pch = 21, cex = 1.25, bg = "#00796B", col = "white", lwd = 0.85)
  if (d$significant) {
    segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.55, col = "#D73027")
    draw_five_point_star(d$peak_lag, d$peak_r)
  }
  text(-3.08, 0.94, labels = d$title, adj = c(0, 1), font = 2, cex = d$title_cex)
  text(-3.08, 0.78, labels = sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), adj = c(0, 1), cex = 1.15)
}

for (d in panels) {
  png(file.path(FIGURE_DIR, d$filename), width = 2400, height = 1800, res = 400, type = "cairo", bg = "white")
  draw_single_panel(d)
  dev.off()
  grDevices::cairo_pdf(file.path(FIGURE_DIR, sub("\\.png$", ".pdf", d$filename)), width = 6, height = 4.5,
    bg = "white")
  draw_single_panel(d)
  dev.off()
  svg(file.path(FIGURE_DIR, sub("\\.png$", ".svg", d$filename)), width = 6, height = 4.5, bg = "white")
  draw_single_panel(d)
  dev.off()
}

message("Styled supplementary CCF single panels saved to: ", FIGURE_DIR)
