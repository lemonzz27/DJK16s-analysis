rm(list = ls())

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()

project_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/28_营养盐叶绿素蓝藻门CCF")

lag_file <- Sys.getenv("CCF_LAG_FILE", unset = file.path(project_dir, "结果表格", "09_原始水平_CCF全部lag.csv"))

peak_file <- Sys.getenv("CCF_PEAK_FILE", unset = file.path(project_dir, "结果表格", "02_CCF峰值与显著性.csv"))

figure_dir <- Sys.getenv("CCF_FIGURE_DIR", unset = file.path(project_dir, "结果图片"))

table_dir <- Sys.getenv("CCF_TABLE_DIR", unset = file.path(project_dir, "结果表格"))

single_figure_dir <- Sys.getenv("CCF_SINGLE_FIGURE_DIR", unset = file.path(figure_dir, "单图"))

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(single_figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

all_lags <- read.csv(lag_file, check.names = FALSE, stringsAsFactors = FALSE)

peaks <- read.csv(peak_file, check.names = FALSE, stringsAsFactors = FALSE)

panel_ids <- c("P01", "P02", "P04", "P09", "P15", "P12")

panel_letters <- letters[seq_along(panel_ids)]

draw_five_point_star <- function(x, y, x_radius = 0.15, y_radius = 0.07, fill = "#F0E442", border = "#7A5C00") {
  angles <- pi/2 + seq(0, 2 * pi, length.out = 11)[1:10]
  radii <- rep(c(1, 0.4), 5)
  polygon(x + x_radius * radii * cos(angles), y + y_radius * radii * sin(angles), col = fill, border = border,
    lwd = 1.35)
}

prepare_panel <- function(comparison_id) {
  z <- all_lags[all_lags$analysis_type == "level_sensitivity" & all_lags$comparison_id == comparison_id,
    , drop = FALSE]
  p <- peaks[peaks$analysis_type == "level_sensitivity" & peaks$comparison_id == comparison_id, , drop = FALSE]
  if (nrow(z) != 7L || nrow(p) != 1L) {
    stop("Invalid panel data for comparison: ", comparison_id)
  }
  synchronous <- p$lead_months == 0L || p$leading_variable == "synchronous"
  reverse_lag <- !synchronous && p$leading_variable == p$var2
  lag_value <- if (reverse_lag)
    -z$var1_lead_months
  else z$var1_lead_months
  peak_lag <- if (reverse_lag)
    -p$peak_var1_lead_months
  else p$peak_var1_lead_months
  ord <- order(lag_value)
  if (synchronous) {
    title <- paste0(p$var1, " leads ", p$var2)
    xlab <- paste0(p$var1, " lead (months)")
  }
  else {
    title <- paste0(p$leading_variable, " leads ", p$following_variable)
    xlab <- paste0(p$leading_variable, " lead (months)")
  }
  if (comparison_id == "P12") {
    lag_value <- -z$var1_lead_months
    peak_lag <- -p$peak_var1_lead_months
    ord <- order(lag_value)
    title <- "Cyanobacteria leads Chla"
    xlab <- "Cyanobacteria lead (months)"
  }
  list(lag = lag_value[ord], correlation = z$correlation[ord], peak_lag = as.numeric(peak_lag), peak_r = as.numeric(p$peak_correlation),
    title = title, xlab = xlab, p_value = as.numeric(p$max_lag_permutation_p), q_value = as.numeric(p$BH_q_global_within_analysis),
    significant = isTRUE(as.logical(p$significant)), raw_peak = p)
}

panel_data <- lapply(panel_ids, prepare_panel)

names(panel_data) <- panel_ids

plot_rows <- do.call(rbind, lapply(seq_along(panel_data), function(i) {
  d <- panel_data[[i]]
  data.frame(panel = panel_letters[i], comparison_id = panel_ids[i], relationship = d$title, oriented_lag_months = d$lag,
    correlation = d$correlation, selected_peak = d$lag == d$peak_lag, peak_permutation_p = d$p_value,
    BH_q = d$q_value, stringsAsFactors = FALSE)
}))

write.csv(plot_rows, file.path(table_dir, "12_CCF六关系统一优化图作图数据.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

draw_six_panels <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.2, 3.7, 4, 0.9), oma = c(0.3, 2, 0.2, 0.2), mgp = c(2.45, 0.68, 0),
    tcl = -0.25, family = "sans", las = 1)
  for (i in seq_along(panel_data)) {
    d <- panel_data[[i]]
    left_column <- i %in% c(1L, 4L)
    plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
      xlab = d$xlab, ylab = if (left_column)
        "CCF of log10 levels"
      else "", main = "", bty = "o", cex.lab = 1.05, cex.axis = 0.95)
    axis(1, at = -3:3, labels = -3:3, cex.axis = 0.95)
    axis(2, at = seq(-1, 1, 0.5), labels = if (left_column)
      sprintf("%.1f", seq(-1, 1, 0.5))
    else FALSE, cex.axis = 0.95)
    box(lwd = 0.85)
    abline(h = 0, v = 0, lty = 2, lwd = 0.7, col = "grey72")
    segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.6, col = "grey82")
    lines(d$lag, d$correlation, lty = 1, lwd = 1.8, col = "#00796B")
    points(d$lag, d$correlation, pch = 21, cex = 1.1, bg = "#00796B", col = "white", lwd = 0.75)
    if (d$significant) {
      segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.35, col = "#D73027")
      draw_five_point_star(d$peak_lag, d$peak_r)
    }
    mtext(d$title, side = 3, line = 2.05, adj = 0, cex = 1.05, font = 2)
    mtext(sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), side = 3, line = 0.75, adj = 0, cex = 0.95)
  }
}

png(file.path(figure_dir, "06_CCF六关系_统一优化.png"), width = 4200, height = 2400, res = 400,
  type = "cairo", bg = "white")

draw_six_panels()

dev.off()

grDevices::cairo_pdf(file.path(figure_dir, "06_CCF六关系_统一优化.pdf"), width = 10.5, height = 6,
  family = "sans")

draw_six_panels()

dev.off()

svg(file.path(figure_dir, "06_CCF六关系_统一优化.svg"), width = 10.5, height = 6, family = "sans",
  bg = "white")

draw_six_panels()

dev.off()

individual_stems <- c("06a_WL_leads_Chla", "06b_WL_leads_TN", "06c_TN_leads_Chla", "06d_WL_leads_Cyanobacteria",
  "06e_TN_leads_Cyanobacteria", "06f_Cyanobacteria_leads_Chla")

draw_single_panel <- function(d) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.8, 5, 1.2, 1), mgp = c(2.9, 0.78, 0), tcl = -0.28, family = "sans", las = 1)
  plot(d$lag, d$correlation, type = "n", xlim = c(-3.25, 3.25), ylim = c(-1, 1), xaxt = "n", yaxt = "n",
    xlab = d$xlab, ylab = "CCF of log10 levels", main = "", bty = "o", cex.lab = 1.25, cex.axis = 1.1)
  axis(1, at = -3:3, labels = -3:3)
  axis(2, at = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5)))
  box(lwd = 1)
  abline(h = 0, v = 0, lty = 2, lwd = 0.8, col = "grey72")
  segments(d$lag, 0, d$lag, d$correlation, lty = 2, lwd = 0.7, col = "grey82")
  lines(d$lag, d$correlation, lty = 1, lwd = 2.1, col = "#00796B")
  points(d$lag, d$correlation, pch = 21, cex = 1.25, bg = "#00796B", col = "white", lwd = 0.85)
  if (d$significant) {
    segments(d$peak_lag, 0, d$peak_lag, d$peak_r, lty = 2, lwd = 1.55, col = "#D73027")
    draw_five_point_star(d$peak_lag, d$peak_r, x_radius = 0.17, y_radius = 0.078)
  }
  text(-3.08, 0.94, labels = d$title, adj = c(0, 1), font = 2, cex = 1.35)
  text(-3.08, 0.78, labels = sprintf("r = %.2f, p = %.4f", d$peak_r, d$p_value), adj = c(0, 1), cex = 1.15)
}

for (i in seq_along(panel_data)) {
  png(file.path(single_figure_dir, paste0(individual_stems[i], ".png")), width = 2400, height = 1800,
    res = 400, type = "cairo", bg = "white")
  draw_single_panel(panel_data[[i]])
  dev.off()
  grDevices::cairo_pdf(file.path(single_figure_dir, paste0(individual_stems[i], ".pdf")), width = 6,
    height = 4.5, bg = "white")
  draw_single_panel(panel_data[[i]])
  dev.off()
  svg(file.path(single_figure_dir, paste0(individual_stems[i], ".svg")), width = 6, height = 4.5, bg = "white")
  draw_single_panel(panel_data[[i]])
  dev.off()
}

message("Six-panel optimized CCF figure saved to: ", figure_dir)
