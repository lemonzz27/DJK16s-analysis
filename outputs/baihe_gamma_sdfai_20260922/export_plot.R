options(stringsAsFactors = FALSE, digits = 17)

Sys.setlocale("LC_ALL", "Chinese_China.utf8")

.libPaths(c(normalizePath("Rlib", winslash = "/"), .libPaths()))

library(jsonlite)

library(grid)

boxed <- "--boxed" %in% commandArgs(trailingOnly = TRUE)

src <- "outputs/baihe_sri_matrix_20260922"

out <- "outputs/baihe_gamma_sdfai_20260922"

a <- fromJSON(file.path(src, "analysis.json"))

z <- a$indices[a$indices$method == "Gamma_MLE", ]

d <- a$data

stopifnot(identical(z$month, d$month), nrow(z) == 820)

r <- data.frame(month = z$month, previous_Q = c(NA, head(d$Q, -1)), Q = d$Q, previous_SRI = z$previous_SRI,
  SRI = z$SRI, SDFAI = z$SDFAI, event = z$SDFAI_event, endpoints = z$endpoints, status = z$status,
  source = d$source, source_cell = d$source_cell, days_valid = d$days_valid, days_expected = d$days_expected)

recalc <- (r$SRI - r$previous_SRI) * (abs(r$SRI) + abs(r$previous_SRI)) * 3.2^(-abs(r$SRI + r$previous_SRI))

stopifnot(isTRUE(all.equal(recalc, r$SDFAI, tolerance = 1e-12)))

labels <- ifelse(is.finite(r$SDFAI), ifelse(r$SDFAI > 1, "枯转丰", ifelse(r$SDFAI < -1, "丰转枯",
  "未识别")), "不可计算")

stopifnot(identical(labels, r$event), sum(is.finite(r$SDFAI)) == 506, sum(r$event == "枯转丰") ==
  19, sum(r$event == "丰转枯") == 21)

ei <- which(r$event %in% c("枯转丰", "丰转枯"))

ev <- data.frame(start_month = r$month[ei - 1], end_month = r$month[ei], previous_Q = r$previous_Q[ei],
  Q = r$Q[ei], previous_SRI = r$previous_SRI[ei], SRI = r$SRI[ei], SDFAI = r$SDFAI[ei], event = r$event[ei],
  endpoints = r$endpoints[ei])

fits <- a$fits[a$fits$method == "Gamma_MLE", ]

blocks <- data.frame(start = c("1958-01", "2006-01", "2014-01", "2024-11"), end = c("1986-12", "2012-12",
  "2018-12", "2026-04"))

mi <- function(x) as.integer(substr(x, 1, 4)) * 12 + as.integer(substr(x, 6, 7)) - 1

r$x <- mi(r$month)

blocks$lo <- mi(blocks$start) - 0.5

blocks$hi <- mi(blocks$end) + 0.5

groups <- lapply(1:4, function(i) r[r$x > blocks$lo[i] & r$x < blocks$hi[i], ])

stopifnot(sum(vapply(groups, function(x) sum(is.finite(x$SDFAI)), 1L)) == 506)

write.csv(r, file.path(out, "白河_GammaMLE_SDFAI_全部月份.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

write.csv(ev, file.path(out, "白河_GammaMLE_SDFAI_事件目录.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

write_json(list(rows = r, events = ev, fits = fits, blocks = blocks, summary = list(valid_pairs = 506,
  dry_to_wet = 19, wet_to_dry = 21), source = file.path(src, "analysis.json")), file.path(out, "selected_result.json"),
  dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 16, pretty = TRUE)

draw_figure <- function(gray = FALSE) {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", col = "#26313D")))
  drycol <- if (gray)
    "#252525"
  else "#0072B2"
  wetcol <- if (gray)
    "#6B6B6B"
  else "#D55E00"
  linecol <- "#303D4A"
  thresholdcol <- "#75828E"
  txt <- function(label, x, y, size = 8, face = "plain", just = "left", col = "#26313D") grid.text(label,
    x = x, y = y, just = just, gp = gpar(fontsize = size, fontface = face, col = col))
  txt("Baihe River: monthly SDFAI", 0.08, 0.95, 11, "bold")
  txt("Gamma MLE SRI-1 | All-data, calendar-month calibration", 0.08, 0.887, 8, col = "#596875")
  txt("506 valid month pairs", 0.945, 0.887, 8, just = "right", col = "#596875")
  xleft <- 0.08
  xright <- 0.945
  gap <- 0.017
  usable <- xright - xleft - 3 * gap
  widths <- usable * (blocks$hi - blocks$lo)/sum(blocks$hi - blocks$lo)
  lefts <- c(xleft, xleft + cumsum(widths)[1:3] + (1:3) * gap)
  bottom <- 0.27
  height <- 0.54
  ylim <- c(-4.6, 7.2)
  ticks <- list(c("1958-01", "1965-01", "1975-01", "1985-01"), c("2006-01", "2012-01"), c("2014-01",
    "2018-01"), "2025-01")
  for (i in 1:4) {
    dat <- groups[[i]]
    ok <- is.finite(dat$SDFAI)
    stopifnot(all(diff(dat$x[ok]) == 1), all(dat$SDFAI[ok] > ylim[1] & dat$SDFAI[ok] < ylim[2]))
    pushViewport(viewport(x = lefts[i], y = bottom, width = widths[i], height = height, just = c("left",
      "bottom"), xscale = c(blocks$lo[i], blocks$hi[i]), yscale = ylim, clip = "off"))
    grid.lines(x = c(0, 1), y = unit(c(0, 0), "native"), gp = gpar(col = "#CAD1D6", lwd = 0.6))
    for (h in c(-1, 1)) grid.lines(x = c(0, 1), y = unit(c(h, h), "native"), gp = gpar(col = thresholdcol,
      lwd = 0.85, lty = "dashed"))
    grid.lines(x = unit(dat$x[ok], "native"), y = unit(dat$SDFAI[ok], "native"), gp = gpar(col = linecol,
      lwd = 0.72))
    for (k in c(1, -1)) {
      e <- which(ok & if (k == 1)
        dat$SDFAI > 1
      else dat$SDFAI < -1)
      if (length(e))
        grid.points(x = unit(dat$x[e], "native"), y = unit(dat$SDFAI[e], "native"), pch = if (k ==
          1)
          24
        else 25, size = unit(1.6, "mm"), gp = gpar(col = if (k == 1)
          drycol
        else wetcol, fill = if (k == 1)
          drycol
        else "white", lwd = 0.7))
    }
    grid.lines(x = c(0, 1), y = c(0, 0), gp = gpar(col = "#26313D", lwd = 0.65))
    if (boxed) {
      grid.lines(x = c(0, 1), y = c(1, 1), gp = gpar(col = "#26313D", lwd = 0.65))
      if (i == 1)
        grid.lines(x = c(0, 0), y = c(0, 1), gp = gpar(col = "#26313D", lwd = 0.65))
      if (i == 4)
        grid.lines(x = c(1, 1), y = c(0, 1), gp = gpar(col = "#26313D", lwd = 0.65))
    }
    grid.xaxis(at = mi(ticks[[i]]), label = substr(ticks[[i]], 1, 4), gp = gpar(fontsize = 7, lwd = 0.65,
      col = "#26313D"))
    if (i == 1)
      grid.yaxis(at = c(-4, -2, 0, 2, 4, 6), label = as.character(c(-4, -2, 0, 2, 4, 6)), gp = gpar(fontsize = 8,
        lwd = 0.65, col = "#26313D"))
    popViewport()
  }
  for (i in 1:3) {
    centre <- lefts[i] + widths[i] + gap/2
    for (offset in c(-0.0025, 0.0025)) grid.lines(x = c(centre + offset - 0.002, centre + offset +
      0.002), y = c(bottom - 0.011, bottom + 0.011), gp = gpar(col = "#26313D", lwd = 0.8))
  }
  if (boxed)
    for (i in 1:3) {
      centre <- lefts[i] + widths[i] + gap/2
      for (offset in c(-0.0025, 0.0025)) grid.lines(x = c(centre + offset - 0.002, centre + offset +
        0.002), y = c(bottom + height - 0.011, bottom + height + 0.011), gp = gpar(col = "#26313D",
        lwd = 0.8))
    }
  grid.text("SDFAI", x = 0.025, y = bottom + height/2, rot = 90, gp = gpar(fontsize = 9))
  for (h in c(-1, 1)) txt(if (h > 0)
    "+1"
  else "-1", 0.955, bottom + height * (h - ylim[1])/diff(ylim), 7.5, col = thresholdcol)
  txt("Year", 0.51, 0.163, 8, just = "center")
  grid.points(x = 0.105, y = 0.108, pch = 24, size = unit(1.8, "mm"), gp = gpar(col = drycol, fill = drycol))
  txt("Dry-to-wet (19)", 0.12, 0.108, 8)
  grid.points(x = 0.35, y = 0.108, pch = 25, size = unit(1.8, "mm"), gp = gpar(col = wetcol, fill = "white"))
  txt("Wet-to-dry (21)", 0.365, 0.108, 8)
  grid.lines(x = c(0.607, 0.642), y = c(0.108, 0.108), gp = gpar(col = thresholdcol, lty = "dashed",
    lwd = 0.85))
  txt("Thresholds: >1 and <-1", 0.655, 0.108, 8)
  txt("Breaks: 1987-2005, 2013, Jan 2019-Oct 2024. Index assigned to the later month; missing data are not joined.",
    0.08, 0.038, 6.8, col = "#596875")
  popViewport()
}

stem <- file.path(out, "白河_GammaMLE_SDFAI_全时期断轴时间曲线")

if (boxed) stem <- paste0(stem, "_封闭边框")

png(paste0(stem, ".png"), width = 9, height = 4.1, units = "in", res = 600, type = "cairo")

draw_figure()

dev.off()

tiff(paste0(stem, ".tiff"), width = 9, height = 4.1, units = "in", res = 600, compression = "lzw", type = "cairo")

draw_figure()

dev.off()

cairo_pdf(paste0(stem, ".pdf"), width = 9, height = 4.1, family = "Arial")

draw_figure()

dev.off()

svglite::svglite(paste0(stem, ".svg"), width = 9, height = 4.1)

draw_figure()

dev.off()

png(file.path(out, if (boxed) "figure_grayscale_boxed_check.png" else "figure_grayscale_check.png"),
  width = 9, height = 4.1, units = "in", res = 180, type = "cairo")

draw_figure(TRUE)

dev.off()

print(range(r$SDFAI, na.rm = TRUE))

print(ev[ev$end_month >= "2024-11", ])
