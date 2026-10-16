rm(list = ls())

options(stringsAsFactors = FALSE, scipen = 999)

try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"), silent = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L) {
  stop("Usage: Rscript script.R <source_csv> <figure_dir> <table_dir>")
}

Encoding(args) <- "UTF-8"

SOURCE_CSV <- args[[1L]]

FIGURE_DIR <- args[[2L]]

TABLE_DIR <- args[[3L]]

utf8_path <- function(x) {
  Encoding(x) <- "UTF-8"
  x
}

if (!file.exists(SOURCE_CSV)) {
  stop("Missing source table: ", SOURCE_CSV)
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

source_data <- read.csv(SOURCE_CSV, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

required_columns <- c("group", "category", "mean_relative_abundance_percent")

if (!all(required_columns %in% names(source_data))) {
  stop("Source table is missing required columns.")
}

plot_data <- source_data[, required_columns]

plot_data$group <- tolower(trimws(as.character(plot_data$group)))

plot_data$category <- trimws(as.character(plot_data$category))

plot_data$mean_relative_abundance_percent <- suppressWarnings(as.numeric(plot_data$mean_relative_abundance_percent))

expected_keys <- c("large|Non-pathogen", "large|Pathogen", "small|Non-pathogen", "small|Pathogen")

observed_keys <- paste(plot_data$group, plot_data$category, sep = "|")

if (nrow(plot_data) != 4L || !setequal(observed_keys, expected_keys) || anyNA(plot_data$mean_relative_abundance_percent) ||
  any(!is.finite(plot_data$mean_relative_abundance_percent))) {
  stop("The relative-abundance source data failed structural validation.")
}

lookup_value <- function(group_name, category_name) {
  plot_data$mean_relative_abundance_percent[plot_data$group == group_name & plot_data$category == category_name][[1L]]
}

group_order <- c("large", "small")

nonpathogen <- vapply(group_order, lookup_value, numeric(1), category_name = "Non-pathogen")

pathogen <- vapply(group_order, lookup_value, numeric(1), category_name = "Pathogen")

if (any(abs(nonpathogen + pathogen - 100) > 1e-06)) {
  stop("Relative-abundance percentages do not sum to 100% within each group.")
}

values <- rbind(`Non-pathogen` = nonpathogen, `Pathogen candidate` = pathogen)

colnames(values) <- c("Large\nn = 60", "Small\nn = 60")

category_colors <- c("#A7B6C2", "#D86648")

text_colors <- c("#1F2933", "white")

width_mm <- 78

height_mm <- 112

width_in <- width_mm/25.4

height_in <- height_mm/25.4

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(family = "Arial", mar = c(4.4, 5, 6, 1), mgp = c(2.8, 0.75, 0), tcl = -0.28, las = 1, xpd = NA,
    bty = "l", col.axis = "#202020", col.lab = "#202020", fg = "#202020")
  mids <- barplot(values, beside = FALSE, col = category_colors, border = "white", lwd = 0.75, width = c(0.31,
    0.31), space = 2.3, ylim = c(0, 106), axes = FALSE, names.arg = colnames(values), cex.names = 0.92,
    font.axis = 1)
  axis(side = 2, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), cex.axis = 0.9, lwd = 0.85,
    lwd.ticks = 0.85)
  axis(side = 1, at = mids, labels = FALSE, lwd = 0.85, lwd.ticks = 0.85)
  mtext("Mean relative abundance (%)", side = 2, line = 3.3, las = 0, cex = 1)
  for (i in seq_along(mids)) {
    text(mids[i], nonpathogen[i]/2, labels = sprintf("%.1f%%", nonpathogen[i]), col = text_colors[1L],
      font = 2, cex = 0.94)
    text(mids[i], nonpathogen[i] + pathogen[i]/2, labels = sprintf("%.1f%%", pathogen[i]), col = text_colors[2L],
      font = 2, cex = 0.94)
    text(mids[i], 102, labels = "100%", col = "#202020", cex = 0.84)
  }
  title(main = "Relative abundance", line = 4.25, font.main = 2, cex.main = 1.12)
  mtext("Rarefied to 21,000 reads per sample", side = 3, line = 2.95, cex = 0.8, col = "#555555")
  legend(x = mean(range(mids)), y = 112, legend = rownames(values), fill = category_colors, border = NA,
    horiz = TRUE, bty = "n", xjust = 0.5, yjust = 0.5, cex = 0.77, x.intersp = 0.45, text.width = c(0.82,
      1.23))
}

file_stem <- utf8_path(file.path(FIGURE_DIR, "pathogen_candidate_relative_abundance_large_small_portrait"))

png_path <- utf8_path(paste0(file_stem, ".png"))

pdf_path <- utf8_path(paste0(file_stem, ".pdf"))

svg_path <- utf8_path(paste0(file_stem, ".svg"))

tiff_path <- utf8_path(paste0(file_stem, ".tiff"))

grDevices::png(png_path, width = width_in, height = height_in, units = "in", res = 600, pointsize = 9,
  type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))

draw_figure()

grDevices::dev.off()

grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial", pointsize = 9,
  onefile = TRUE)

draw_figure()

grDevices::dev.off()

grDevices::svg(svg_path, width = width_in, height = height_in, family = "Arial", pointsize = 9, onefile = TRUE)

draw_figure()

grDevices::dev.off()

grDevices::tiff(tiff_path, width = width_in, height = height_in, units = "in", res = 600, pointsize = 9,
  compression = "lzw", type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))

draw_figure()

grDevices::dev.off()

plot_data <- plot_data[order(match(plot_data$group, group_order), match(plot_data$category, c("Non-pathogen",
  "Pathogen"))), ]

plot_data_path <- utf8_path(file.path(TABLE_DIR, "pathogen_candidate_relative_abundance_large_small_plot_data.csv"))

write.csv(plot_data, plot_data_path, row.names = FALSE, fileEncoding = "UTF-8")

qa_path <- utf8_path(file.path(TABLE_DIR, "pathogen_candidate_relative_abundance_large_small_QA.txt"))

capture.output({
  cat("Single-panel pathogen-candidate relative-abundance figure\n")
  cat("========================================================\n")
  cat("Source: ", normalizePath(SOURCE_CSV, winslash = "/"), "\n", sep = "")
  cat("Source MD5: ", unname(tools::md5sum(SOURCE_CSV)), "\n", sep = "")
  cat("Rarefaction depth: 21,000 reads per sample\n")
  cat("Large samples: n = 60\n")
  cat("Small samples: n = 60\n")
  cat("Pathogen candidate mean relative abundance: large = ", sprintf("%.4f%%", pathogen[[1L]]), "; small = ",
    sprintf("%.4f%%", pathogen[[2L]]), "\n", sep = "")
  cat("Figure size: 78 mm x 112 mm; raster resolution: 600 dpi\n")
  cat("Bars start at zero and sum to 100% within each size group.\n")
}, file = qa_path)

message("Completed: ", normalizePath(FIGURE_DIR, winslash = "/"))
