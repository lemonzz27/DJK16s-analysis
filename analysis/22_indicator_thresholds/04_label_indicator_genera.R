Sys.setlocale("LC_ALL", "Chinese_China.utf8")

root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

.libPaths(c(file.path(root, "99_运行环境/R_library"), .libPaths()))

library(ggplot2)

library(grid)

stopifnot(requireNamespace("ragg", quietly = TRUE), requireNamespace("svglite", quietly = TRUE))

project <- file.path(root, "01_分析项目/36_TITAN2致病菌ASV水平_水位")

datadir <- file.path(project, "结果表格/WL/ASV/plotdata_figure2_top20_unique_species_ridgeline_latest_annotation_sunrise")

source <- file.path(Sys.getenv("DJK_CODE_ROOT", unset = "."), "analysis/22_indicator_thresholds/03_plot_indicator_ridgelines.R")

unambiguous <- "--unambiguous" %in% commandArgs(TRUE)

out <- file.path(root, if (unambiguous) "TITAN_属名缩写无歧义版_20260928" else "TITAN_属名缩写_20260928")

dir.create(out, showWarnings = FALSE)

ridge_summary <- read.csv(file.path(datadir, "07_top20_ridgeline_summary.csv"), check.names = FALSE)

curve_data <- read.csv(file.path(datadir, "08_top20_density_curve_plotdata.csv"), check.names = FALSE)

polygon_data <- read.csv(file.path(datadir, "09_top20_density_polygon_plotdata.csv"), check.names = FALSE)

original_summary <- ridge_summary

original_curve <- curve_data

original_polygon <- polygon_data

x_limits <- c(149.5, 170.5)

bubble_x <- x_limits[1] - 0.4

code <- readLines(source, encoding = "UTF-8", warn = FALSE)

start <- grep("^direction_colors <-", code)

end <- grep("^figure_width <-", code) - 1L

stopifnot(length(start) == 1L, length(end) == 1L)

eval(parse(text = code[start:end]))

original_plot <- ridge_plot + theme(legend.position = "none")

names_clean <- sub("^Candidatus\\s+", "", trimws(ridge_summary$species), ignore.case = TRUE)

parts <- strsplit(names_clean, "\\s+")

genus_names <- vapply(parts, `[`, character(1), 1)

genus_short <- paste0(substr(genus_names, 1, 1), ".")

if (unambiguous) {
  expanded <- c(Pseudomonas = "Pseud.", Paracoccus = "Par.", Mycobacterium = "Myco.", Microbacterium = "Micro.",
    Comamonas = "Com.", Curtobacterium = "Curt.")
  matched <- genus_names %in% names(expanded)
  genus_short[matched] <- unname(expanded[genus_names[matched]])
  genus_map <- unique(data.frame(genus = genus_names, abbreviation = genus_short))
  stopifnot(!anyDuplicated(genus_map$abbreviation))
  write.csv(genus_map, file.path(out, "genus_abbreviations.csv"), row.names = FALSE)
}

new_labels <- paste(genus_short, vapply(parts, `[`, character(1), 2))

stopifnot(length(new_labels) == 20L, !anyDuplicated(new_labels))

mapping <- data.frame(original_species = ridge_summary$species, original_label = ridge_summary$plot_label,
  new_label = new_labels, annotation_genus = ridge_summary$genus, label_genus = vapply(parts, `[`,
    character(1), 1))

write.csv(mapping, file.path(out, "label_mapping.csv"), row.names = FALSE, fileEncoding = "UTF-8")

revised_plot <- original_plot + scale_y_continuous(breaks = ridge_summary$Ridge_index, labels = new_labels,
  expand = expansion(add = c(0.68, 1.58)))

ragg::agg_png(file.path(out, "layout_measure.png"), width = 8, height = 11.2, units = "in", res = 100)

old_grob <- ggplotGrob(original_plot)

new_grob <- ggplotGrob(revised_plot)

extra_width <- convertWidth(sum(new_grob$widths) - sum(old_grob$widths), "in", valueOnly = TRUE)

dev.off()

w <- 8 + extra_width

h <- 11.2

prefix <- file.path(out, if (unambiguous) "TITAN_ridgeline_unambiguous" else "TITAN_ridgeline_genus_abbreviated")

ggsave(paste0(prefix, ".png"), revised_plot, width = w, height = h, dpi = 600, bg = "white", device = ragg::agg_png)

ggsave(paste0(prefix, ".pdf"), revised_plot, width = w, height = h, bg = "white", device = cairo_pdf)

ggsave(paste0(prefix, ".svg"), revised_plot, width = w, height = h, bg = "white", device = svglite::svglite)

ggsave(file.path(out, "preview.png"), revised_plot, width = w, height = h, dpi = 120, bg = "white", device = ragg::agg_png)

stopifnot(identical(original_summary, ridge_summary), identical(original_curve, curve_data), identical(original_polygon,
  polygon_data))

old_build <- ggplot_build(original_plot)

new_build <- ggplot_build(revised_plot)

stopifnot(identical(old_build$data, new_build$data))

svg <- paste(readLines(paste0(prefix, ".svg"), warn = FALSE), collapse = "\n")

stopifnot(all(vapply(new_labels, function(s) grepl(paste0(">", s, "</text>"), svg, fixed = TRUE), logical(1))))

writeLines(c("Changed: y-axis labels only; genus initials taken from original species binomials.", "Preserved: all data layers, point sizes, colors, density shapes, thresholds and order.",
  "Legends omitted to match supplied panel crop; extra width reserved only for longer labels.", paste("Size in inches:",
    round(w, 4), h), "PNG: 600 dpi; PDF/SVG: vector output."), file.path(out, "QA.txt"))

cat("PASS: 20 labels; unchanged data layers. Output:", prefix, "\n")
