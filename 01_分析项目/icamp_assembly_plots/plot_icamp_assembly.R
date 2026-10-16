options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

all_args <- commandArgs(trailingOnly = FALSE)

script_arg <- grep("^--file=", all_args, value = TRUE)

script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_name <- "DJK16s_bin48_pairwise_process_importance.csv"

input_candidates <- c(file.path(dirname(script_dir), input_name), file.path(dirname(dirname(script_dir)),
  input_name), file.path(getwd(), input_name))

existing_input <- input_candidates[file.exists(input_candidates)]

default_input <- if (length(existing_input) > 0L) existing_input[1L] else input_candidates[1L]

input_file <- if (length(args) >= 1L) args[1L] else default_input

output_dir <- if (length(args) >= 2L) args[2L] else script_dir

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

input_file <- normalizePath(input_file, winslash = "/")

output_dir <- normalizePath(output_dir, winslash = "/")

process_cols <- c("Heterogeneous.Selection", "Homogeneous.Selection", "Dispersal.Limitation", "Homogenizing.Dispersal",
  "Drift.and.Others")

process_labels <- c("Heterogeneous Selection", "Homogeneous Selection", "Dispersal Limitation", "Homogenizing Dispersal",
  "Drift and Others")

process_colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")

size_order <- c("ConF", "USF")

parse_sample_metadata <- function(samples) {
  cleaned <- gsub("_", "", trimws(samples), fixed = TRUE)
  pattern <- "^([A-Za-z]+)([0-9]{2})([0-9]{2})(0[12])$"
  matches <- regexec(pattern, cleaned)
  pieces <- regmatches(cleaned, matches)
  valid <- lengths(pieces) == 5L
  if (!all(valid)) {
    bad <- unique(samples[!valid])
    stop("Could not parse sample name(s): ", paste(utils::head(bad, 8L), collapse = ", "), ". Expected <site><YY><MM><01|02>.")
  }
  values <- do.call(rbind, lapply(pieces, function(x) x[2:5]))
  month_number <- as.integer(values[, 3L])
  if (any(month_number < 1L | month_number > 12L)) {
    stop("Invalid month code found in sample names.")
  }
  size_code <- values[, 4L]
  data.frame(Sample = samples, Site = values[, 1L], Month = paste0("20", values[, 2L], "-", values[,
    3L]), MonthDate = as.Date(paste0("20", values[, 2L], "-", values[, 3L], "-01")), Size = ifelse(size_code ==
    "01", "USF", "ConF"), SizeCode = size_code, check.names = FALSE)
}

aggregate_processes <- function(data, group_cols) {
  means <- aggregate(data[process_cols], data[group_cols], mean, na.rm = TRUE)
  counts <- aggregate(data[[process_cols[1L]]], data[group_cols], length)
  names(counts)[ncol(counts)] <- "n_pairs"
  merge(means, counts, by = group_cols, sort = FALSE)
}

write_summary <- function(data, filename) {
  out <- data
  names(out)[match(process_cols, names(out))] <- process_labels
  utils::write.csv(out, file.path(output_dir, filename), row.names = FALSE, fileEncoding = "UTF-8")
}

draw_stacked_bars <- function(data, labels, main = NULL, xlab = NULL, show_legend = TRUE, name_cex = 0.82,
  label_angle = 45, xlab_line = 4.8) {
  values <- as.matrix(data[process_cols])
  storage.mode(values) <- "double"
  values <- values/rowSums(values)
  display_labels <- labels
  if (length(display_labels) > 8L && abs(label_angle)%%180 != 90) {
    display_labels[seq_along(display_labels)%%2L == 0L] <- ""
  }
  mids <- barplot(t(values), col = process_colors, border = "black", lwd = 0.45, space = 0.28, axes = FALSE,
    names.arg = rep("", nrow(values)), ylim = c(0, 1), main = main, cex.main = 1.05, font.main = 2)
  if (label_angle%%180 == 0) {
    axis(1, at = mids, labels = display_labels, las = 1, tick = FALSE, line = -0.35, cex.axis = name_cex)
  }
  else if (abs(label_angle)%%180 == 90) {
    axis(1, at = mids, labels = display_labels, las = 2, tick = FALSE, line = -0.35, cex.axis = name_cex)
  }
  else {
    axis(1, at = mids, labels = FALSE, tick = FALSE)
    user_coords <- par("usr")
    label_y <- user_coords[3L] - 0.035 * (user_coords[4L] - user_coords[3L])
    text(x = mids, y = label_y, labels = display_labels, srt = label_angle, adj = 1, xpd = NA, cex = name_cex)
  }
  axis(2, at = seq(0, 1, by = 0.2), labels = paste0(seq(0, 100, by = 20), "%"), las = 1, tcl = -0.25,
    cex.axis = 1.05)
  box(bty = "l", lwd = 0.8)
  mtext("Relative importance", side = 2, line = 3.2, cex = 1.12, las = 0)
  if (!is.null(xlab)) {
    mtext(xlab, side = 1, line = xlab_line, cex = 1.12)
  }
  if (show_legend) {
    user_coords <- par("usr")
    plot_width <- user_coords[2L] - user_coords[1L]
    legend(x = user_coords[2L] + 0.025 * plot_width, y = user_coords[4L], xjust = 0, yjust = 1, xpd = NA,
      legend = process_labels, fill = process_colors, border = "black", bty = "n", title = "Assembly mechanisms",
      cex = 1, y.intersp = 1.1)
  }
  invisible(mids)
}

save_figure <- function(basename, width, height, draw_function) {
  pdf_path <- file.path(output_dir, paste0(basename, ".pdf"))
  png_path <- file.path(output_dir, paste0(basename, ".png"))
  grDevices::pdf(pdf_path, width = width, height = height, family = "Helvetica", useDingbats = FALSE,
    bg = "white")
  draw_function()
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = 600, bg = "white")
  draw_function()
  grDevices::dev.off()
}

message("Reading: ", input_file)

pairwise <- utils::read.csv(input_file, check.names = FALSE)

required_cols <- c("sample1", "sample2", process_cols)

missing_cols <- setdiff(required_cols, names(pairwise))

if (length(missing_cols) > 0L) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

if (anyNA(pairwise[required_cols])) {
  stop("Missing values were found in required columns.")
}

for (column in process_cols) {
  pairwise[[column]] <- as.numeric(pairwise[[column]])
}

process_matrix <- as.matrix(pairwise[process_cols])

if (any(!is.finite(process_matrix))) {
  stop("Non-finite process importance values were found.")
}

if (any(process_matrix < -1e-10 | process_matrix > 1 + 1e-10)) {
  stop("Process importance values must lie between 0 and 1.")
}

row_total_error <- max(abs(rowSums(process_matrix) - 1))

if (row_total_error > 1e-08) {
  stop("Process importance rows do not sum to 1; max error = ", row_total_error)
}

meta1 <- parse_sample_metadata(pairwise$sample1)

meta2 <- parse_sample_metadata(pairwise$sample2)

month_pairs <- data.frame(Month1 = meta1$Month, Month2 = meta2$Month, pairwise[process_cols], check.names = FALSE)

month_pairs <- month_pairs[month_pairs$Month1 == month_pairs$Month2, , drop = FALSE]

month_pairs$Month <- month_pairs$Month1

monthly <- aggregate_processes(month_pairs, "Month")

monthly$MonthDate <- as.Date(paste0(monthly$Month, "-01"))

monthly <- monthly[order(monthly$MonthDate), c("Month", process_cols, "n_pairs")]

month_size_pairs <- data.frame(Month1 = meta1$Month, Month2 = meta2$Month, Size1 = meta1$Size, Size2 = meta2$Size,
  pairwise[process_cols], check.names = FALSE)

month_size_pairs <- month_size_pairs[month_size_pairs$Month1 == month_size_pairs$Month2 & month_size_pairs$Size1 ==
  month_size_pairs$Size2, , drop = FALSE]

month_size_pairs$Month <- month_size_pairs$Month1

month_size_pairs$Size <- month_size_pairs$Size1

monthly_size <- aggregate_processes(month_size_pairs, c("Month", "Size"))

monthly_size$MonthDate <- as.Date(paste0(monthly_size$Month, "-01"))

monthly_size$SizeOrder <- match(monthly_size$Size, size_order)

monthly_size <- monthly_size[order(monthly_size$SizeOrder, monthly_size$MonthDate), c("Month", "Size",
  process_cols, "n_pairs")]

size_pairs <- data.frame(Size1 = meta1$Size, Size2 = meta2$Size, pairwise[process_cols], check.names = FALSE)

size_pairs <- size_pairs[size_pairs$Size1 == size_pairs$Size2, , drop = FALSE]

size_pairs$Size <- size_pairs$Size1

by_size <- aggregate_processes(size_pairs, "Size")

by_size$SizeOrder <- match(by_size$Size, size_order)

by_size <- by_size[order(by_size$SizeOrder), c("Size", process_cols, "n_pairs")]

write_summary(monthly, "plot_data_month.csv")

write_summary(monthly_size, "plot_data_month_by_size.csv")

write_summary(by_size, "plot_data_size.csv")

save_figure("iCAMP_month", width = 12, height = 5.5, function() {
  par(mar = c(6.3, 5, 2.5, 13), family = "sans", las = 1, mgp = c(2.5, 0.7, 0), xaxs = "r", yaxs = "i")
  draw_stacked_bars(monthly, labels = monthly$Month, main = "Monthly Assembly Mechanisms (All Sizes)",
    xlab = "Month", show_legend = TRUE, name_cex = 1.08)
})

save_figure("iCAMP_month_by_size", width = 12, height = 8.4, function() {
  par(mfrow = c(2, 1), mar = c(5.8, 5, 2.3, 13), oma = c(0.2, 0, 0.3, 0), family = "sans", las = 1,
    mgp = c(2.5, 0.7, 0), xaxs = "r", yaxs = "i")
  for (i in seq_along(size_order)) {
    current_size <- size_order[i]
    panel <- monthly_size[monthly_size$Size == current_size, , drop = FALSE]
    panel <- panel[order(as.Date(paste0(panel$Month, "-01"))), , drop = FALSE]
    draw_stacked_bars(panel, labels = panel$Month, main = current_size, xlab = "Month", show_legend = (i ==
      1L), name_cex = 1.02)
  }
})

save_figure("iCAMP_size", width = 7.4, height = 5.2, function() {
  par(mar = c(4.2, 5, 2.5, 13), family = "sans", las = 1, mgp = c(2.5, 0.7, 0), xaxs = "r", yaxs = "i")
  draw_stacked_bars(by_size, labels = by_size$Size, main = "Assembly Mechanisms by Size", xlab = "Size",
    show_legend = TRUE, name_cex = 0.9, label_angle = 0, xlab_line = 2.7)
})

all_samples <- unique(c(pairwise$sample1, pairwise$sample2))

all_metadata <- parse_sample_metadata(all_samples)

message("Completed successfully.")

message("  Pairwise rows: ", nrow(pairwise))

message("  Unique samples: ", length(all_samples))

message("  Months: ", length(unique(all_metadata$Month)), " (", min(all_metadata$Month), " to ", max(all_metadata$Month),
  ")")

message("  Sizes: ", paste(size_order, collapse = ", "))

message("  Maximum row-sum error: ", format(row_total_error, scientific = TRUE))

message("  Output directory: ", output_dir)
