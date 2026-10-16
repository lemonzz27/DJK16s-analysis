min_total_abundance <- 10

alpha_value <- 0.05

lda_threshold <- 2

max_bar_taxa <- 50

max_tree_taxa <- 150

max_tree_features <- 50

random_seed <- 1234

group_order <- c("small", "large")

group_colors <- c(large = "#D95F02", small = "#1B9E77")

input_override <- Sys.getenv("DJK16S_ACTIVE_ROOT", unset = "")

all_args <- commandArgs(trailingOnly = FALSE)

file_arg <- grep("^--file=", all_args, value = TRUE)

if (nzchar(input_override)) {
  input_dir <- normalizePath(input_override, winslash = "/", mustWork = TRUE)
} else if (length(file_arg) > 0) {
  script_file <- sub("^--file=", "", file_arg[1])
  input_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
} else {
  source_file <- NULL
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile) && length(frame$ofile) == 1 && file.exists(frame$ofile)) {
      source_file <- frame$ofile
      break
    }
  }
  if (!is.null(source_file)) {
    input_dir <- dirname(normalizePath(source_file, winslash = "/", mustWork = TRUE))
  }
  else {
    input_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
}

project_library <- file.path(input_dir, "R_library")

if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c("microeco", "ggplot2", "purrr", "tidytree", "ggtree")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "), "\n", "Install them first with:\n",
    "install.packages(c(", paste(sprintf("'%s'", missing_packages), collapse = ", "), "), repos = 'https://cloud.r-project.org')",
    call. = FALSE)
}

suppressPackageStartupMessages({
  library(microeco)
  library(ggplot2)
})

output_dir <- file.path(input_dir, "LEfSe_results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

abundance_file <- file.path(input_dir, "tpm.csv")

taxonomy_file <- file.path(input_dir, "tax.csv")

metadata_file <- file.path(input_dir, "groupsize.csv")

input_files <- c(abundance_file, taxonomy_file, metadata_file)

if (any(!file.exists(input_files))) {
  stop("Missing input file(s): ", paste(basename(input_files[!file.exists(input_files)]), collapse = ", "),
    call. = FALSE)
}

read_input_csv <- function(path) {
  x <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(x) <- trimws(sub("^\\ufeff", "", names(x)))
  x
}

set_first_column_as_rownames <- function(x, table_name) {
  ids <- trimws(as.character(x[[1]]))
  if (anyNA(ids) || any(ids == "")) {
    stop(table_name, " contains an empty ID in its first column.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    duplicate_ids <- unique(ids[duplicated(ids)])
    stop(table_name, " contains duplicated IDs: ", paste(head(duplicate_ids, 10), collapse = ", "),
      call. = FALSE)
  }
  x <- x[-1]
  rownames(x) <- ids
  x
}

asv <- set_first_column_as_rownames(read_input_csv(abundance_file), "tpm.csv")

tax <- set_first_column_as_rownames(read_input_csv(taxonomy_file), "tax.csv")

group <- read_input_csv(metadata_file)

names(group) <- tolower(trimws(names(group)))

if (!all(c("sample", "group") %in% names(group))) {
  stop("groupsize.csv must contain columns named 'sample' and 'group'.", call. = FALSE)
}

group$sample <- trimws(as.character(group$sample))

group$group <- trimws(as.character(group$group))

if (anyNA(group$sample) || any(group$sample == "") || anyNA(group$group) || any(group$group == "")) {
  stop("groupsize.csv contains an empty sample name or group name.", call. = FALSE)
}

if (anyDuplicated(group$sample)) {
  stop("groupsize.csv contains duplicated sample names.", call. = FALSE)
}

rownames(group) <- group$sample

group$sample <- NULL

input_taxonomy_columns <- names(tax)

names(tax) <- tolower(trimws(names(tax)))

rank_name_map <- c(domain = "Kingdom", kingdom = "Kingdom", phylum = "Phylum", class = "Class", order = "Order",
  family = "Family", genus = "Genus", species = "Species")

known_ranks <- intersect(names(tax), names(rank_name_map))

names(tax)[match(known_ranks, names(tax))] <- unname(rank_name_map[known_ranks])

taxonomy_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

if (anyDuplicated(names(tax))) {
  duplicate_columns <- unique(names(tax)[duplicated(names(tax))])
  stop("tax.csv maps multiple columns to the same taxonomy rank: ", paste(duplicate_columns, collapse = ", "),
    call. = FALSE)
}

missing_ranks <- setdiff(taxonomy_ranks, names(tax))

if (length(missing_ranks) > 0) {
  stop("tax.csv is missing required taxonomy rank(s): ", paste(missing_ranks, collapse = ", "), call. = FALSE)
}

excluded_taxonomy_columns <- setdiff(names(tax), taxonomy_ranks)

if (length(excluded_taxonomy_columns) > 0) {
  message("Excluding non-taxonomic column(s) from LEfSe tax_table: ", paste(excluded_taxonomy_columns,
    collapse = ", "))
}

taxonomy_columns_audit <- data.frame(input_column = input_taxonomy_columns, normalized_column = names(tax),
  LEfSe_role = ifelse(names(tax) %in% taxonomy_ranks, "taxonomy_rank", "excluded_non_taxonomy"), stringsAsFactors = FALSE)

write.csv(taxonomy_columns_audit, file.path(output_dir, "LEfSe_taxonomy_columns_audit.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

tax <- tax[, taxonomy_ranks, drop = FALSE]

tax[] <- lapply(tax, function(z) trimws(as.character(z)))

tax <- microeco::tidy_taxonomy(tax)

asv[] <- lapply(asv, function(z) suppressWarnings(as.numeric(z)))

if (anyNA(as.matrix(asv))) {
  stop("tpm.csv contains missing or non-numeric abundance values.", call. = FALSE)
}

if (any(as.matrix(asv) < 0)) {
  stop("tpm.csv contains negative abundance values.", call. = FALSE)
}

missing_samples <- setdiff(rownames(group), colnames(asv))

if (length(missing_samples) > 0) {
  stop("These metadata samples are absent from tpm.csv: ", paste(head(missing_samples, 20), collapse = ", "),
    call. = FALSE)
}

extra_samples <- setdiff(colnames(asv), rownames(group))

if (length(extra_samples) > 0) {
  message("Ignoring ", length(extra_samples), " abundance sample(s) absent from groupsize.csv.")
}

asv <- asv[, rownames(group), drop = FALSE]

missing_taxa <- setdiff(rownames(asv), rownames(tax))

if (length(missing_taxa) > 0) {
  stop("These ASVs are absent from tax.csv: ", paste(head(missing_taxa, 20), collapse = ", "), call. = FALSE)
}

keep_asv <- rowSums(asv) >= min_total_abundance

asv <- asv[keep_asv, , drop = FALSE]

tax <- tax[rownames(asv), , drop = FALSE]

if (nrow(asv) == 0) {
  stop("No ASVs remain after abundance filtering.", call. = FALSE)
}

if (any(colSums(asv) == 0)) {
  stop("At least one sample has a zero total count after filtering.", call. = FALSE)
}

observed_groups <- unique(group$group)

if (!setequal(group_order, observed_groups)) {
  message("Configured group_order does not exactly match the data; using data order: ", paste(observed_groups,
    collapse = ", "))
  group_order <- observed_groups
}

group$group <- factor(group$group, levels = group_order)

if (!all(group_order %in% names(group_colors))) {
  palette_fallback <- (grDevices::colorRampPalette(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9")))(length(group_order))
  group_colors <- stats::setNames(palette_fallback, group_order)
} else {
  group_colors <- group_colors[group_order]
}

message("Samples: ", ncol(asv), " | retained ASVs: ", nrow(asv))

message("Groups: ", paste(names(table(group$group)), as.integer(table(group$group)), sep = "=", collapse = ", "))

set.seed(random_seed)

dataset <- microtable$new(otu_table = asv, tax_table = tax, sample_table = group, auto_tidy = FALSE)

lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", taxa_level = "all", alpha = alpha_value,
  p_adjust_method = "none", lefse_norm = 1e+06)

if (is.null(lefse$res_diff) || nrow(lefse$res_diff) == 0) {
  stop("LEfSe found no significant taxa at alpha = ", alpha_value, ".", call. = FALSE)
}

write.csv(lefse$res_diff, file.path(output_dir, "LEfSe_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(lefse, file.path(output_dir, "LEfSe_microeco_object.rds"))

bar_result <- lefse$res_diff[is.finite(lefse$res_diff$LDA) & lefse$res_diff$LDA >= lda_threshold, , drop = FALSE]

if (nrow(bar_result) == 0) {
  stop("No LEfSe feature has LDA >= ", lda_threshold, ".", call. = FALSE)
}

bar_result <- bar_result[order(bar_result$LDA, decreasing = TRUE), , drop = FALSE]

bar_result <- head(bar_result, max_bar_taxa)

lefse_bar <- lefse$clone(deep = TRUE)

lefse_bar$res_diff <- bar_result

p_bar <- lefse_bar$plot_diff_bar(use_number = seq_len(nrow(bar_result)), threshold = NULL, group_order = group_order,
  color_values = unname(group_colors), color_group_map = TRUE, group_two_sep = FALSE, keep_prefix = FALSE,
  width = 0.72, ytext_size = 8) + scale_fill_manual(values = group_colors, breaks = group_order, labels = c(small = "USF",
  large = "ConF")[group_order]) + labs(x = "Taxon", y = "LDA score", fill = "Group") + theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color = "black"),
    legend.position = "right", legend.title = element_text(face = "bold"), plot.margin = margin(8,
      12, 8, 8))

ggsave(file.path(output_dir, "LEfSe_LDA_bar.png"), p_bar, width = 10, height = 11, dpi = 600, bg = "white")

ggsave(file.path(output_dir, "LEfSe_LDA_bar.pdf"), p_bar, width = 10, height = 11, device = grDevices::cairo_pdf,
  bg = "white")

p_clado <- lefse$plot_diff_cladogram(use_taxa_num = min(max_tree_taxa, nrow(tax)), use_feature_num = min(max_tree_features,
  nrow(lefse$res_diff)), clade_label_level = 1, group_order = group_order, select_show_labels = NULL,
  only_select_show = FALSE, color = unname(group_colors), sep = "|", branch_size = 0.4, alpha = 0.2,
  clade_label_size = 2.2) + theme(legend.position = "right", legend.title = element_text(face = "bold"),
  plot.margin = margin(12, 30, 12, 12))

tree_label_source <- head(lefse$res_diff, min(max_tree_features, nrow(lefse$res_diff)))

tree_label_source$node_label <- sub(".*\\|", "", tree_label_source$Taxa)

tree_label_source$display_label <- sub("^.__", "", tree_label_source$node_label)

tree_label_source <- tree_label_source[!duplicated(tree_label_source$node_label), c("node_label", "display_label",
  "Group"), drop = FALSE]

tree_node_data <- as.data.frame(p_clado$data)[, c("label", "x", "y", "angle"), drop = FALSE]

clade_labels <- merge(tree_node_data, tree_label_source, by.x = "label", by.y = "node_label", all = FALSE,
  sort = FALSE)

if (nrow(clade_labels) > 0) {
  reverse_text <- clade_labels$angle > 90 & clade_labels$angle < 270
  clade_labels$label_angle <- clade_labels$angle
  clade_labels$label_angle[reverse_text] <- clade_labels$label_angle[reverse_text] + 180
  clade_labels$label_hjust <- ifelse(reverse_text, 1, 0)
  clade_labels$label_x <- clade_labels$x + 0.12
  p_clado <- p_clado + geom_text(data = clade_labels, aes(x = label_x, y = y, label = display_label,
    angle = label_angle, hjust = label_hjust), inherit.aes = FALSE, color = "black", size = 2.2,
    family = "sans", lineheight = 0.9, show.legend = FALSE)
}

ggsave(file.path(output_dir, "LEfSe_cladogram.tiff"), p_clado, width = 18, height = 12, dpi = 600, bg = "white",
  compression = "lzw")

ggsave(file.path(output_dir, "LEfSe_cladogram.png"), p_clado, width = 18, height = 12, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "LEfSe_cladogram.pdf"), p_clado, width = 18, height = 12, device = grDevices::cairo_pdf,
  bg = "white")

summary_lines <- c(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), paste0("Input directory: ",
  input_dir), paste0("Samples: ", ncol(asv)), paste0("Retained ASVs: ", nrow(asv)), paste0("Groups: ",
  paste(names(table(group$group)), as.integer(table(group$group)), sep = "=", collapse = ", ")), paste0("Alpha: ",
  alpha_value), paste0("P-value adjustment: none"), paste0("Taxonomy ranks used: ", paste(taxonomy_ranks,
  collapse = ", ")), paste0("Excluded non-taxonomy columns: ", if (length(excluded_taxonomy_columns) >
  0) {
  paste(excluded_taxonomy_columns, collapse = ", ")
} else {
  "none"
}), paste0("LDA threshold for bar plot: ", lda_threshold), paste0("Biomarkers shown in LDA bar plot: ",
  nrow(bar_result)), paste0("Significant LEfSe rows: ", nrow(lefse$res_diff)), paste0("microeco version: ",
  as.character(utils::packageVersion("microeco"))))

writeLines(summary_lines, file.path(output_dir, "LEfSe_run_summary.txt"), useBytes = TRUE)

message("LEfSe analysis completed. Results: ", output_dir)
