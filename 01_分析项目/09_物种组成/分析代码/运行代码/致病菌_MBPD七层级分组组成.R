options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c("dplyr", "ggplot2", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

project_root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

hits_file <- file.path(project_root, "00_输入文件/04_分类注释/MBPD/02_filtered_hits_with_taxonomy.tsv")

abundance_file <- file.path(project_root, "00_输入文件/05_致病菌候选/完整群落抽平21000后候选ASV子表.csv")

complete_rarefied_file <- file.path(project_root, "00_输入文件/01_ASV与序列/ASV抽平21000_36097ASV_120样本.csv")

group_files <- list(Size = file.path(project_root, "00_输入文件/02_样本分组与站点/孔径分组_groupsize.csv"),
  Location = file.path(project_root, "00_输入文件/02_样本分组与站点/地点分组_grouplocation.csv"),
  Season = file.path(project_root, "00_输入文件/02_样本分组与站点/季节分组_groupseason.csv"),
  Rainfall = file.path(project_root, "00_输入文件/02_样本分组与站点/降雨分组_grouprain.csv"))

required_files <- c(hits_file, abundance_file, complete_rarefied_file, unlist(group_files))

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

output_dir <- file.path(project_root, "outputs/20260918_candidate_rarefied_unified/analysis")

plot_dir <- file.path(output_dir, "figures")

table_dir <- file.path(output_dir, "tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 15L

taxonomic_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")

group_contract <- list(Size = list(order = c("large", "small"), labels = c(large = "ConF", small = "USF")),
  Location = list(order = c("BS", "KX", "LHKX", "TC"), labels = c(BS = "BS", KX = "KX", LHKX = "LHKX",
    TC = "TC")), Season = list(order = c("spring", "summer", "autumn", "winter"), labels = c(spring = "Spring",
    summer = "Summer", autumn = "Autumn", winter = "Winter")), Rainfall = list(order = c("dry", "wet"),
    labels = c(dry = "Dry", wet = "Wet")))

read_group_file <- function(path, category) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if (ncol(x) < 2L)
    stop("Invalid group file: ", path)
  names(x)[1:2] <- c("SampleID", "GroupRaw")
  x <- x[, c("SampleID", "GroupRaw"), drop = FALSE]
  x$SampleID <- trimws(as.character(x$SampleID))
  x$GroupRaw <- trimws(as.character(x$GroupRaw))
  x$Category <- category
  contract <- group_contract[[category]]
  unknown <- setdiff(unique(x$GroupRaw), names(contract$labels))
  if (length(unknown) > 0L) {
    stop("Unexpected group labels in ", path, ": ", paste(unknown, collapse = ", "))
  }
  x$Group <- unname(contract$labels[x$GroupRaw])
  x
}

group_metadata <- lapply(names(group_files), function(category) {
  read_group_file(group_files[[category]], category)
})

names(group_metadata) <- names(group_files)

category_levels <- names(group_contract)

group_levels <- unlist(lapply(category_levels, function(category) {
  contract <- group_contract[[category]]
  unname(contract$labels[contract$order])
}), use.names = FALSE)

hits <- read.delim(hits_file, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
  stringsAsFactors = FALSE)

required_hit_columns <- c("qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "qcovs",
  "taxonomy")

if (!all(required_hit_columns %in% names(hits))) {
  stop("MBPD hit table is missing columns: ", paste(setdiff(required_hit_columns, names(hits)), collapse = ", "))
}

numeric_columns <- c("pident", "length", "qlen", "evalue", "bitscore", "qcovs")

for (column in numeric_columns) hits[[column]] <- as.numeric(hits[[column]])

if (anyNA(hits[, numeric_columns])) stop("MBPD hit metrics contain missing/non-numeric values.")

numeric_best_ties <- hits %>%
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>%
  filter(pident == max(pident)) %>%
  filter(qcovs == max(qcovs)) %>%
  filter(length == max(length)) %>%
  filter(evalue == min(evalue)) %>%
  ungroup()

top_bitscore_audit <- hits %>%
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>%
  summarise(Top_Bitscore_Hit_Rows = n(), Top_Bitscore_Distinct_Taxonomies = n_distinct(taxonomy), .groups = "drop")

numeric_tie_audit <- numeric_best_ties %>%
  group_by(qseqid) %>%
  summarise(Final_Numeric_Tie_Rows = n(), Final_Numeric_Tie_Distinct_Taxonomies = n_distinct(taxonomy),
    .groups = "drop")

all_hit_audit <- hits %>%
  count(qseqid, name = "All_Filtered_Hit_Rows")

best_hits <- numeric_best_ties %>%
  arrange(qseqid, sseqid) %>%
  group_by(qseqid) %>%
  slice_head(n = 1L) %>%
  ungroup() %>%
  left_join(all_hit_audit, by = "qseqid") %>%
  left_join(top_bitscore_audit, by = "qseqid") %>%
  left_join(numeric_tie_audit, by = "qseqid")

rank_prefix <- c(domain = "d__", phylum = "p__", class = "c__", order = "o__", family = "f__", genus = "g__",
  species = "s__")

parse_taxonomy <- function(text) {
  tokens <- trimws(strsplit(text, ";", fixed = TRUE)[[1]])
  result <- setNames(rep(NA_character_, length(rank_prefix)), names(rank_prefix))
  for (rank in names(rank_prefix)) {
    token <- tokens[startsWith(tokens, rank_prefix[[rank]])]
    if (length(token) > 0L) {
      value <- trimws(sub(paste0("^", rank_prefix[[rank]]), "", token[1]))
      if (nzchar(value))
        result[[rank]] <- value
    }
  }
  result
}

parsed <- t(vapply(best_hits$taxonomy, parse_taxonomy, FUN.VALUE = setNames(rep(NA_character_, length(rank_prefix)),
  names(rank_prefix))))

parsed <- as.data.frame(parsed, stringsAsFactors = FALSE)

best_hits_ranked <- bind_cols(best_hits, parsed)

is_informative <- function(x) {
  !is.na(x) & nzchar(trimws(x)) & !grepl("uncultured|unclassified|unassigned|unknown|metagenome|environmental sample|^bacterium$",
    x, ignore.case = TRUE)
}

taxonomy <- data.frame(ASV_ID = best_hits_ranked$qseqid, All_Filtered_Hit_Rows = best_hits_ranked$All_Filtered_Hit_Rows,
  Top_Bitscore_Hit_Rows = best_hits_ranked$Top_Bitscore_Hit_Rows, Top_Bitscore_Distinct_Taxonomies = best_hits_ranked$Top_Bitscore_Distinct_Taxonomies,
  Final_Numeric_Tie_Rows = best_hits_ranked$Final_Numeric_Tie_Rows, Final_Numeric_Tie_Distinct_Taxonomies = best_hits_ranked$Final_Numeric_Tie_Distinct_Taxonomies,
  Best_Hit_Rows = 1L, Best_Distinct_Taxonomies = 1L, Selected_sseqid = best_hits_ranked$sseqid, Selection_Rule = paste0("max_bitscore_then_max_pident_then_max_qcovs_then_max_length_",
    "then_min_evalue_then_lexical_sseqid"), Best_Bitscore = best_hits_ranked$bitscore, Best_Pident = best_hits_ranked$pident,
  Best_Qcovs = best_hits_ranked$qcovs, Best_Alignment_Length = best_hits_ranked$length, stringsAsFactors = FALSE,
  check.names = FALSE)

for (rank in taxonomic_levels) {
  values <- trimws(as.character(best_hits_ranked[[rank]]))
  taxonomy[[rank]] <- ifelse(is_informative(values), values, NA_character_)
}

taxonomy <- taxonomy %>%
  arrange(ASV_ID)

abundance_raw <- read.csv(abundance_file, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

names(abundance_raw)[1] <- "ASV_ID"

abundance_raw$ASV_ID <- as.character(abundance_raw$ASV_ID)

if (anyDuplicated(abundance_raw$ASV_ID)) stop("Duplicated ASV IDs in abundance table.")

abundance_matrix <- matrix(as.numeric(as.matrix(abundance_raw[, -1, drop = FALSE])), nrow = nrow(abundance_raw),
  ncol = ncol(abundance_raw) - 1L, dimnames = list(abundance_raw$ASV_ID, names(abundance_raw)[-1]))

if (anyNA(abundance_matrix) || any(abundance_matrix < 0)) {
  stop("Abundance table contains missing, non-numeric, or negative values.")
}

complete_rarefied <- read.csv(complete_rarefied_file, row.names = 1, check.names = FALSE)

stopifnot(all(colSums(complete_rarefied) == 21000))

expected_subset <- as.matrix(complete_rarefied[rownames(abundance_matrix), colnames(abundance_matrix),
  drop = FALSE])

stopifnot(identical(dimnames(abundance_matrix), dimnames(expected_subset)))

stopifnot(all(abundance_matrix == expected_subset))

hit_ids <- sort(unique(hits$qseqid))

abundance_ids <- sort(abundance_raw$ASV_ID)

if (!identical(hit_ids, abundance_ids)) {
  stop("MBPD query IDs and candidate abundance ASV IDs do not match. Missing taxonomy=", length(setdiff(abundance_ids,
    hit_ids)), "; extra taxonomy=", length(setdiff(hit_ids, abundance_ids)))
}

if (length(hit_ids) != 2755L) stop("Expected 2,755 candidate ASVs; found ", length(hit_ids))

prepare_level_data <- function(taxon_level) {
  tax_index <- match(rownames(abundance_matrix), taxonomy$ASV_ID)
  taxon <- taxonomy[[taxon_level]][tax_index]
  keep <- !is.na(taxon) & is_informative(taxon)
  if (!any(keep)) {
    return(list(data = NULL, audit = data.frame(Taxonomic_Level = taxon_level, Input_ASVs = nrow(abundance_matrix),
      Samples = ncol(abundance_matrix), Taxonomy_Matched_ASVs = sum(!is.na(tax_index)), Classified_ASVs = 0L,
      Classified_Taxa = 0L, Group_Bars = 0L, Max_Group_Sum_Error = NA_real_, stringsAsFactors = FALSE)))
  }
  taxon_sum <- rowsum(abundance_matrix[keep, , drop = FALSE], group = as.character(taxon[keep]), reorder = FALSE)
  sample_totals <- colSums(taxon_sum)
  valid_samples <- sample_totals > 0
  taxon_relative <- sweep(taxon_sum[, valid_samples, drop = FALSE], 2, sample_totals[valid_samples],
    "/") * 100
  plot_parts <- list()
  part_index <- 0L
  for (category in category_levels) {
    metadata <- group_metadata[[category]]
    contract <- group_contract[[category]]
    for (raw_group in contract$order) {
      target_samples <- metadata$SampleID[metadata$GroupRaw == raw_group]
      samples <- intersect(colnames(taxon_relative), target_samples)
      if (length(samples) == 0L)
        next
      group_mean <- rowMeans(taxon_relative[, samples, drop = FALSE])
      top_taxa <- names(head(sort(group_mean, decreasing = TRUE), top_n))
      display_taxon <- ifelse(names(group_mean) %in% top_taxa, names(group_mean), "Others")
      summary <- aggregate(group_mean, by = list(DisplayTaxon = display_taxon), FUN = sum)
      names(summary)[2] <- "Abundance"
      summary$Group <- unname(contract$labels[[raw_group]])
      summary$GroupRaw <- raw_group
      summary$Category <- category
      summary$SampleN <- length(samples)
      summary$Taxonomic_Level <- taxon_level
      part_index <- part_index + 1L
      plot_parts[[part_index]] <- summary
    }
  }
  plot_data <- bind_rows(plot_parts)
  plot_data$Category <- factor(plot_data$Category, levels = category_levels)
  plot_data$Group <- factor(plot_data$Group, levels = group_levels)
  sums <- aggregate(Abundance ~ Taxonomic_Level + Category + Group, data = plot_data, FUN = sum)
  audit <- data.frame(Taxonomic_Level = taxon_level, Input_ASVs = nrow(abundance_matrix), Samples = ncol(abundance_matrix),
    Taxonomy_Matched_ASVs = sum(!is.na(tax_index)), Classified_ASVs = sum(keep), Classified_Taxa = nrow(taxon_sum),
    Group_Bars = nrow(sums), Top_N_Per_Group = top_n, Max_Group_Sum_Error = max(abs(sums$Abundance -
      100)), stringsAsFactors = FALSE)
  list(data = plot_data, audit = audit)
}

base_palette <- c("#E18727FF", "#FFDC91FF", "#91D1C2FF", "#B09C85FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF",
  "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#7E6148FF", "#008B45FF", "#631879FF", "#A20056FF", "#5F559BFF",
  "#BC3C29FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#DC0000FF", "#1f77b4", "#ff7f0e", "#2ca02c",
  "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

make_colour_map <- function(plot_data) {
  scores <- tapply(plot_data$Abundance[plot_data$DisplayTaxon != "Others"], plot_data$DisplayTaxon[plot_data$DisplayTaxon !=
    "Others"], mean)
  ordered_taxa <- names(sort(scores, decreasing = TRUE))
  colours <- if (length(ordered_taxa) <= length(base_palette)) {
    base_palette[seq_along(ordered_taxa)]
  }
  else {
    (grDevices::colorRampPalette(base_palette))(length(ordered_taxa))
  }
  setNames(c(colours, "#E0E0E0"), c(ordered_taxa, "Others"))
}

title_case <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

make_plot <- function(plot_data, taxon_level, colour_map) {
  breaks <- c(names(colour_map)[names(colour_map) != "Others"], "Others")
  plot_data$DisplayTaxon <- factor(plot_data$DisplayTaxon, levels = names(colour_map))
  italic_legend <- taxon_level %in% c("genus", "species")
  legend_columns <- if (length(breaks) <= 8L)
    2L
  else 3L
  ggplot(plot_data, aes(x = Group, y = Abundance, fill = DisplayTaxon)) + geom_col(width = 0.56, colour = "black",
    linewidth = 0.2, position = position_stack(reverse = TRUE)) + facet_grid(~Category, scales = "free_x",
    space = "free_x", switch = "x", drop = TRUE) + scale_fill_manual(values = colour_map, breaks = breaks,
    drop = FALSE, guide = guide_legend(title = title_case(taxon_level), ncol = legend_columns, byrow = TRUE,
      override.aes = list(colour = NA))) + scale_y_continuous(limits = c(0, 100.1), breaks = c(0,
    25, 50, 75, 100), labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0,
    0))) + labs(x = NULL, y = "Mean relative abundance", title = paste0("Pathogen-candidate ASVs (MBPD top1) · ",
    title_case(taxon_level), " composition"), subtitle = paste0("Top 15 taxa within each group; Others pooled; classified subset renormalized to 100%")) +
    theme_classic(base_size = 9.2, base_family = "Arial") + theme(plot.background = element_rect(fill = "white",
    colour = NA), panel.background = element_rect(fill = "white", colour = NA), panel.grid.major.y = element_line(colour = "#E8E8E8",
    linewidth = 0.28), panel.grid.minor = element_blank(), axis.line.x = element_blank(), axis.line.y = element_line(colour = "#252525",
    linewidth = 0.42), axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour = "#252525",
    linewidth = 0.35), axis.text.y = element_text(colour = "#252525", size = 8.5, margin = margin(r = 3)),
    axis.text.x = element_text(colour = "#252525", size = 8.1, face = "bold", angle = 48, hjust = 1,
      vjust = 1), axis.title.y = element_text(colour = "#202020", size = 9.3, face = "bold", margin = margin(r = 8)),
    strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(colour = "#202020",
      size = 9.3, face = "bold", margin = margin(t = 6)), panel.spacing.x = grid::unit(4.5, "mm"),
    plot.title = element_text(colour = "#171717", size = 11.2, face = "bold"), plot.subtitle = element_text(colour = "#5A5A5A",
      size = 8.4, margin = margin(b = 9)), legend.position = "bottom", legend.direction = "horizontal",
    legend.title.position = "top", legend.title = element_text(size = 8.4, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 7.8, face = if (italic_legend)
      "italic"
    else "plain", colour = "#303030"), legend.key.width = grid::unit(3.6, "mm"), legend.key.height = grid::unit(3.2,
      "mm"), legend.spacing.x = grid::unit(2.5, "mm"), legend.box = "vertical", legend.box.just = "center",
    legend.box.margin = margin(t = 3), plot.margin = margin(t = 8, r = 8, b = 7, l = 8))
}

make_size_plot <- function(plot_data, taxon_level, colour_map) {
  present_taxa <- unique(as.character(plot_data$DisplayTaxon))
  colour_map <- colour_map[names(colour_map) %in% present_taxa]
  breaks <- c(names(colour_map)[names(colour_map) != "Others"], "Others")
  plot_data$DisplayTaxon <- factor(plot_data$DisplayTaxon, levels = names(colour_map))
  ggplot(plot_data, aes(x = Group, y = Abundance, fill = DisplayTaxon)) + geom_col(width = 0.52, colour = "black",
    linewidth = 0.2, position = position_stack(reverse = TRUE)) + scale_x_discrete(expand = expansion(add = 0.75)) +
    scale_fill_manual(values = colour_map, breaks = breaks, drop = FALSE, guide = guide_legend(title = title_case(taxon_level),
      ncol = 2L, byrow = FALSE, override.aes = list(colour = NA))) + scale_y_continuous(limits = c(0,
    100.1), breaks = c(0, 25, 50, 75, 100), labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0,
    0))) + labs(x = "Size", y = "Mean relative abundance", title = paste0("Pathogen-candidate ASVs (MBPD top1) · ",
    title_case(taxon_level), " composition"), subtitle = paste0("Size groups only; Top 15 taxa within each size group; ",
    "Others pooled; classified subset renormalized to 100%")) + theme_classic(base_size = 9.2, base_family = "Arial") +
    theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white",
      colour = NA), panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.28), panel.grid.minor = element_blank(),
      axis.line.x = element_blank(), axis.line.y = element_line(colour = "#252525", linewidth = 0.42),
      axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour = "#252525", linewidth = 0.35),
      axis.text.y = element_text(colour = "#252525", size = 8.5, margin = margin(r = 3)), axis.text.x = element_text(colour = "#252525",
        size = 8.5, face = "bold"), axis.title.x = element_text(colour = "#202020", size = 9.3,
        face = "bold", margin = margin(t = 7)), axis.title.y = element_text(colour = "#202020",
        size = 9.3, face = "bold", margin = margin(r = 8)), plot.title = element_text(colour = "#171717",
        size = 11.2, face = "bold"), plot.subtitle = element_text(colour = "#5A5A5A", size = 8.4,
        margin = margin(b = 9)), legend.position = "right", legend.direction = "vertical", legend.title.position = "top",
      legend.title = element_text(size = 8.4, face = "bold", hjust = 0.5), legend.text = element_text(size = 7.8,
        face = "italic", colour = "#303030"), legend.key.width = grid::unit(3.6, "mm"), legend.key.height = grid::unit(3.2,
        "mm"), legend.spacing.x = grid::unit(2.2, "mm"), legend.spacing.y = grid::unit(0.5, "mm"),
      legend.box = "vertical", legend.box.just = "top", legend.box.margin = margin(l = 5), plot.margin = margin(t = 8,
        r = 8, b = 7, l = 8))
}

figure_height_mm <- function(n_legend_taxa) {
  legend_columns <- if (n_legend_taxa <= 8L)
    2L
  else 3L
  legend_rows <- ceiling(n_legend_taxa/legend_columns)
  min(215, max(145, 121 + legend_rows * 5.4))
}

save_plot <- function(plot, file_base, height_mm, key_level = FALSE) {
  grDevices::png(paste0(file_base, ".png"), width = 250, height = height_mm, units = "mm", res = 400,
    type = if (capabilities("cairo"))
      "cairo-png"
    else "windows", bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(file_base, ".pdf"), width = 250/25.4, height = height_mm/25.4, family = "Arial",
    onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  if (key_level) {
    grDevices::svg(paste0(file_base, ".svg"), width = 250/25.4, height = height_mm/25.4, family = "sans",
      bg = "white")
    print(plot)
    grDevices::dev.off()
    grDevices::tiff(paste0(file_base, ".tiff"), width = 250, height = height_mm, units = "mm", res = 600,
      compression = "lzw", type = if (capabilities("cairo"))
        "cairo"
      else "windows", bg = "white")
    print(plot)
    grDevices::dev.off()
  }
}

results <- list()

audits <- list()

plots <- list()

colour_tables <- list()

for (index in seq_along(taxonomic_levels)) {
  taxon_level <- taxonomic_levels[index]
  message("Preparing MBPD ", taxon_level, " composition...")
  result <- prepare_level_data(taxon_level)
  results[[taxon_level]] <- result$data
  audits[[taxon_level]] <- result$audit
  if (is.null(result$data) || nrow(result$data) == 0L)
    next
  colour_map <- make_colour_map(result$data)
  colour_tables[[taxon_level]] <- data.frame(Taxonomic_Level = taxon_level, Taxon = names(colour_map),
    Colour = unname(colour_map), stringsAsFactors = FALSE)
  figure_data <- result$data
  if (identical(taxon_level, "genus")) {
    original_top_n <- top_n
    top_n <- 1000000L
    complete_genus <- prepare_level_data(taxon_level)$data
    top_n <- original_top_n
    write.csv(complete_genus[complete_genus$Category == "Size", , drop = FALSE], file.path(table_dir,
      "06_genus_size_complete_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    figure_data <- result$data[result$data$Category == "Size", , drop = FALSE]
    plot <- make_size_plot(figure_data, taxon_level, colour_map)
  }
  else {
    plot <- make_plot(figure_data, taxon_level, colour_map)
  }
  plots[[taxon_level]] <- plot
  height_mm <- if (identical(taxon_level, "genus")) {
    145
  }
  else {
    figure_height_mm(length(unique(figure_data$DisplayTaxon)))
  }
  file_base <- file.path(plot_dir, paste0(sprintf("%02d", index), "_", taxon_level))
  save_plot(plot, file_base, height_mm, key_level = taxon_level %in% c("phylum", "genus"))
  write.csv(result$data, file.path(table_dir, paste0(sprintf("%02d", index), "_", taxon_level, "_plot_data.csv")),
    row.names = FALSE, fileEncoding = "UTF-8")
  if (identical(taxon_level, "genus")) {
    write.csv(figure_data, file.path(table_dir, "06_genus_size_plot_data.csv"), row.names = FALSE,
      fileEncoding = "UTF-8")
  }
}

audit <- bind_rows(audits)

if (any(audit$Taxonomy_Matched_ASVs != 2755L)) stop("Taxonomy matching audit failed.")

if (any(audit$Group_Bars != 12L)) stop("Expected 12 group bars per taxonomic level.")

if (any(audit$Max_Group_Sum_Error > 1e-08, na.rm = TRUE)) stop("Stacked bars do not sum to 100%.")

hit_audit <- data.frame(Metric = c("Raw_hit_rows", "Unique_query_ASVs", "Candidate_ASVs_matched", "Minimum_pident",
  "Minimum_qcovs", "ASVs_with_top_bitscore_ties", "ASVs_with_conflicting_top_bitscore_taxonomies",
  "ASVs_with_final_numeric_ties", "ASVs_with_conflicting_final_numeric_taxonomies", "ASVs_with_single_selected_annotation"),
  Value = c(nrow(hits), n_distinct(hits$qseqid), length(intersect(hit_ids, abundance_ids)), min(hits$pident),
    min(hits$qcovs), sum(taxonomy$Top_Bitscore_Hit_Rows > 1L), sum(taxonomy$Top_Bitscore_Distinct_Taxonomies >
      1L), sum(taxonomy$Final_Numeric_Tie_Rows > 1L), sum(taxonomy$Final_Numeric_Tie_Distinct_Taxonomies >
      1L), sum(taxonomy$Best_Hit_Rows == 1L)), stringsAsFactors = FALSE)

write.csv(taxonomy, file.path(table_dir, "MBPD_ASV_taxonomy_LCA_2755.csv"), row.names = FALSE)

write.csv(audit, file.path(table_dir, "MBPD_composition_audit.csv"), row.names = FALSE)

write.csv(hit_audit, file.path(table_dir, "MBPD_hit_selection_audit.csv"), row.names = FALSE)

write.csv(bind_rows(colour_tables), file.path(table_dir, "MBPD_colour_map.csv"), row.names = FALSE)

write.csv(bind_rows(group_metadata), file.path(table_dir, "MBPD_group_metadata_audit.csv"), row.names = FALSE)

if (length(plots) > 0L) {
  grDevices::cairo_pdf(file.path(plot_dir, "MBPD_pathogen_candidate_all_taxonomic_levels.pdf"), width = 250/25.4,
    height = 170/25.4, family = "Arial", onefile = TRUE)
  for (plot in plots) print(plot)
  grDevices::dev.off()
}

readme <- c("MBPD注释的候选致病菌七层级分组组成图", "", "输入：02_filtered_hits_with_taxonomy.tsv；2,755-ASV候选丰度表；四套分组文件。",
  "命中表57,964行，覆盖全部2,755个候选ASV；最低pident=97%，最低qcovs=95%。", "",
  "注释规则：", "1. 每个ASV只选择一个命中：最高bitscore；并列时依次取最高pident、最高qcovs、最长alignment、最低E-value。",
  "2. 数值指标仍完全并列时按sseqid字典序选择一条，保证结果可重复；不再汇总全部命中或做LCA。",
  "3. 所选记录中的uncultured、unclassified、unassigned、unknown、metagenome等不作为有效分类名。",
  "", "组成规则：先在每个样本内将该层级已分类ASV换算为相对丰度，再按Size、Location、Season、Rainfall分组求均值；每组Top15，其余为Others。",
  "注意：每个层级只对该层级已分类部分重新归一化到100%，不同层级的已分类ASV数见MBPD_composition_audit.csv。",
  "SILVA与Greengenes结果保留，本结果单独存放，不覆盖。")

readme <- c("MBPD注释的候选致病菌七层级分组组成图", "", "输入：02_filtered_hits_with_taxonomy.tsv、2,755-ASV候选丰度表、四套分组文件。",
  "命中表57,964行，覆盖全部2,755个候选ASV；最低pident=97%，最低qcovs=95%。", "",
  "注释规则：", "1. 每个ASV只选择一个命中：最高bitscore；并列时依次取最高pident、最高qcovs、最长alignment、最低E-value。",
  "2. 数值指标仍完全并列时按sseqid字典序选择一条，保证结果可重复；不再汇总全部命中或做LCA。",
  "3. 所选记录中的uncultured、unclassified、unassigned、unknown、metagenome等不作为有效分类名。",
  "", "组成规则：先在每个样本内将该层级已分类ASV换算为相对丰度，再按Size、Location、Season、Rainfall分组求均值；每组Top15，其余合并为Others。",
  "展示调整：属水平图仅展示Size组（Large和Small），并仅列出这两个柱中实际出现的属；完整四类分组数据仍保留在06_genus_plot_data.csv。",
  "注意：每个层级只对该层级已分类部分重新归一化到100%；不同层级的已分类ASV数见MBPD_composition_audit.csv。",
  "SILVA与Greengenes结果保留；本结果单独存放，不覆盖。")

writeLines(readme, file.path(table_dir, "README_CN.txt"), useBytes = TRUE)

writeLines(c("Input counts: candidate ASV rows extracted from the complete community rarefied to 21,000 reads per sample.",
  "The candidate subset is NOT rarefied again. Exact count equality with the complete-community subset is checked before analysis.",
  "Taxonomy, sample-wise classified-subset normalization and group averaging are unchanged.", paste0("Candidate input: ",
    abundance_file), paste0("Complete input: ", complete_rarefied_file)), file.path(table_dir, "rarefied_input_provenance.txt"))

capture.output(sessionInfo(), file = file.path(table_dir, "sessionInfo.txt"))

cat("\nMBPD pathogen composition completed.\n")

cat("Output directory: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")

print(hit_audit)

print(audit)
