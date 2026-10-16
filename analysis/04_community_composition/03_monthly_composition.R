local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
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

  analysis_root <- file.path(project_root, "01_分析项目", "10_物种组成年际变化")

  rarefied_file <- file.path(project_root, "00_输入文件", "01_ASV与序列", "ASV抽平21000_36097ASV_120样本.csv")

  month_file <- file.path(project_root, "00_输入文件", "02_样本分组与站点", "月份分组_groupmonth.csv")

  size_file <- file.path(project_root, "00_输入文件", "02_样本分组与站点", "孔径分组_groupsize.csv")

  candidate_file <- file.path(project_root, "00_输入文件", "05_致病菌候选", "候选致病菌ASV_ID清单_2755.csv")

  silva_tax_file <- file.path(project_root, "00_输入文件", "08_派生输入", "分类注释标准化",
    "SILVA", "tax.csv")

  greengenes_tax_file <- file.path(project_root, "00_输入文件", "08_派生输入", "分类注释标准化",
    "Greengenes", "tax.csv")

  mbpd_tax_file <- file.path(project_root, "00_输入文件", "04_分类注释", "MBPD", "MBPD_ASV_taxonomy_LCA_2755.csv")

  required_files <- c(rarefied_file, month_file, size_file, candidate_file, silva_tax_file, greengenes_tax_file,
    mbpd_tax_file)

  missing_files <- required_files[!file.exists(required_files)]

  if (length(missing_files) > 0L) {
    stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
  }

  top_n <- 15L

  taxonomic_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  level_labels <- c(domain = "Domain", phylum = "Phylum", class = "Class", order = "Order", family = "Family",
    genus = "Genus", species = "Species")

  panel_levels <- c("Overall", "ConF", "USF")

  read_abundance_matrix <- function(path) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    names(x)[1] <- "ASV_ID"
    x$ASV_ID <- gsub("-", "_", trimws(as.character(x$ASV_ID)))
    if (anyDuplicated(x$ASV_ID))
      stop("Duplicated ASV IDs in rarefied table.")
    m <- matrix(as.numeric(as.matrix(x[, -1, drop = FALSE])), nrow = nrow(x), ncol = ncol(x) - 1L, dimnames = list(x$ASV_ID,
      names(x)[-1]))
    if (anyNA(m) || any(m < 0)) {
      stop("Rarefied abundance table contains NA, non-numeric, or negative values.")
    }
    m
  }

  read_two_column_group <- function(path, value_name) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    if (ncol(x) < 2L)
      stop("Invalid group file: ", path)
    names(x)[1:2] <- c("SampleID", value_name)
    x <- x[, c("SampleID", value_name), drop = FALSE]
    x$SampleID <- trimws(as.character(x$SampleID))
    x[[value_name]] <- trimws(as.character(x[[value_name]]))
    x
  }

  read_taxonomy <- function(path) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    id_column <- intersect(c("ASV_ID", "mags"), names(x))[1]
    if (is.na(id_column))
      stop("Taxonomy table lacks ASV_ID/mags: ", path)
    names(x)[names(x) == id_column] <- "ASV_ID"
    x$ASV_ID <- gsub("-", "_", trimws(as.character(x$ASV_ID)))
    for (level in intersect(taxonomic_levels, names(x))) {
      x[[level]] <- trimws(as.character(x[[level]]))
    }
    if (anyDuplicated(x$ASV_ID))
      stop("Duplicated taxonomy ASV IDs: ", path)
    x
  }

  read_candidate_ids <- function(path) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    ids <- gsub("-", "_", trimws(as.character(x[[1]])))
    ids <- unique(ids[nzchar(ids)])
    ids
  }

  abundance_all <- read_abundance_matrix(rarefied_file)

  month_meta <- read_two_column_group(month_file, "MonthCode")

  size_meta <- read_two_column_group(size_file, "SizeRaw")

  candidate_ids <- read_candidate_ids(candidate_file)

  if (nrow(abundance_all) != 36097L || ncol(abundance_all) != 120L) {
    stop("Expected 36,097 ASVs x 120 samples in the rarefied table; found ", nrow(abundance_all), " x ",
      ncol(abundance_all), ".")
  }

  sample_totals <- colSums(abundance_all)

  if (any(sample_totals != 21000)) {
    stop("Rarefaction audit failed: not every sample sums to 21,000 reads.")
  }

  if (!setequal(colnames(abundance_all), month_meta$SampleID)) {
    stop("Month metadata sample IDs do not match the rarefied table.")
  }

  if (!setequal(colnames(abundance_all), size_meta$SampleID)) {
    stop("Size metadata sample IDs do not match the rarefied table.")
  }

  if (length(candidate_ids) != 2755L) {
    stop("Expected 2,755 pathogen-candidate ASV IDs; found ", length(candidate_ids), ".")
  }

  missing_candidates <- setdiff(candidate_ids, rownames(abundance_all))

  if (length(missing_candidates) > 0L) {
    stop("Candidate ASVs missing from rarefied table: ", length(missing_candidates))
  }

  size_value <- tolower(size_meta$SizeRaw)

  if (!setequal(unique(size_value), c("large", "small"))) {
    stop("Unexpected size labels: ", paste(sort(unique(size_value)), collapse = ", "))
  }

  size_meta$SizePanel <- ifelse(size_value == "large", "ConF", "USF")

  month_meta$MonthDate <- as.Date(paste0("20", month_meta$MonthCode, "01"), format = "%Y%m%d")

  if (anyNA(month_meta$MonthDate)) stop("Failed to parse month codes.")

  observed_months <- sort(unique(month_meta$MonthDate))

  if (length(observed_months) != 15L) {
    stop("Expected 15 sampled months; found ", length(observed_months), ".")
  }

  sample_metadata <- month_meta %>%
    left_join(size_meta[, c("SampleID", "SizePanel")], by = "SampleID")

  if (anyNA(sample_metadata$SizePanel)) stop("Missing size metadata after joining.")

  silva_tax <- read_taxonomy(silva_tax_file)

  greengenes_tax <- read_taxonomy(greengenes_tax_file)

  mbpd_tax <- read_taxonomy(mbpd_tax_file)

  dataset_definitions <- list(SILVA_all = list(database = "SILVA", dataset_key = "all", dataset_label = "All ASVs",
    taxonomy = silva_tax, image_dir = file.path(analysis_root, "结果图片", "SILVA", "整体_物种组成年际变化_抽平逐月版"),
    table_dir = file.path(analysis_root, "结果表格", "SILVA", "整体_物种组成年际变化_抽平逐月版")),
    SILVA_pathogen = list(database = "SILVA", dataset_key = "pathogen", dataset_label = "Pathogen-candidate ASVs",
      taxonomy = silva_tax, image_dir = file.path(analysis_root, "结果图片", "SILVA", "致病菌_物种组成年际变化_抽平逐月版"),
      table_dir = file.path(analysis_root, "结果表格", "SILVA", "致病菌_物种组成年际变化_抽平逐月版")),
    Greengenes_all = list(database = "Greengenes", dataset_key = "all", dataset_label = "All ASVs", taxonomy = greengenes_tax,
      image_dir = file.path(analysis_root, "结果图片", "Greengenes", "整体_物种组成年际变化_抽平逐月版"),
      table_dir = file.path(analysis_root, "结果表格", "Greengenes", "整体_物种组成年际变化_抽平逐月版")),
    Greengenes_pathogen = list(database = "Greengenes", dataset_key = "pathogen", dataset_label = "Pathogen-candidate ASVs",
      taxonomy = greengenes_tax, image_dir = file.path(analysis_root, "结果图片", "Greengenes",
        "致病菌_物种组成年际变化_抽平逐月版"), table_dir = file.path(analysis_root,
        "结果表格", "Greengenes", "致病菌_物种组成年际变化_抽平逐月版")), MBPD_pathogen = list(database = "MBPD top1",
      dataset_key = "pathogen", dataset_label = "Pathogen-candidate ASVs", taxonomy = mbpd_tax, image_dir = file.path(analysis_root,
        "结果图片", "MBPD", "致病菌_物种组成年际变化_抽平逐月版"), table_dir = file.path(analysis_root,
        "结果表格", "MBPD", "致病菌_物种组成年际变化_抽平逐月版")))

  for (definition in dataset_definitions) {
    dir.create(definition$image_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(definition$table_dir, recursive = TRUE, showWarnings = FALSE)
  }

  is_informative <- function(x) {
    x <- trimws(as.character(x))
    !is.na(x) & nzchar(x) & !grepl(paste0("uncultured|unclassified|unassigned|unknown|metagenome|", "environmental sample|incertae[ _-]*sedis|^bacterium$|^na$"),
      x, ignore.case = TRUE)
  }

  get_dataset_matrix <- function(dataset_key) {
    if (dataset_key == "all")
      return(abundance_all)
    abundance_all[candidate_ids, , drop = FALSE]
  }

  prepare_monthly_level <- function(definition, taxon_level) {
    abundance <- get_dataset_matrix(definition$dataset_key)
    taxonomy <- definition$taxonomy
    tax_index <- match(rownames(abundance), taxonomy$ASV_ID)
    taxon <- taxonomy[[taxon_level]][tax_index]
    keep <- !is.na(tax_index) & is_informative(taxon)
    empty_audit <- data.frame(Database = definition$database, Dataset = definition$dataset_label, Taxonomic_Level = taxon_level,
      Input_ASVs = nrow(abundance), Taxonomy_Matched_ASVs = sum(!is.na(tax_index)), Classified_ASVs = sum(keep),
      Classified_Taxa = 0L, Samples = ncol(abundance), Sampled_Months = length(observed_months), Panel_Month_Bars = 0L,
      Top_N = top_n, Max_Bar_Sum_Error = NA_real_, Status = "skipped_no_informative_taxonomy", stringsAsFactors = FALSE)
    if (!any(keep))
      return(list(data = NULL, audit = empty_audit))
    taxon_sum <- rowsum(abundance[keep, , drop = FALSE], group = as.character(taxon[keep]), reorder = FALSE)
    classified_totals <- colSums(taxon_sum)
    valid_samples <- classified_totals > 0
    if (!any(valid_samples))
      return(list(data = NULL, audit = empty_audit))
    relative <- sweep(taxon_sum[, valid_samples, drop = FALSE], 2, classified_totals[valid_samples],
      "/") * 100
    global_mean <- rowMeans(relative)
    top_taxa <- names(head(sort(global_mean, decreasing = TRUE), top_n))
    plot_parts <- list()
    part_index <- 0L
    for (panel in panel_levels) {
      panel_sample_ids <- if (panel == "Overall") {
        sample_metadata$SampleID
      }
      else {
        sample_metadata$SampleID[sample_metadata$SizePanel == panel]
      }
      for (month_date in observed_months) {
        month_date <- as.Date(month_date, origin = "1970-01-01")
        month_sample_ids <- sample_metadata$SampleID[sample_metadata$MonthDate == month_date]
        samples <- intersect(colnames(relative), intersect(panel_sample_ids, month_sample_ids))
        if (length(samples) == 0L)
          next
        month_mean <- rowMeans(relative[, samples, drop = FALSE])
        display_taxon <- ifelse(names(month_mean) %in% top_taxa, names(month_mean), "Others")
        summary <- aggregate(month_mean, by = list(DisplayTaxon = display_taxon), FUN = sum)
        names(summary)[2] <- "Abundance"
        summary$Panel <- panel
        summary$MonthDate <- month_date
        summary$Month <- format(month_date, "%Y-%m")
        summary$SampleN <- length(samples)
        summary$Taxonomic_Level <- taxon_level
        summary$Database <- definition$database
        summary$Dataset <- definition$dataset_label
        part_index <- part_index + 1L
        plot_parts[[part_index]] <- summary
      }
    }
    plot_data <- bind_rows(plot_parts)
    plot_data$Panel <- factor(plot_data$Panel, levels = panel_levels)
    plot_data$Month <- factor(plot_data$Month, levels = format(observed_months, "%Y-%m"))
    sums <- aggregate(Abundance ~ Panel + Month, data = plot_data, FUN = sum)
    audit <- empty_audit
    audit$Classified_Taxa <- nrow(taxon_sum)
    audit$Panel_Month_Bars <- nrow(sums)
    audit$Max_Bar_Sum_Error <- max(abs(sums$Abundance - 100))
    audit$Status <- "success"
    if (audit$Panel_Month_Bars != length(panel_levels) * length(observed_months)) {
      stop("Unexpected number of panel-month bars for ", definition$database, " ", definition$dataset_label,
        " ", taxon_level, ".")
    }
    if (audit$Max_Bar_Sum_Error > 1e-08) {
      stop("Monthly stacked bars do not sum to 100%.")
    }
    list(data = plot_data, audit = audit)
  }

  base_palette <- c("#E18727FF", "#FFDC91FF", "#91D1C2FF", "#B09C85FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF",
    "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#7E6148FF", "#008B45FF", "#631879FF", "#A20056FF", "#5F559BFF",
    "#BC3C29FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#DC0000FF")

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

  prepare_panel_geometry <- function(plot_data, panel, colour_map, bar_height = 0.58) {
    month_levels <- levels(plot_data$Month)
    taxon_levels <- names(colour_map)
    panel_abundance <- plot_data %>%
      filter(Panel == panel) %>%
      transmute(Month = as.character(Month), DisplayTaxon = as.character(DisplayTaxon), Abundance = Abundance)
    bar_data <- expand.grid(Month = month_levels, DisplayTaxon = taxon_levels, stringsAsFactors = FALSE) %>%
      left_join(panel_abundance, by = c("Month", "DisplayTaxon"))
    bar_data$Abundance[is.na(bar_data$Abundance)] <- 0
    bar_data$Month <- factor(bar_data$Month, levels = month_levels, ordered = TRUE)
    bar_data$DisplayTaxon <- factor(bar_data$DisplayTaxon, levels = taxon_levels, ordered = TRUE)
    bar_data <- bar_data %>%
      arrange(Month, DisplayTaxon) %>%
      group_by(Month) %>%
      mutate(MonthIndex = as.integer(Month), StackMin = cumsum(Abundance) - Abundance, StackMax = cumsum(Abundance),
        RowMin = MonthIndex - bar_height/2, RowMax = MonthIndex + bar_height/2) %>%
      ungroup()
    connector_parts <- list()
    connector_index <- 0L
    month_indices <- sort(unique(bar_data$MonthIndex))
    if (length(month_indices) >= 2L) {
      for (month_position in seq_len(length(month_indices) - 1L)) {
        upper_index <- month_indices[[month_position]]
        lower_index <- month_indices[[month_position + 1L]]
        upper_bar <- bar_data %>%
          filter(MonthIndex == upper_index) %>%
          transmute(DisplayTaxon, UpperY = RowMax, UpperXMin = StackMin, UpperXMax = StackMax,
            UpperAbundance = Abundance)
        lower_bar <- bar_data %>%
          filter(MonthIndex == lower_index) %>%
          transmute(DisplayTaxon, LowerY = RowMin, LowerXMin = StackMin, LowerXMax = StackMax,
            LowerAbundance = Abundance)
        joined <- inner_join(upper_bar, lower_bar, by = "DisplayTaxon")
        joined <- joined[joined$UpperAbundance > 0 | joined$LowerAbundance > 0, , drop = FALSE]
        if (nrow(joined) == 0L)
          next
        for (row_index in seq_len(nrow(joined))) {
          connector_index <- connector_index + 1L
          row <- joined[row_index, , drop = FALSE]
          connector_parts[[connector_index]] <- data.frame(ConnectorID = paste0(upper_index, "_",
            lower_index, "_", as.character(row$DisplayTaxon)), DisplayTaxon = row$DisplayTaxon,
            X = c(row$UpperXMin, row$LowerXMin, row$LowerXMax, row$UpperXMax), Y = c(row$UpperY,
            row$LowerY, row$LowerY, row$UpperY), VertexOrder = 1:4, stringsAsFactors = FALSE)
        }
      }
    }
    connector_data <- if (length(connector_parts) == 0L) {
      data.frame(ConnectorID = character(), DisplayTaxon = character(), X = numeric(), Y = numeric(),
        VertexOrder = integer(), stringsAsFactors = FALSE)
    }
    else {
      bind_rows(connector_parts)
    }
    connector_data$DisplayTaxon <- factor(connector_data$DisplayTaxon, levels = taxon_levels, ordered = TRUE)
    list(bars = bar_data, connectors = connector_data, month_levels = month_levels)
  }

  make_plot <- function(plot_data, panel, definition, taxon_level, colour_map) {
    breaks <- c(names(colour_map)[names(colour_map) != "Others"], "Others")
    italic_legend <- taxon_level %in% c("genus", "species")
    legend_columns <- 1L
    geometry <- prepare_panel_geometry(plot_data = plot_data, panel = panel, colour_map = colour_map)
    ggplot() + geom_polygon(data = geometry$connectors, aes(x = X, y = Y, group = ConnectorID, fill = DisplayTaxon),
      colour = NA, alpha = 0.7) + geom_rect(data = geometry$bars, aes(xmin = StackMin, xmax = StackMax,
      ymin = RowMin, ymax = RowMax, fill = DisplayTaxon), colour = "black", linewidth = 0.18) + scale_fill_manual(values = colour_map,
      breaks = breaks, drop = FALSE, guide = guide_legend(title = level_labels[[taxon_level]], ncol = legend_columns,
        byrow = TRUE, override.aes = list(colour = NA))) + scale_x_continuous(limits = c(0, 100.1),
      breaks = c(0, 25, 50, 75, 100), labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0,
        0))) + scale_y_reverse(breaks = seq_along(geometry$month_levels), labels = geometry$month_levels,
      limits = c(length(geometry$month_levels) + 0.55, 0.45), expand = expansion(mult = c(0, 0))) +
      labs(x = "Mean relative abundance", y = "Month", title = paste0(panel, " · ", definition$dataset_label,
        " (", definition$database, ") · Monthly ", tolower(level_labels[[taxon_level]]), " composition"),
        subtitle = paste0("21,000-read rarefied ASV table; one horizontal bar per sampled month; ",
          "global Top 15 taxa; Others pooled; incertae sedis excluded;\n", "connecting bands track the same taxon between adjacent months")) +
      theme_classic(base_size = 11.8, base_family = "Arial") + theme(plot.background = element_rect(fill = "white",
      colour = NA), panel.background = element_rect(fill = "white", colour = NA), panel.grid.major.x = element_line(colour = "#E8E8E8",
      linewidth = 0.28), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = "#252525", linewidth = 0.42), axis.line.y = element_line(colour = "#252525",
        linewidth = 0.42), axis.ticks.x = element_line(colour = "#252525", linewidth = 0.3), axis.ticks.y = element_line(colour = "#252525",
        linewidth = 0.35), axis.text.y = element_text(colour = "#252525", size = 10.5, face = "bold",
        margin = margin(r = 4)), axis.text.x = element_text(colour = "#252525", size = 10.5), axis.title = element_text(colour = "#202020",
        size = 11.8, face = "bold"), axis.title.y = element_text(margin = margin(r = 8)), axis.title.x = element_text(margin = margin(t = 8)),
      plot.title = element_text(colour = "#171717", size = 14, face = "bold"), plot.subtitle = element_text(colour = "#5A5A5A",
        size = 10.2, margin = margin(b = 9)), legend.position = "bottom", legend.direction = "vertical",
      legend.title.position = "top", legend.title = element_text(size = 10.2, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 9.5, face = if (italic_legend)
        "italic"
      else "plain", colour = "#303030"), legend.key.width = grid::unit(3.6, "mm"), legend.key.height = grid::unit(3.2,
        "mm"), legend.spacing.x = grid::unit(2.5, "mm"), legend.box = "vertical", legend.box.just = "center",
      legend.box.margin = margin(t = 3), plot.margin = margin(t = 8, r = 8, b = 7, l = 8))
  }

  make_pyramid_plot <- function(plot_data, definition, taxon_level, colour_map) {
    breaks <- c(names(colour_map)[names(colour_map) != "Others"], "Others")
    italic_legend <- taxon_level %in% c("genus", "species")
    legend_columns <- 1L
    center_half_gap <- 10
    small_geometry <- prepare_panel_geometry(plot_data = plot_data, panel = "USF", colour_map = colour_map)
    large_geometry <- prepare_panel_geometry(plot_data = plot_data, panel = "ConF", colour_map = colour_map)
    small_bars <- small_geometry$bars %>%
      transmute(Side = "USF", Month, MonthIndex, DisplayTaxon, Abundance, PlotXMin = -center_half_gap -
        StackMax, PlotXMax = -center_half_gap - StackMin, RowMin, RowMax)
    large_bars <- large_geometry$bars %>%
      transmute(Side = "ConF", Month, MonthIndex, DisplayTaxon, Abundance, PlotXMin = center_half_gap +
        StackMin, PlotXMax = center_half_gap + StackMax, RowMin, RowMax)
    pyramid_bars <- bind_rows(small_bars, large_bars)
    small_connectors <- small_geometry$connectors %>%
      transmute(ConnectorID = paste0("Small_", ConnectorID), DisplayTaxon, X = -center_half_gap - X,
        Y, VertexOrder)
    large_connectors <- large_geometry$connectors %>%
      transmute(ConnectorID = paste0("Large_", ConnectorID), DisplayTaxon, X = center_half_gap + X,
        Y, VertexOrder)
    pyramid_connectors <- bind_rows(small_connectors, large_connectors)
    month_levels <- small_geometry$month_levels
    month_labels <- data.frame(MonthIndex = seq_along(month_levels), Month = month_levels, stringsAsFactors = FALSE)
    header_labels <- data.frame(X = c(-center_half_gap - 50, 0, center_half_gap + 50), Y = c(0.3, 0.3,
      0.3), Label = c("USF", "Month", "ConF"), stringsAsFactors = FALSE)
    pyramid_breaks <- c(-center_half_gap - c(100, 75, 50, 25, 0), center_half_gap + c(0, 25, 50, 75,
      100))
    pyramid_labels <- c("100%", "75%", "50%", "25%", "0%", "0%", "25%", "50%", "75%", "100%")
    ggplot() + geom_polygon(data = pyramid_connectors, aes(x = X, y = Y, group = ConnectorID, fill = DisplayTaxon),
      colour = NA, alpha = 0.7) + geom_vline(xintercept = c(-center_half_gap, center_half_gap), colour = "#252525",
      linewidth = 0.45) + geom_rect(data = pyramid_bars, aes(xmin = PlotXMin, xmax = PlotXMax, ymin = RowMin,
      ymax = RowMax, fill = DisplayTaxon), colour = "black", linewidth = 0.18) + geom_label(data = month_labels,
      aes(x = 0, y = MonthIndex, label = Month), inherit.aes = FALSE, fill = "white", colour = "#202020",
      fontface = "bold", size = 3.6, label.padding = grid::unit(0.65, "mm"), label.r = grid::unit(0,
        "mm"), linewidth = 0) + geom_text(data = header_labels, aes(x = X, y = Y, label = Label),
      inherit.aes = FALSE, colour = "#202020", fontface = "bold", size = 3.35) + scale_fill_manual(values = colour_map,
      breaks = breaks, drop = FALSE, guide = guide_legend(title = level_labels[[taxon_level]], ncol = legend_columns,
        byrow = TRUE, override.aes = list(colour = NA))) + scale_x_continuous(limits = c(-center_half_gap -
      103, center_half_gap + 103), breaks = pyramid_breaks, labels = pyramid_labels, expand = expansion(mult = c(0,
      0))) + scale_y_reverse(limits = c(length(month_levels) + 0.55, 0.05), expand = expansion(mult = c(0,
      0))) + labs(x = "Mean relative abundance", y = NULL, title = paste0(definition$dataset_label,
      " (", definition$database, ")\n", "USF (left) vs ConF (right) · Monthly ", tolower(level_labels[[taxon_level]]),
      " composition"), subtitle = paste0("21,000-read rarefied ASV table; inner edges = 0%; abundance increases outward; ",
      "months shown in the center;\nglobal Top 15 taxa; Others pooled; incertae sedis excluded; ",
      "connecting bands track the same taxon between adjacent months")) + theme_classic(base_size = 11.8,
      base_family = "Arial") + theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white",
      colour = NA), panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.28), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "#252525", linewidth = 0.42),
      axis.line.y = element_blank(), axis.ticks.x = element_line(colour = "#252525", linewidth = 0.3),
      axis.ticks.y = element_blank(), axis.text.x = element_text(colour = "#252525", size = 10.5),
      axis.text.y = element_blank(), axis.title.x = element_text(colour = "#202020", size = 11.8, face = "bold",
        margin = margin(t = 8)), plot.title = element_text(colour = "#171717", size = 14, face = "bold",
        lineheight = 1.08), plot.subtitle = element_text(colour = "#5A5A5A", size = 10.2, lineheight = 1.08,
        margin = margin(b = 9)), legend.position = "bottom", legend.direction = "vertical", legend.title.position = "top",
      legend.title = element_text(size = 10.2, face = "bold", hjust = 0.5), legend.text = element_text(size = 9.5,
        face = if (italic_legend)
          "italic"
        else "plain", colour = "#303030"), legend.key.width = grid::unit(3.6, "mm"), legend.key.height = grid::unit(3.2,
        "mm"), legend.spacing.x = grid::unit(2.5, "mm"), legend.box = "vertical", legend.box.just = "center",
      legend.box.margin = margin(t = 3), plot.margin = margin(t = 8, r = 8, b = 7, l = 8))
  }

  save_plot <- function(plot, file_base, legend_taxa, key_level = FALSE) {
    legend_rows <- legend_taxa
    height_mm <- min(330, max(260, 225 + legend_rows * 5.6))
    grDevices::png(paste0(file_base, ".png"), width = 230, height = height_mm, units = "mm", res = 400,
      type = if (capabilities("cairo"))
        "cairo-png"
      else "windows", bg = "white")
    print(plot)
    grDevices::dev.off()
    if (key_level) {
      grDevices::svg(paste0(file_base, ".svg"), width = 230/25.4, height = height_mm/25.4, family = "sans",
        bg = "white")
      print(plot)
      grDevices::dev.off()
      grDevices::tiff(paste0(file_base, ".tiff"), width = 230, height = height_mm, units = "mm", res = 600,
        compression = "lzw", type = if (capabilities("cairo"))
          "cairo"
        else "windows", bg = "white")
      print(plot)
      grDevices::dev.off()
    }
  }

  all_audits <- list()

  audit_index <- 0L

  for (dataset_name in names(dataset_definitions)) {
    definition <- dataset_definitions[[dataset_name]]
    message("Running monthly composition: ", dataset_name)
    dataset_audits <- list()
    colour_tables <- list()
    for (level_index in seq_along(taxonomic_levels)) {
      taxon_level <- taxonomic_levels[[level_index]]
      message("  Preparing ", taxon_level, "...")
      result <- prepare_monthly_level(definition, taxon_level)
      dataset_audits[[taxon_level]] <- result$audit
      audit_index <- audit_index + 1L
      all_audits[[audit_index]] <- result$audit
      if (is.null(result$data) || nrow(result$data) == 0L)
        next
      colour_map <- make_colour_map(result$data)
      colour_tables[[taxon_level]] <- data.frame(Taxonomic_Level = taxon_level, Taxon = names(colour_map),
        Colour = unname(colour_map), stringsAsFactors = FALSE)
      file_stem <- paste0(sprintf("%02d", level_index), "_", taxon_level, "_monthly")
      overall_file_stem <- paste0(file_stem, "_overall")
      overall_plot <- make_plot(plot_data = result$data, panel = "Overall", definition = definition,
        taxon_level = taxon_level, colour_map = colour_map)
      save_plot(overall_plot, file.path(definition$image_dir, overall_file_stem), legend_taxa = length(colour_map),
        key_level = taxon_level %in% c("phylum", "genus"))
      write.csv(result$data %>%
        filter(Panel == "Overall"), file.path(definition$table_dir, paste0(overall_file_stem, "_plot_data.csv")),
        row.names = FALSE, fileEncoding = "UTF-8")
      pyramid_file_stem <- paste0(file_stem, "_small_large_pyramid")
      pyramid_plot <- make_pyramid_plot(plot_data = result$data, definition = definition, taxon_level = taxon_level,
        colour_map = colour_map)
      save_plot(pyramid_plot, file.path(definition$image_dir, pyramid_file_stem), legend_taxa = length(colour_map),
        key_level = taxon_level %in% c("phylum", "genus"))
      write.csv(result$data %>%
        filter(Panel %in% c("USF", "ConF")), file.path(definition$table_dir, paste0(pyramid_file_stem,
        "_plot_data.csv")), row.names = FALSE, fileEncoding = "UTF-8")
      for (panel in c("USF", "ConF")) {
        write.csv(result$data %>%
          filter(Panel == panel), file.path(definition$table_dir, paste0(file_stem, "_", tolower(panel),
          "_plot_data.csv")), row.names = FALSE, fileEncoding = "UTF-8")
      }
      write.csv(result$data, file.path(definition$table_dir, paste0(file_stem, "_plot_data.csv")),
        row.names = FALSE, fileEncoding = "UTF-8")
    }
    write.csv(bind_rows(dataset_audits), file.path(definition$table_dir, "monthly_composition_audit.csv"),
      row.names = FALSE, fileEncoding = "UTF-8")
    if (length(colour_tables) > 0L) {
      write.csv(bind_rows(colour_tables), file.path(definition$table_dir, "monthly_colour_map.csv"),
        row.names = FALSE, fileEncoding = "UTF-8")
    }
    capture.output(sessionInfo(), file = file.path(definition$table_dir, "sessionInfo.txt"))
  }

  sampling_audit <- bind_rows(lapply(panel_levels, function(panel) {
    ids <- if (panel == "Overall") {
      sample_metadata$SampleID
    }
    else {
      sample_metadata$SampleID[sample_metadata$SizePanel == panel]
    }
    bind_rows(lapply(observed_months, function(month_date) {
      month_ids <- sample_metadata$SampleID[sample_metadata$MonthDate == month_date]
      data.frame(Panel = panel, Month = format(month_date, "%Y-%m"), SampleN = length(intersect(ids,
        month_ids)), stringsAsFactors = FALSE)
    }))
  }))

  combined_table_dir <- file.path(analysis_root, "结果表格", "抽平逐月版_总审计")

  dir.create(combined_table_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(bind_rows(all_audits), file.path(combined_table_dir, "01_all_dataset_level_audit.csv"), row.names = FALSE,
    fileEncoding = "UTF-8")

  write.csv(sampling_audit, file.path(combined_table_dir, "02_month_panel_sample_counts.csv"), row.names = FALSE,
    fileEncoding = "UTF-8")

  write.csv(data.frame(Metric = c("Rarefied_ASVs", "Samples", "Rarefaction_depth_min", "Rarefaction_depth_max",
    "Candidate_ASVs", "Sampled_months"), Value = c(nrow(abundance_all), ncol(abundance_all), min(sample_totals),
    max(sample_totals), length(candidate_ids), length(observed_months)), stringsAsFactors = FALSE), file.path(combined_table_dir,
    "03_rarefaction_and_scope_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  cat("\nAll rarefied monthly taxonomic-composition plots completed.\n")

  cat("Observed months: ", length(observed_months), "\n", sep = "")

  cat("Sample depth: ", min(sample_totals), " reads for every sample\n", sep = "")

  cat("Output root: ", normalizePath(analysis_root, winslash = "/"), "\n", sep = "")
})
