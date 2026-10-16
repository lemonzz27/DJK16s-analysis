local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  library(ggplot2)

  library(dplyr)

  library(tidyr)

  library(stringr)

  library(tibble)

  library(grid)

  TPM_FILE <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "00_输入文件/01_ASV与序列/ASV抽平21000_36097ASV_120样本.csv")

  GROUP_FILE <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "00_输入文件/02_样本分组与站点/孔径分组_groupsize.csv")

  TAXONOMY_FILE <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "01_分析项目/01_分类注释_SILVA_Greengenes/结果表格/SILVA_tax.csv")

  ID_COL <- "mags"

  PHYLUM_COL <- "phylum"

  GENUS_COL <- "genus"

  TOP_PHYLUM <- 5

  TOP_GENUS <- 3

  R_ROOT <- 0.8

  R_PHYLUM <- 2.25

  R_GENUS <- 3.75

  BASE_FONT <- "sans"

  GROUP_COLORS <- c(small = "#1B9E77", large = "#D95F02")

  PHYLUM_COLORS <- c("#3C5488", "#7666A5", "#A95676", "#A64F4A", "#85705A", "#66727E", "#58789C", "#8A6F9E",
    "#9D607B", "#59636E")

  tpm_data <- read.csv(TPM_FILE, row.names = 1, check.names = FALSE)

  group_data <- read.csv(GROUP_FILE, check.names = FALSE)

  taxonomy_raw <- read.csv(TAXONOMY_FILE, check.names = FALSE)

  colnames(tpm_data) <- gsub("\\.", "_", colnames(tpm_data))

  group_data$sample <- gsub("\\.", "_", group_data$sample)

  group_col <- intersect(c("group", "size_group", "Group", "Size_group"), colnames(group_data))[1]

  if (is.na(group_col)) {
    stop("groupsize.csv 中没有找到 group 或 size_group 列。")
  }

  group_data <- group_data %>%
    rename(Group = all_of(group_col)) %>%
    mutate(Group = tolower(Group)) %>%
    filter(Group %in% c("small", "large"))

  taxonomy <- taxonomy_raw %>%
    transmute(MAG = as.character(.data[[ID_COL]]), Phylum = as.character(.data[[PHYLUM_COL]]), Genus = as.character(.data[[GENUS_COL]])) %>%
    distinct(MAG, .keep_all = TRUE) %>%
    mutate(Phylum = str_remove(Phylum, "^p__"), Genus = str_remove(Genus, "^g__"), Phylum = str_trim(Phylum),
      Genus = str_trim(Genus))

  bad_taxa <- c("", "NA", "NaN", "unknown", "Unknown", "UNKNOWN", "unclassified", "Unclassified", "UNCLASSIFIED",
    "uncultured", "Uncultured")

  taxonomy <- taxonomy %>%
    mutate(Phylum = ifelse(Phylum %in% bad_taxa, NA, Phylum), Genus = ifelse(Genus %in% bad_taxa, NA,
      Genus))

  tpm_long <- tpm_data %>%
    rownames_to_column("MAG") %>%
    pivot_longer(cols = -MAG, names_to = "sample", values_to = "TPM") %>%
    mutate(TPM = as.numeric(TPM))

  taxon_data <- tpm_long %>%
    inner_join(group_data %>%
      select(sample, Group), by = "sample") %>%
    inner_join(taxonomy, by = "MAG")

  phylum_abundance <- taxon_data %>%
    filter(!is.na(Phylum)) %>%
    group_by(Group, Phylum) %>%
    summarise(Abundance = sum(TPM, na.rm = TRUE)/n_distinct(sample), .groups = "drop")

  genus_abundance <- taxon_data %>%
    filter(!is.na(Phylum), !is.na(Genus)) %>%
    filter(!str_detect(Genus, regex(paste(c("incertae[ _-]*sedis", "uncultured", "unclassified", "unidentified",
      "unknown", "norank", "no[ _-]*rank", "metagenome", "environmental[ _-]*sample", "candidate[ _-]*division",
      "uncultivated"), collapse = "|"), ignore_case = TRUE))) %>%
    group_by(Group, Phylum, Genus) %>%
    summarise(Abundance = sum(TPM, na.rm = TRUE)/n_distinct(sample), .groups = "drop")

  top_phylum <- phylum_abundance %>%
    group_by(Group) %>%
    slice_max(order_by = Abundance, n = TOP_PHYLUM, with_ties = FALSE) %>%
    ungroup()

  top_genus <- genus_abundance %>%
    semi_join(top_phylum, by = c("Group", "Phylum")) %>%
    group_by(Group, Phylum) %>%
    slice_max(order_by = Abundance, n = TOP_GENUS, with_ties = FALSE) %>%
    ungroup()

  write.csv(top_phylum, "Top5_Phylum_by_Group.csv", row.names = FALSE)

  write.csv(top_genus, "Top3_Genus_within_Top5_Phylum.csv", row.names = FALSE)

  phylum_levels <- unique(top_phylum$Phylum)

  if (length(phylum_levels) > length(PHYLUM_COLORS)) {
    stop("Top Phylum 数量超过当前颜色数量，请扩充 PHYLUM_COLORS。")
  }

  phylum_palette <- setNames(PHYLUM_COLORS[seq_along(phylum_levels)], phylum_levels)

  fill_palette <- c(GROUP_COLORS, phylum_palette)

  write.csv(data.frame(Phylum = names(phylum_palette), Colour = unname(phylum_palette), stringsAsFactors = FALSE),
    "Phylum_color_palette_CNS.csv", row.names = FALSE)

  root_nodes <- tibble(Group = c("small", "large"), node_id = c("root::small", "root::large"), label = c("USF",
    "ConF"), angle_deg = c(180, 0), radius = R_ROOT) %>%
    mutate(angle_rad = angle_deg * pi/180, x = radius * cos(angle_rad), y = radius * sin(angle_rad),
      Level = "Group")

  make_phylum_position <- function(df, group_name, angle_start, angle_end) {
    temp <- df %>%
      filter(Group == group_name) %>%
      arrange(desc(Abundance))
    if (nrow(temp) == 0) {
      return(temp)
    }
    temp$angle_deg <- seq(angle_start, angle_end, length.out = nrow(temp))
    temp
  }

  small_phylum <- make_phylum_position(top_phylum, "small", 125, 235)

  large_phylum <- make_phylum_position(top_phylum, "large", -55, 55)

  phylum_nodes <- bind_rows(small_phylum, large_phylum) %>%
    mutate(radius = R_PHYLUM, angle_rad = angle_deg * pi/180, x = radius * cos(angle_rad), y = radius *
      sin(angle_rad), node_id = paste("phy", Group, Phylum, sep = "::"), label = str_replace_all(Phylum,
      "_", " "), Level = "Phylum")

  genus_nodes <- top_genus %>%
    left_join(phylum_nodes %>%
      select(Group, Phylum, parent_angle = angle_deg), by = c("Group", "Phylum")) %>%
    group_by(Group, Phylum) %>%
    arrange(desc(Abundance), .by_group = TRUE) %>%
    mutate(genus_offset = if (n() == 1) {
      0
    } else {
      seq(-11, 11, length.out = n())
    }, angle_deg = parent_angle + genus_offset) %>%
    ungroup() %>%
    mutate(radius = R_GENUS, angle_rad = angle_deg * pi/180, x = radius * cos(angle_rad), y = radius *
      sin(angle_rad), node_id = paste("gen", Group, Phylum, Genus, sep = "::"), label = str_replace_all(Genus,
      "_", " "), Level = "Genus")

  edge_root_phylum <- phylum_nodes %>%
    transmute(from = paste0("root::", Group), to = node_id, Phylum = Phylum)

  edge_phylum_genus <- genus_nodes %>%
    transmute(from = paste("phy", Group, Phylum, sep = "::"), to = node_id, Phylum = Phylum)

  edges <- bind_rows(edge_root_phylum, edge_phylum_genus)

  all_nodes <- bind_rows(root_nodes %>%
    transmute(node_id, x, y), phylum_nodes %>%
    transmute(node_id, x, y), genus_nodes %>%
    transmute(node_id, x, y))

  edges_plot <- edges %>%
    left_join(all_nodes %>%
      rename(from = node_id, x_from = x, y_from = y), by = "from") %>%
    left_join(all_nodes %>%
      rename(to = node_id, x_to = x, y_to = y), by = "to")

  genus_labels <- genus_nodes %>%
    mutate(label_radius = R_GENUS + 0.32, label_x = label_radius * cos(angle_rad), label_y = label_radius *
      sin(angle_rad), angle_norm = (angle_deg%%360 + 360)%%360, text_angle = ifelse(angle_norm > 90 &
      angle_norm < 270, angle_norm + 180, angle_norm), hjust = ifelse(angle_norm > 90 & angle_norm <
      270, 1, 0), label = iconv(label, from = "", to = "UTF-8", sub = ""))

  p <- ggplot() + geom_segment(data = edges_plot, aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
    color = Phylum), linewidth = 0.65, alpha = 0.35, lineend = "round", show.legend = FALSE) + geom_point(data = phylum_nodes,
    aes(x = x, y = y, size = Abundance, fill = Phylum), shape = 21, color = "black", stroke = 0.85, alpha = 1) +
    geom_point(data = genus_nodes, aes(x = x, y = y, size = Abundance, fill = Phylum), shape = 21, color = "black",
      stroke = 0.7, alpha = 0.95) + geom_point(data = root_nodes, aes(x = x, y = y, fill = Group),
    shape = 21, size = 14, color = "black", stroke = 1.1, show.legend = FALSE) + geom_text(data = root_nodes,
    aes(x = x, y = y, label = label), family = BASE_FONT, fontface = "bold", size = 5.2, color = "white") +
    geom_text(data = genus_labels, aes(x = label_x, y = label_y, label = label, angle = text_angle, hjust = hjust),
      family = BASE_FONT, fontface = "italic", size = 3.6, color = "grey15") + scale_fill_manual(values = fill_palette,
    breaks = phylum_levels, name = "Phylum") + scale_color_manual(values = phylum_palette) + scale_size(range = c(4,
    20), name = "Mean abundance") + guides(fill = guide_legend(title = "Phylum", override.aes = list(size = 8,
    shape = 21, color = "black", stroke = 0.8, alpha = 1), order = 1), size = guide_legend(title = "Mean abundance",
    order = 2)) + coord_fixed(xlim = c(-4.8, 4.8), ylim = c(-4.6, 4.6), clip = "off") + labs(title = "Taxonomic composition of USF and ConF") +
    theme_void(base_family = BASE_FONT) + theme(plot.background = element_rect(fill = "transparent",
    color = NA), panel.background = element_rect(fill = "transparent", color = NA), plot.title = element_text(size = 20,
    face = "bold", hjust = 0.5, color = "black", margin = margin(b = 55)), legend.title = element_text(size = 13,
    face = "bold", color = "black"), legend.text = element_text(size = 12, color = "grey15"), legend.key.height = unit(0.75,
    "cm"), legend.key.width = unit(0.75, "cm"), legend.key = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA), legend.box.spacing = unit(2.6,
      "cm"), legend.margin = margin(l = 10), legend.position = "right", plot.margin = margin(20, 50,
      20, 40))

  graphics.off()

  print(p)

  ggsave(filename = "Taxonomic_Radial_Network_CNS.png", plot = p, width = 14.5, height = 10, dpi = 600,
    bg = "transparent")

  ggsave(filename = "Taxonomic_Radial_Network_CNS.pdf", plot = p, width = 14.5, height = 10, device = cairo_pdf,
    bg = "transparent")

  ggsave(filename = "Taxonomic_Radial_Network_CNS.svg", plot = p, width = 14.5, height = 10, device = svglite::svglite,
    bg = "transparent")

  ggsave(filename = "Taxonomic_Radial_Network_CNS.tiff", plot = p, width = 14.5, height = 10, dpi = 600,
    compression = "lzw", device = ragg::agg_tiff, bg = "transparent")

  cat("\n绘图完成：\n", "Taxonomic_Radial_Network_CNS.png\n", "Taxonomic_Radial_Network_CNS.pdf\n",
    "Taxonomic_Radial_Network_CNS.svg\n", "Taxonomic_Radial_Network_CNS.tiff\n")
})
