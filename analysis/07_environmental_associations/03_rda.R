local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  set.seed(123)

  .libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/figure2_hydro_core/Rlib"), .libPaths()))

  required_packages <- c("vegan", "ggplot2", "patchwork")

  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }

  suppressPackageStartupMessages({
    library(vegan)
    library(ggplot2)
    library(patchwork)
  })

  work_dir <- file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/rda_cca_noQ")

  input_dir <- file.path(work_dir, "input")

  output_root <- file.path(work_dir, "output")

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  whole_file <- file.path(input_dir, "whole.csv")

  candidate_file <- file.path(input_dir, "candidate.csv")

  env_file <- file.path(input_dir, "env.csv")

  group_file <- file.path(input_dir, "groupsize.csv")

  PERMUTATIONS <- 999

  DETAIL_PERMUTATIONS <- 199

  MIN_PREVALENCE_PROP <- 0.1

  RUN_SETS <- c("All", "small", "large")

  GROUP_COLORS <- c(All = "#2C7FB8", small = "#2A9D8F", large = "#E67E22")

  read_table <- function(path) {
    read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }

  whole_raw <- read_table(whole_file)

  candidate_raw <- read_table(candidate_file)

  env_raw <- read_table(env_file)

  group_raw <- read_table(group_file)

  whole_matrix <- as.matrix(whole_raw[, -1, drop = FALSE])

  mode(whole_matrix) <- "numeric"

  whole_depths <- colSums(whole_matrix, na.rm = TRUE)

  if (ncol(whole_matrix) != 120L || any(whole_depths != 21000)) {
    stop("Whole-community rarefaction audit failed.")
  }

  candidate_ids <- as.character(candidate_raw[[1]])

  whole_ids <- as.character(whole_raw[[1]])

  candidate_samples <- colnames(candidate_raw)[-1]

  whole_samples <- colnames(whole_raw)[-1]

  if (!identical(candidate_samples, whole_samples)) {
    stop("Candidate and whole-community sample columns differ.")
  }

  candidate_index <- match(candidate_ids, whole_ids)

  if (anyNA(candidate_index)) {
    stop("Candidate ASVs are not all present in the whole rarefied table.")
  }

  candidate_matrix <- as.matrix(candidate_raw[, candidate_samples, drop = FALSE])

  whole_subset <- as.matrix(whole_raw[candidate_index, whole_samples, drop = FALSE])

  mode(candidate_matrix) <- "numeric"

  mode(whole_subset) <- "numeric"

  if (!identical(unname(candidate_matrix), unname(whole_subset))) {
    stop("Candidate counts do not match the whole rarefied table.")
  }

  if (!all(c("sample", "group") %in% colnames(group_raw))) {
    stop("Group file must contain sample and group columns.")
  }

  if (!"sample" %in% colnames(env_raw)) {
    stop("Environmental file must contain a sample column.")
  }

  sig_stars <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  }

  make_hull <- function(df) {
    pieces <- lapply(split(df, df$plot_group), function(x) {
      x <- unique(x[, c("Axis1", "Axis2", "plot_group")])
      if (nrow(x) < 3)
        return(NULL)
      x[chull(x$Axis1, x$Axis2), , drop = FALSE]
    })
    do.call(rbind, pieces)
  }

  scale_arrows <- function(site_df, arrow_df, factor = 0.82) {
    if (nrow(arrow_df) == 0)
      return(arrow_df)
    sx <- diff(range(site_df$Axis1, na.rm = TRUE))
    sy <- diff(range(site_df$Axis2, na.rm = TRUE))
    ax <- diff(range(arrow_df$Axis1, na.rm = TRUE))
    ay <- diff(range(arrow_df$Axis2, na.rm = TRUE))
    mx <- ifelse(is.finite(ax) && ax > 0, sx/ax, 1)
    my <- ifelse(is.finite(ay) && ay > 0, sy/ay, 1)
    mult <- min(mx, my) * factor
    arrow_df$Axis1 <- arrow_df$Axis1 * mult
    arrow_df$Axis2 <- arrow_df$Axis2 * mult
    arrow_df
  }

  axis_percent <- function(ord) {
    eig <- vegan::eigenvals(ord, model = "constrained")
    if (length(eig) < 2 || sum(eig) <= 0)
      return(c(NA_real_, NA_real_))
    round(100 * eig[1:2]/sum(eig), 1)
  }

  extract_p <- function(test) {
    as.numeric(test$`Pr(>F)`[1])
  }

  format_p <- function(p) {
    ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  }

  save_plot_all_formats <- function(plot, stem, width, height) {
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(paste0(stem, ".svg"), plot, width = width, height = height)
    }
  }

  plot_ordination <- function(ord, env_use, meta_use, method, set_name, output_dir) {
    site <- as.data.frame(vegan::scores(ord, display = "sites", choices = 1:2, scaling = 2))
    colnames(site)[1:2] <- c("Axis1", "Axis2")
    site$sample <- rownames(site)
    site$size_group <- meta_use[site$sample, "group"]
    site$plot_group <- factor(set_name, levels = names(GROUP_COLORS))
    bp <- as.data.frame(vegan::scores(ord, display = "bp", choices = 1:2, scaling = 2))
    if (nrow(bp) > 0) {
      colnames(bp)[1:2] <- c("Axis1", "Axis2")
      bp$variable <- rownames(bp)
    }
    else {
      bp <- data.frame(Axis1 = numeric(), Axis2 = numeric(), variable = character())
    }
    fit <- vegan::envfit(ord, env_use, permutations = PERMUTATIONS, choices = 1:2)
    fit_df <- data.frame(variable = names(fit$vectors$r), r2 = as.numeric(fit$vectors$r), p_value = as.numeric(fit$vectors$pvals),
      stringsAsFactors = FALSE)
    fit_df$p_adj_BH <- p.adjust(fit_df$p_value, method = "BH")
    fit_df$significance <- sig_stars(fit_df$p_adj_BH)
    bp <- merge(bp, fit_df, by = "variable", all.x = TRUE)
    bp$label <- paste0(bp$variable, bp$significance)
    bp <- scale_arrows(site, bp)
    overall_test <- vegan::anova.cca(ord, permutations = PERMUTATIONS)
    term_test <- vegan::anova.cca(ord, by = "term", permutations = DETAIL_PERMUTATIONS)
    axis_test <- vegan::anova.cca(ord, by = "axis", permutations = DETAIL_PERMUTATIONS)
    r2 <- vegan::RsquareAdj(ord)
    pct <- axis_percent(ord)
    hull <- make_hull(site)
    plot_color <- GROUP_COLORS[[set_name]]
    subtitle <- paste0("R2 = ", sprintf("%.3f", r2$r.squared), "; adjusted R2 = ", sprintf("%.3f", r2$adj.r.squared),
      "; permutation P ", ifelse(extract_p(overall_test) < 0.001, "< 0.001", paste0("= ", format_p(extract_p(overall_test)))))
    p <- ggplot(site, aes(Axis1, Axis2)) + geom_hline(yintercept = 0, linewidth = 0.3, color = "grey72") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey72") + geom_polygon(data = hull, aes(fill = plot_group,
      group = plot_group), alpha = 0.14, color = NA, inherit.aes = TRUE) + geom_point(aes(fill = plot_group,
      color = plot_group), shape = 21, size = 3.4, stroke = 0.65, alpha = 0.88) + geom_segment(data = bp,
      aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, linewidth = 0.55, color = "black",
      arrow = grid::arrow(length = grid::unit(0.18, "cm"))) + geom_text(data = bp, aes(Axis1, Axis2,
      label = label), inherit.aes = FALSE, size = 3.2, color = "black", check_overlap = TRUE, vjust = -0.55) +
      scale_fill_manual(values = GROUP_COLORS, drop = FALSE) + scale_color_manual(values = GROUP_COLORS,
      drop = FALSE) + labs(title = paste0(ifelse(set_name == "small", "Small", ifelse(set_name == "large",
      "Large", "All")), " - ", method), subtitle = subtitle, x = paste0(method, "1 (", pct[1], "%)"),
      y = paste0(method, "2 (", pct[2], "%)"), fill = NULL, color = NULL) + theme_classic(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13), plot.subtitle = element_text(hjust = 0.5,
        size = 8.5), axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),
        legend.position = "none", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
        plot.margin = margin(6, 7, 6, 6))
    prefix <- file.path(output_dir, paste0(set_name, "_", method))
    save_plot_all_formats(p, prefix, 7.2, 5.4)
    write.csv(site, paste0(prefix, "_sample_scores.csv"), row.names = FALSE)
    write.csv(bp, paste0(prefix, "_env_scores.csv"), row.names = FALSE)
    write.csv(fit_df, paste0(prefix, "_envfit.csv"), row.names = FALSE)
    write.csv(as.data.frame(overall_test), paste0(prefix, "_anova_overall.csv"))
    write.csv(as.data.frame(term_test), paste0(prefix, "_anova_terms.csv"))
    write.csv(as.data.frame(axis_test), paste0(prefix, "_anova_axis.csv"))
    list(plot = p, r2 = r2$r.squared, adj_r2 = r2$adj.r.squared, p = extract_p(overall_test), axis1 = pct[1],
      axis2 = pct[2])
  }

  prepare_dataset <- function(raw_table) {
    taxa <- as.character(raw_table[[1]])
    mat <- raw_table[, -1, drop = FALSE]
    rownames(mat) <- taxa
    mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
    common <- Reduce(intersect, list(colnames(mat), env_raw$sample, group_raw$sample))
    common <- colnames(mat)[colnames(mat) %in% common]
    comm <- as.data.frame(t(mat[, common, drop = FALSE]))
    comm[is.na(comm)] <- 0
    env <- env_raw[match(common, env_raw$sample), , drop = FALSE]
    rownames(env) <- env$sample
    env <- env[, setdiff(colnames(env), c("sample", "Longitude", "Latitude", "Q")), drop = FALSE]
    env[] <- lapply(env, function(x) as.numeric(as.character(x)))
    env <- env[, colSums(!is.na(env)) > 0, drop = FALSE]
    env <- env[, vapply(env, function(x) sd(x, na.rm = TRUE) > 0, logical(1)), drop = FALSE]
    meta <- group_raw[match(common, group_raw$sample), , drop = FALSE]
    rownames(meta) <- meta$sample
    meta$group <- factor(meta$group, levels = c("small", "large"))
    keep <- complete.cases(env)
    list(comm = comm[keep, , drop = FALSE], env = env[keep, , drop = FALSE], meta = meta[keep, , drop = FALSE])
  }

  run_dataset <- function(raw_table, dataset_name, audit_row) {
    output_dir <- file.path(output_root, dataset_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    prepared <- prepare_dataset(raw_table)
    comm_all <- prepared$comm
    env_all <- prepared$env
    meta_all <- prepared$meta
    if ("Q" %in% colnames(env_all))
      stop("Q was not removed.")
    write.csv(data.frame(variable = colnames(env_all)), file.path(output_dir, "environment_variables_used.csv"),
      row.names = FALSE)
    write.csv(data.frame(test = c("overall_model", "envfit", "term_test", "axis_test"), permutations = c(PERMUTATIONS,
      PERMUTATIONS, DETAIL_PERMUTATIONS, DETAIL_PERMUTATIONS)), file.path(output_dir, "permutation_settings.csv"),
      row.names = FALSE)
    write.csv(audit_row, file.path(output_dir, "rarefied_input_audit.csv"), row.names = FALSE)
    results <- list()
    plots <- list()
    for (set_name in RUN_SETS) {
      if (set_name == "All") {
        samples <- rownames(comm_all)
      }
      else {
        samples <- rownames(meta_all)[meta_all$group == set_name]
      }
      comm <- comm_all[samples, , drop = FALSE]
      env_use <- env_all[samples, , drop = FALSE]
      meta_use <- meta_all[samples, , drop = FALSE]
      keep_sample <- rowSums(comm, na.rm = TRUE) > 0
      comm <- comm[keep_sample, , drop = FALSE]
      env_use <- env_use[keep_sample, , drop = FALSE]
      meta_use <- meta_use[keep_sample, , drop = FALSE]
      min_prevalence <- ceiling(nrow(comm) * MIN_PREVALENCE_PROP)
      keep_taxa <- colSums(comm > 0, na.rm = TRUE) >= min_prevalence & colSums(comm, na.rm = TRUE) >
        0
      if (sum(keep_taxa) < 10) {
        keep_taxa <- colSums(comm > 0, na.rm = TRUE) >= 2 & colSums(comm, na.rm = TRUE) > 0
      }
      comm <- comm[, keep_taxa, drop = FALSE]
      comm <- comm[, vapply(comm, function(x) sd(x, na.rm = TRUE) > 0, logical(1)), drop = FALSE]
      env_scaled <- as.data.frame(scale(env_use))
      rownames(env_scaled) <- rownames(env_use)
      comm_hel <- vegan::decostand(comm, method = "hellinger")
      comm_rel <- vegan::decostand(comm, method = "total")
      comm_rel[is.na(comm_rel)] <- 0
      dca <- vegan::decorana(comm_hel)
      dca_scores <- vegan::scores(dca, display = "sites", choices = 1:4)
      dca_lengths <- apply(dca_scores, 2, function(x) diff(range(x, na.rm = TRUE)))
      max_dca <- max(dca_lengths, na.rm = TRUE)
      write.csv(data.frame(axis = paste0("DCA", seq_along(dca_lengths)), gradient_length = dca_lengths),
        file.path(output_dir, paste0(set_name, "_DCA_gradient_length.csv")), row.names = FALSE)
      rda_model <- vegan::rda(comm_hel ~ ., data = env_scaled)
      cca_model <- vegan::cca(comm_rel ~ ., data = env_scaled)
      rda_result <- plot_ordination(rda_model, env_scaled, meta_use, "RDA", set_name, output_dir)
      cca_result <- plot_ordination(cca_model, env_scaled, meta_use, "CCA", set_name, output_dir)
      plots[[paste0(set_name, "_RDA")]] <- rda_result$plot
      plots[[paste0(set_name, "_CCA")]] <- cca_result$plot
      suggested <- ifelse(max_dca < 3, "RDA", ifelse(max_dca > 4, "CCA", "RDA_or_CCA"))
      results[[set_name]] <- rbind(data.frame(dataset = dataset_name, set = set_name, method = "RDA",
        n_samples = nrow(comm), n_taxa = ncol(comm), max_DCA_gradient = max_dca, suggested_method = suggested,
        R2 = rda_result$r2, adjusted_R2 = rda_result$adj_r2, p_value = rda_result$p, axis1_percent = rda_result$axis1,
        axis2_percent = rda_result$axis2), data.frame(dataset = dataset_name, set = set_name, method = "CCA",
        n_samples = nrow(comm), n_taxa = ncol(comm), max_DCA_gradient = max_dca, suggested_method = suggested,
        R2 = cca_result$r2, adjusted_R2 = cca_result$adj_r2, p_value = cca_result$p, axis1_percent = cca_result$axis1,
        axis2_percent = cca_result$axis2))
    }
    summary_table <- do.call(rbind, results)
    write.csv(summary_table, file.path(output_dir, "RDA_CCA_summary_noQ.csv"), row.names = FALSE)
    combined <- (plots$All_RDA + plots$All_CCA + plots$small_RDA + plots$small_CCA + plots$large_RDA +
      plots$large_CCA) + plot_layout(ncol = 2) + plot_annotation(title = paste0(ifelse(dataset_name ==
      "whole", "Whole community", "Pathogen candidates"), ": rarefied RDA/CCA without Q"), theme = theme(plot.title = element_text(hjust = 0.5,
      face = "bold", size = 16)))
    save_plot_all_formats(combined, file.path(output_dir, "RDA_CCA_combined_noQ"), 13.5, 16.2)
    summary_table
  }

  whole_audit <- data.frame(dataset = "whole", input = normalizePath(whole_file, winslash = "/"), n_taxa = nrow(whole_matrix),
    n_samples = ncol(whole_matrix), min_depth = min(whole_depths), median_depth = median(whole_depths),
    max_depth = max(whole_depths), rarefied_to_21000 = all(whole_depths == 21000), candidate_exact_subset = NA)

  candidate_depths <- colSums(candidate_matrix, na.rm = TRUE)

  candidate_audit <- data.frame(dataset = "candidate", input = normalizePath(candidate_file, winslash = "/"),
    n_taxa = nrow(candidate_matrix), n_samples = ncol(candidate_matrix), min_depth = min(candidate_depths),
    median_depth = median(candidate_depths), max_depth = max(candidate_depths), rarefied_to_21000 = FALSE,
    candidate_exact_subset = TRUE)

  whole_summary <- run_dataset(whole_raw, "whole", whole_audit)

  candidate_summary <- run_dataset(candidate_raw, "candidate", candidate_audit)

  all_summary <- rbind(whole_summary, candidate_summary)

  write.csv(all_summary, file.path(output_root, "RDA_CCA_all_results_noQ.csv"), row.names = FALSE)

  print(all_summary)
})

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  .libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/r_packages_extra"), file.path(Sys.getenv("DJK_PROJECT_ROOT",
    unset = "."), "work/figure2_hydro_core/Rlib"), .libPaths()))

  required_packages <- c("ggplot2", "ggrepel", "patchwork")

  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }

  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
  })

  get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1L) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
    }
    normalizePath(getwd(), winslash = "/")
  }

  script_dir <- get_script_dir()

  project_root <- Sys.getenv("DJK_PROJECT_ROOT", unset = ".")

  default_input_dir <- file.path(project_root, "01_分析项目", "14_RDA_CCA", "结果表格", "基础_注释无关",
    "整体_RDA_CCA_抽平去Q版", "outputs", "RDA_CCA_rarefied_noQ")

  default_output_dir <- file.path(project_root, "01_分析项目", "14_RDA_CCA", "结果图片", "基础_注释无关",
    "整体_RDA_CCA_抽平去Q版", "outputs", "RDA_CCA_rarefied_noQ", "正方形三组RDA")

  input_dir <- Sys.getenv("RDA_SQUARE_INPUT_DIR", unset = default_input_dir)

  output_dir <- Sys.getenv("RDA_SQUARE_OUTPUT_DIR", unset = default_output_dir)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  GROUP_COLORS <- c(All = "#2C7FB8", small = "#2A9D8F", large = "#E67E22")

  GROUP_LABELS <- c(All = "All", small = "USF", large = "ConF")

  SET_ORDER <- c("All", "large", "small")

  summary_table <- read.csv(file.path(input_dir, "RDA_CCA_summary_noQ.csv"), check.names = FALSE, stringsAsFactors = FALSE)

  summary_table <- summary_table[summary_table$method == "RDA", , drop = FALSE]

  make_hull <- function(df) {
    unique_xy <- unique(df[, c("Axis1", "Axis2")])
    if (nrow(unique_xy) < 3L)
      return(unique_xy)
    unique_xy[chull(unique_xy$Axis1, unique_xy$Axis2), , drop = FALSE]
  }

  make_rda_plot <- function(set_name, show_legend = TRUE) {
    site <- read.csv(file.path(input_dir, paste0(set_name, "_RDA_sample_scores.csv")), check.names = FALSE,
      stringsAsFactors = FALSE)
    env <- read.csv(file.path(input_dir, paste0(set_name, "_RDA_env_scores.csv")), check.names = FALSE,
      stringsAsFactors = FALSE)
    stats <- summary_table[summary_table$set == set_name, , drop = FALSE]
    if (nrow(stats) != 1L)
      stop("Missing RDA summary for set: ", set_name)
    site$plot_group <- factor(set_name, levels = names(GROUP_COLORS))
    hull <- make_hull(site)
    hull$plot_group <- factor(set_name, levels = names(GROUP_COLORS))
    x_range <- range(c(site$Axis1, env$Axis1), na.rm = TRUE)
    y_range <- range(c(site$Axis2, env$Axis2), na.rm = TRUE)
    x_pad <- diff(x_range) * 0.2
    y_pad <- diff(y_range) * 0.2
    x_limits <- x_range + c(-x_pad, x_pad)
    y_limits <- y_range + c(-y_pad, y_pad)
    env$arrow_length <- sqrt(env$Axis1^2 + env$Axis2^2)
    env$highlight <- FALSE
    eligible <- which(is.finite(env$p_adj_BH) & env$p_adj_BH < 0.05)
    if (length(eligible) > 0L) {
      highlight_order <- eligible[order(env$p_adj_BH[eligible], -env$arrow_length[eligible], env$variable[eligible])]
      env$highlight[head(highlight_order, 5L)] <- TRUE
    }
    shaft_fraction <- seq(0, 1.04, length.out = 50L)
    arrow_obstacles <- do.call(rbind, lapply(seq_len(nrow(env)), function(i) {
      data.frame(repel_x = env$Axis1[i] * shaft_fraction, repel_y = env$Axis2[i] * shaft_fraction,
        repel_label = "", stringsAsFactors = FALSE)
    }))
    label_points <- data.frame(repel_x = env$Axis1, repel_y = env$Axis2, repel_label = env$label, stringsAsFactors = FALSE)
    repel_data <- rbind(arrow_obstacles, label_points)
    stat_x <- x_limits[1] + 0.035 * diff(x_limits)
    stat_y <- y_limits[2] - 0.035 * diff(y_limits)
    stat_label <- sprintf("R^2 == %.3f*','~~P == %.3f", stats$R2, stats$p_value)
    legend_df <- data.frame(Axis1 = NA_real_, Axis2 = NA_real_, plot_group = factor(SET_ORDER, levels = names(GROUP_COLORS)))
    p <- ggplot(site, aes(Axis1, Axis2)) + geom_hline(yintercept = 0, linewidth = 0.3, color = "grey76") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey76") + geom_polygon(data = hull, aes(fill = plot_group,
      group = plot_group), alpha = 0.14, color = NA, inherit.aes = TRUE) + geom_point(aes(fill = plot_group,
      color = plot_group), shape = 21, size = 3.2, stroke = 0.55, alpha = 0.88) + geom_point(data = legend_df,
      aes(Axis1, Axis2, fill = plot_group, color = plot_group), inherit.aes = FALSE, shape = 21, size = 3.2,
      alpha = 0, na.rm = TRUE, show.legend = show_legend) + geom_segment(data = env[!env$highlight,
      , drop = FALSE], aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, linewidth = 0.55,
      color = "black", arrow = grid::arrow(length = grid::unit(0.17, "cm"))) + geom_segment(data = env[env$highlight,
      , drop = FALSE], aes(x = 0, y = 0, xend = Axis1, yend = Axis2), inherit.aes = FALSE, linewidth = 1.05,
      color = "#E6A700", arrow = grid::arrow(length = grid::unit(0.19, "cm"))) + ggrepel::geom_text_repel(data = repel_data,
      aes(repel_x, repel_y, label = repel_label), inherit.aes = FALSE, seed = 123, size = 3, color = "black",
      box.padding = grid::unit(0.34, "lines"), point.padding = grid::unit(0.18, "lines"), min.segment.length = 0,
      segment.color = "grey42", segment.size = 0.25, segment.alpha = 0.78, force = 2.2, force_pull = 0.58,
      max.overlaps = Inf, max.time = 18, max.iter = 2e+05, direction = "both", xlim = x_limits, ylim = y_limits) +
      annotate("text", x = stat_x, y = stat_y, label = stat_label, hjust = 0, vjust = 1, size = 4.1,
        fontface = "plain", parse = TRUE) + scale_fill_manual(values = GROUP_COLORS, breaks = SET_ORDER,
      labels = GROUP_LABELS[SET_ORDER], drop = FALSE) + scale_color_manual(values = GROUP_COLORS, breaks = SET_ORDER,
      labels = GROUP_LABELS[SET_ORDER], drop = FALSE) + scale_x_continuous(limits = x_limits, expand = c(0,
      0)) + scale_y_continuous(limits = y_limits, expand = c(0, 0)) + labs(x = sprintf("RDA1 (%.1f%%)",
      stats$axis1_percent), y = sprintf("RDA2 (%.1f%%)", stats$axis2_percent), fill = "Group", color = "Group") +
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1)), color = "none") +
      theme_classic(base_size = 11) + theme(aspect.ratio = 1, axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"), legend.position = if (show_legend)
        "bottom"
      else "none", legend.direction = "horizontal", legend.title = element_text(face = "bold"), panel.border = element_rect(color = "black",
        fill = NA, linewidth = 0.65), plot.margin = margin(7, 7, 7, 7))
    p
  }

  plots_no_legend <- lapply(SET_ORDER, make_rda_plot, show_legend = FALSE)

  names(plots_no_legend) <- SET_ORDER

  for (set_name in SET_ORDER) {
    ggsave(file.path(output_dir, paste0(set_name, "_RDA_square.png")), plots_no_legend[[set_name]], width = 6,
      height = 6, dpi = 600, bg = "white")
    ggsave(file.path(output_dir, paste0(set_name, "_RDA_square.pdf")), plots_no_legend[[set_name]], width = 6,
      height = 6, bg = "white")
  }

  plots_combined <- list(make_rda_plot("All", show_legend = TRUE), make_rda_plot("large", show_legend = FALSE),
    make_rda_plot("small", show_legend = FALSE))

  combined <- wrap_plots(plots_combined, nrow = 1, guides = "collect")

  ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.png"), combined, width = 15.5, height = 5.7,
    dpi = 600, bg = "white")

  ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.pdf"), combined, width = 15.5, height = 5.7,
    bg = "white")

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(output_dir, "RDA_three_groups_square_shared_legend.svg"), combined, width = 15.5,
      height = 5.7, bg = "white")
  }

  message("Square RDA panels saved to: ", output_dir)
})
