rm(list = ls())

options(stringsAsFactors = FALSE, warn = 1)

r_library <- Sys.getenv("NETWORK_R_LIB", unset = "")

if (nzchar(r_library)) {
  .libPaths(c(r_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(igraph)
  library(Hmisc)
  library(ggplot2)
})

get_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value))
    return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed))
    stop("Environment variable ", name, " must be numeric.")
  parsed
}

PREVALENCE_CUTOFF <- get_env_numeric("NETWORK_PREVALENCE", 0.2)

COR_R_CUTOFF <- get_env_numeric("NETWORK_R_CUTOFF", 0.7)

COR_Q_CUTOFF <- get_env_numeric("NETWORK_Q_CUTOFF", 0.001)

REMOVAL_REPS <- as.integer(get_env_numeric("NETWORK_REMOVAL_REPS", 100))

REMOVAL_STEP <- get_env_numeric("NETWORK_REMOVAL_STEP", 0.01)

RANDOM_SEED <- as.integer(get_env_numeric("NETWORK_RANDOM_SEED", 20260826))

OUTPUT_DIR <- Sys.getenv("NETWORK_OUTPUT_DIR", unset = "Results_prevalence20_enhanced")

DATASET_LABEL <- Sys.getenv("NETWORK_DATASET_LABEL", unset = "Overall_ASV")

if (PREVALENCE_CUTOFF <= 0 || PREVALENCE_CUTOFF > 1) {
  stop("NETWORK_PREVALENCE must be in (0, 1].")
}

if (COR_R_CUTOFF < 0 || COR_R_CUTOFF > 1) {
  stop("NETWORK_R_CUTOFF must be in [0, 1].")
}

if (COR_Q_CUTOFF <= 0 || COR_Q_CUTOFF > 1) {
  stop("NETWORK_Q_CUTOFF must be in (0, 1].")
}

if (REMOVAL_REPS < 1) stop("NETWORK_REMOVAL_REPS must be at least 1.")

if (REMOVAL_STEP <= 0 || REMOVAL_STEP > 1) {
  stop("NETWORK_REMOVAL_STEP must be in (0, 1].")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(OUTPUT_DIR, "Plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nStatic network analysis started: ", format(Sys.time()), "\n", sep = "")

cat("Output directory: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")

required_files <- c("tpm.csv", "groupsize.csv", "tax.csv", "Ultramicrobacteria.csv")

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

tpm_data <- read.csv("tpm.csv", row.names = 1, check.names = FALSE)

group_info <- read.csv("groupsize.csv", check.names = FALSE)

tax_info <- read.csv("tax.csv", check.names = FALSE)

ultra_info <- read.csv("Ultramicrobacteria.csv", check.names = FALSE)

names(group_info)[1:2] <- c("sample", "group")

names(tax_info) <- trimws(names(tax_info))

names(ultra_info) <- trimws(names(ultra_info))

if (!"mags" %in% names(tax_info)) stop("tax.csv must contain a 'mags' column.")

if (!all(c("mags", "Ultra_Group") %in% names(ultra_info))) {
  stop("Ultramicrobacteria.csv must contain 'mags' and 'Ultra_Group' columns.")
}

ultra_match <- match(tax_info$mags, ultra_info$mags)

tax_info$is_Ultramicrobacteria <- ultra_info$Ultra_Group[ultra_match]

tax_info$is_Ultramicrobacteria[is.na(tax_info$is_Ultramicrobacteria) | tax_info$is_Ultramicrobacteria ==
  ""] <- "Others"

rownames(tax_info) <- tax_info$mags

common_ids <- intersect(rownames(tpm_data), rownames(tax_info))

if (length(common_ids) < 2) stop("Fewer than two ASVs remain after taxonomy matching.")

tpm_data <- tpm_data[common_ids, , drop = FALSE]

tax_info <- tax_info[common_ids, , drop = FALSE]

numeric_tpm <- as.matrix(tpm_data)

storage.mode(numeric_tpm) <- "double"

if (any(!is.finite(numeric_tpm))) stop("tpm.csv contains non-finite abundance values.")

if (any(numeric_tpm < 0)) stop("tpm.csv contains negative abundance values.")

tpm_data <- numeric_tpm

groups_available <- unique(tolower(trimws(group_info$group)))

target_groups <- c("small", "large")

group_display <- c(small = "USF", large = "ConF")

relabel_group_for_output <- function(data) {
  if ("Group" %in% names(data)) {
    data$Group <- unname(group_display[as.character(data$Group)])
  }
  data
}

missing_groups <- setdiff(target_groups, groups_available)

if (length(missing_groups) > 0) {
  stop("groupsize.csv is missing group(s): ", paste(missing_groups, collapse = ", "))
}

group_info$group <- tolower(trimws(group_info$group))

safe_number <- function(x) {
  ifelse(is.finite(x), as.numeric(x), NA_real_)
}

calc_natural_connectivity <- function(g, weighted = FALSE) {
  n <- vcount(g)
  if (n == 0)
    return(0)
  if (n == 1)
    return(0)
  if (weighted) {
    adjacency <- as.matrix(as_adjacency_matrix(g, attr = "Weight", sparse = FALSE))
  }
  else {
    adjacency <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
    adjacency[adjacency != 0] <- 1
  }
  eigenvalues <- eigen(adjacency, symmetric = TRUE, only.values = TRUE)$values
  shift <- max(eigenvalues)
  safe_number(shift + log(mean(exp(eigenvalues - shift))))
}

calc_global_efficiency <- function(g) {
  n <- vcount(g)
  if (n < 2)
    return(0)
  distance_matrix <- distances(g, mode = "all", weights = NA)
  inverse_distance <- 1/distance_matrix
  inverse_distance[!is.finite(inverse_distance)] <- 0
  diag(inverse_distance) <- 0
  safe_number(sum(inverse_distance)/(n * (n - 1)))
}

calc_zi_pi <- function(g, modules) {
  node_names <- V(g)$name
  modules <- modules[node_names]
  total_degree <- degree(g, mode = "all")
  result <- data.frame(ID = node_names, Zi = 0, Pi = 0, Role = "Peripheral", row.names = node_names,
    check.names = FALSE)
  module_internal_degrees <- lapply(unique(modules), function(module_id) {
    nodes_in_module <- names(modules)[modules == module_id]
    values <- vapply(nodes_in_module, function(node) {
      sum(names(neighbors(g, node, mode = "all")) %in% nodes_in_module)
    }, numeric(1))
    names(values) <- nodes_in_module
    values
  })
  names(module_internal_degrees) <- as.character(unique(modules))
  for (node in node_names) {
    module_id <- modules[[node]]
    internal_degrees <- module_internal_degrees[[as.character(module_id)]]
    module_sd <- stats::sd(internal_degrees)
    result[node, "Zi"] <- if (!is.na(module_sd) && module_sd > 0) {
      (internal_degrees[[node]] - mean(internal_degrees))/module_sd
    }
    else {
      0
    }
    node_degree <- total_degree[[node]]
    if (node_degree > 0) {
      neighbor_modules <- modules[names(neighbors(g, node, mode = "all"))]
      result[node, "Pi"] <- 1 - sum((table(neighbor_modules)/node_degree)^2)
    }
  }
  result$Role <- ifelse(result$Zi < 2.5 & result$Pi < 0.62, "Peripheral", ifelse(result$Zi < 2.5 &
    result$Pi >= 0.62, "Connector", ifelse(result$Zi >= 2.5 & result$Pi < 0.62, "Module Hub", "Network Hub")))
  result
}

calculate_global_metrics <- function(g, group_name, sample_count, eligible_count, modules) {
  n <- vcount(g)
  m <- ecount(g)
  component_info <- components(g, mode = "weak")
  largest_component_id <- which.max(component_info$csize)
  largest_nodes <- V(g)$name[component_info$membership == largest_component_id]
  largest_graph <- induced_subgraph(g, vids = largest_nodes)
  degrees <- degree(g)
  edge_r <- E(g)$r
  module_count <- length(unique(modules))
  global_transitivity <- suppressWarnings(transitivity(g, type = "global", isolates = "zero"))
  average_transitivity <- suppressWarnings(transitivity(g, type = "average", isolates = "zero"))
  degree_assortativity <- suppressWarnings(assortativity_degree(g, directed = FALSE))
  modularity_value <- suppressWarnings(modularity(g, modules, weights = E(g)$Weight))
  mean_path <- if (vcount(largest_graph) >= 2) {
    mean_distance(largest_graph, directed = FALSE, unconnected = FALSE, weights = NA)
  }
  else {
    0
  }
  graph_diameter <- if (vcount(largest_graph) >= 2) {
    diameter(largest_graph, directed = FALSE, unconnected = FALSE, weights = NA)
  }
  else {
    0
  }
  data.frame(Dataset = DATASET_LABEL, Group = group_name, Samples = sample_count, Input_ASVs_after_taxonomy_match = nrow(tpm_data),
    Prevalence_eligible_ASVs = eligible_count, Network_nodes_nonisolated = n, Isolated_eligible_ASVs_excluded = eligible_count -
      n, Edges = m, Positive_edges = sum(edge_r > 0), Negative_edges = sum(edge_r < 0), Positive_edge_fraction = if (m >
      0)
      sum(edge_r > 0)/m
    else NA_real_, Negative_edge_fraction = if (m > 0)
      sum(edge_r < 0)/m
    else NA_real_, Average_degree = if (n > 0)
      mean(degrees)
    else 0, Median_degree = if (n > 0)
      stats::median(degrees)
    else 0, Maximum_degree = if (n > 0)
      max(degrees)
    else 0, Density_network_nodes = if (n >= 2)
      edge_density(g, loops = FALSE)
    else 0, Density_all_prevalence_eligible_ASVs = if (eligible_count >= 2) {
      2 * m/(eligible_count * (eligible_count - 1))
    }
    else {
      0
    }, Connected_components = component_info$no, Largest_component_nodes = max(component_info$csize),
    Largest_component_fraction = max(component_info$csize)/n, Number_of_modules = module_count, Modularity_abs_r = safe_number(modularity_value),
    Global_clustering_coefficient = safe_number(global_transitivity), Average_local_clustering = safe_number(average_transitivity),
    Degree_assortativity = safe_number(degree_assortativity), Mean_shortest_path_largest_component = safe_number(mean_path),
    Diameter_largest_component = safe_number(graph_diameter), Global_efficiency = calc_global_efficiency(g),
    Natural_connectivity_unweighted = calc_natural_connectivity(g, weighted = FALSE), Natural_connectivity_abs_r_weighted = calc_natural_connectivity(g,
      weighted = TRUE), Mean_absolute_r = if (m > 0)
      mean(abs(edge_r))
    else NA_real_, Degree_centralization = safe_number(centr_degree(g, mode = "all", normalized = TRUE)$centralization),
    stringsAsFactors = FALSE, check.names = FALSE)
}

build_group_network <- function(group_name) {
  samples_requested <- group_info$sample[group_info$group == group_name]
  samples <- intersect(samples_requested, colnames(tpm_data))
  if (length(samples) < 5)
    stop("Too few matched samples for group ", group_name, ".")
  abundance <- tpm_data[, samples, drop = FALSE]
  minimum_samples <- ceiling(length(samples) * PREVALENCE_CUTOFF)
  prevalence_count <- rowSums(abundance > 0)
  eligible <- prevalence_count >= minimum_samples
  filtered <- abundance[eligible, , drop = FALSE]
  if (nrow(filtered) < 2)
    stop("Fewer than two ASVs pass prevalence filtering for ", group_name, ".")
  cat("\n[", group_name, "] samples=", length(samples), ", minimum detections=", minimum_samples, ", prevalence-passing ASVs=",
    nrow(filtered), "\n", sep = "")
  cat("[", group_name, "] calculating Spearman correlations and BH q-values...\n", sep = "")
  correlation <- Hmisc::rcorr(t(filtered), type = "spearman")
  r_matrix <- correlation$r
  p_matrix <- correlation$P
  r_matrix[is.na(r_matrix)] <- 0
  p_matrix[is.na(p_matrix)] <- 1
  upper_mask <- upper.tri(r_matrix)
  q_matrix <- matrix(1, nrow = nrow(r_matrix), ncol = ncol(r_matrix), dimnames = dimnames(r_matrix))
  q_matrix[upper_mask] <- p.adjust(p_matrix[upper_mask], method = "BH")
  q_matrix[lower.tri(q_matrix)] <- t(q_matrix)[lower.tri(q_matrix)]
  diag(q_matrix) <- 1
  edge_index <- which(upper_mask & abs(r_matrix) >= COR_R_CUTOFF & q_matrix < COR_Q_CUTOFF, arr.ind = TRUE)
  if (nrow(edge_index) == 0)
    stop("No significant edges remain for group ", group_name, ".")
  edge_values <- r_matrix[edge_index]
  edges <- data.frame(Source = rownames(r_matrix)[edge_index[, 1]], Target = colnames(r_matrix)[edge_index[,
    2]], Weight = abs(edge_values), r = edge_values, p_raw = p_matrix[edge_index], q_BH = q_matrix[edge_index],
    Direction = ifelse(edge_values > 0, "Positive", "Negative"), Group = group_name, stringsAsFactors = FALSE,
    check.names = FALSE)
  eligible_ids <- rownames(filtered)
  graph_all <- graph_from_data_frame(edges[, c("Source", "Target", "Weight", "r", "p_raw", "q_BH",
    "Direction", "Group")], directed = FALSE, vertices = data.frame(name = eligible_ids))
  graph_network <- delete_vertices(graph_all, which(degree(graph_all) == 0))
  if (vcount(graph_network) < 2)
    stop("Fewer than two non-isolated nodes remain for ", group_name, ".")
  community <- cluster_fast_greedy(graph_network, weights = E(graph_network)$Weight)
  modules <- membership(community)
  zi_pi <- calc_zi_pi(graph_network, modules)
  node_ids <- V(graph_network)$name
  component_info <- components(graph_network)
  local_clustering <- transitivity(graph_network, type = "local", isolates = "zero")
  names(local_clustering) <- node_ids
  eigenvector_values <- eigen_centrality(graph_network, directed = FALSE, weights = E(graph_network)$Weight)$vector
  pagerank_values <- page_rank(graph_network, directed = FALSE, weights = E(graph_network)$Weight)$vector
  nodes <- data.frame(Id = node_ids, Label = node_ids, Degree = degree(graph_network)[node_ids], Betweenness = betweenness(graph_network,
    v = node_ids, directed = FALSE, normalized = TRUE, weights = NA), Closeness = closeness(graph_network,
    vids = node_ids, normalized = TRUE, weights = NA), Module = as.numeric(modules[node_ids]), Zi = zi_pi[node_ids,
    "Zi"], Pi = zi_pi[node_ids, "Pi"], Keystone_Role = zi_pi[node_ids, "Role"], mean = rowMeans(filtered)[node_ids],
    Group = group_name, StrengthAbs = strength(graph_network, vids = node_ids, weights = E(graph_network)$Weight),
    HarmonicCentrality = harmonic_centrality(graph_network, vids = node_ids, normalized = TRUE, weights = NA),
    Eigenvector = eigenvector_values[node_ids], PageRank = pagerank_values[node_ids], Coreness = coreness(graph_network,
      mode = "all")[node_ids], LocalClustering = local_clustering[node_ids], Component = component_info$membership[node_ids],
    PrevalenceCount = prevalence_count[node_ids], PrevalenceFraction = prevalence_count[node_ids]/length(samples),
    stringsAsFactors = FALSE, check.names = FALSE)
  taxonomy_columns <- setdiff(names(tax_info), "mags")
  nodes <- cbind(nodes, tax_info[node_ids, taxonomy_columns, drop = FALSE])
  edge_order <- c("Source", "Target", "Weight", "r", "Direction", "Group", "p_raw", "q_BH")
  edges <- edges[, edge_order, drop = FALSE]
  global_metrics <- calculate_global_metrics(graph_network, group_name, length(samples), nrow(filtered),
    modules)
  audit <- data.frame(Dataset = DATASET_LABEL, Group = group_name, Samples = length(samples), Prevalence_cutoff = PREVALENCE_CUTOFF,
    Minimum_positive_samples = minimum_samples, Input_ASVs_after_taxonomy_match = nrow(tpm_data),
    Prevalence_eligible_ASVs = nrow(filtered), Network_nodes_nonisolated = vcount(graph_network),
    Isolated_eligible_ASVs_excluded = nrow(filtered) - vcount(graph_network), Edges = ecount(graph_network),
    stringsAsFactors = FALSE)
  cat("[", group_name, "] network nodes=", vcount(graph_network), ", edges=", ecount(graph_network),
    "\n", sep = "")
  list(group = group_name, samples = samples, eligible_ids = eligible_ids, graph = graph_network, nodes = nodes,
    edges = edges, global = global_metrics, audit = audit)
}

removal_metrics <- function(g, removed_nodes, baseline_nodes) {
  total_nodes <- vcount(g)
  remove_count <- length(removed_nodes)
  if (remove_count >= total_nodes) {
    return(c(RemainingNodes = 0, RemainingEdges = 0, AverageDegree = 0, NaturalConnectivity = 0,
      LCCNodes = 0, LCCFraction = 0))
  }
  if (remove_count == 0) {
    reduced <- g
  }
  else {
    reduced <- delete_vertices(g, removed_nodes)
  }
  remaining_nodes <- vcount(reduced)
  remaining_edges <- ecount(reduced)
  component_info <- components(reduced)
  lcc_nodes <- if (remaining_nodes > 0)
    max(component_info$csize)
  else 0
  c(RemainingNodes = remaining_nodes, RemainingEdges = remaining_edges, AverageDegree = if (remaining_nodes >
    0) mean(degree(reduced)) else 0, NaturalConnectivity = calc_natural_connectivity(reduced, weighted = FALSE),
    LCCNodes = lcc_nodes, LCCFraction = lcc_nodes/baseline_nodes)
}

simulate_node_removal <- function(network_result, group_seed) {
  g <- network_result$graph
  fractions <- unique(c(seq(0, 1, by = REMOVAL_STEP), 0.5, 1))
  fractions <- sort(round(fractions[fractions >= 0 & fractions <= 1], 10))
  set.seed(group_seed)
  output <- vector("list", length(fractions) * REMOVAL_REPS)
  position <- 1L
  total_nodes <- vcount(g)
  baseline_result <- removal_metrics(g, character(0), total_nodes)
  fully_removed_result <- removal_metrics(g, V(g)$name, total_nodes)
  cat("[", network_result$group, "] random node-removal simulation: ", length(fractions), " fractions x ",
    REMOVAL_REPS, " repetitions...\n", sep = "")
  for (replicate_id in seq_len(REMOVAL_REPS)) {
    removal_order <- sample(V(g)$name, size = total_nodes, replace = FALSE)
    for (fraction in fractions) {
      remove_count <- min(total_nodes, round(fraction * total_nodes))
      one_result <- if (remove_count == 0) {
        baseline_result
      }
      else if (remove_count >= total_nodes) {
        fully_removed_result
      }
      else {
        removal_metrics(g, removed_nodes = removal_order[seq_len(remove_count)], baseline_nodes = total_nodes)
      }
      output[[position]] <- data.frame(Dataset = DATASET_LABEL, Group = network_result$group, RemovalFraction = fraction,
        RemovalPercent = fraction * 100, Replicate = replicate_id, as.list(one_result), stringsAsFactors = FALSE,
        check.names = FALSE)
      position <- position + 1L
    }
    if (replicate_id%%10L == 0L || replicate_id == REMOVAL_REPS) {
      cat("  ", network_result$group, ": ", replicate_id, "/", REMOVAL_REPS, " sequential trajectories completed\n",
        sep = "")
    }
  }
  do.call(rbind, output)
}

summarize_removal <- function(raw_data) {
  split_key <- interaction(raw_data$Group, raw_data$RemovalFraction, drop = TRUE)
  pieces <- lapply(split(seq_len(nrow(raw_data)), split_key), function(index) {
    block <- raw_data[index, , drop = FALSE]
    summarize_metric <- function(metric) {
      values <- block[[metric]]
      c(Mean = mean(values), SD = stats::sd(values), Lower95 = unname(stats::quantile(values, 0.025,
        names = FALSE)), Upper95 = unname(stats::quantile(values, 0.975, names = FALSE)))
    }
    average_degree <- summarize_metric("AverageDegree")
    natural_connectivity <- summarize_metric("NaturalConnectivity")
    lcc_fraction <- summarize_metric("LCCFraction")
    data.frame(Dataset = block$Dataset[1], Group = block$Group[1], RemovalFraction = block$RemovalFraction[1],
      RemovalPercent = block$RemovalPercent[1], Repetitions = nrow(block), AverageDegree_Mean = average_degree[["Mean"]],
      AverageDegree_SD = average_degree[["SD"]], AverageDegree_Lower95 = average_degree[["Lower95"]],
      AverageDegree_Upper95 = average_degree[["Upper95"]], NaturalConnectivity_Mean = natural_connectivity[["Mean"]],
      NaturalConnectivity_SD = natural_connectivity[["SD"]], NaturalConnectivity_Lower95 = natural_connectivity[["Lower95"]],
      NaturalConnectivity_Upper95 = natural_connectivity[["Upper95"]], LCCFraction_Mean = lcc_fraction[["Mean"]],
      LCCFraction_SD = lcc_fraction[["SD"]], LCCFraction_Lower95 = lcc_fraction[["Lower95"]], LCCFraction_Upper95 = lcc_fraction[["Upper95"]],
      stringsAsFactors = FALSE, check.names = FALSE)
  })
  result <- do.call(rbind, pieces)
  result[order(result$Group, result$RemovalFraction), , drop = FALSE]
}

make_node_long_table <- function(node_table) {
  metric_columns <- c("Degree", "StrengthAbs", "Betweenness", "Closeness", "HarmonicCentrality", "Eigenvector",
    "PageRank", "Coreness", "LocalClustering", "Zi", "Pi", "mean", "PrevalenceCount", "PrevalenceFraction")
  do.call(rbind, lapply(metric_columns, function(metric) {
    data.frame(Dataset = DATASET_LABEL, Group = node_table$Group, Id = node_table$Id, Metric = metric,
      Value = node_table[[metric]], stringsAsFactors = FALSE)
  }))
}

compare_50_percent <- function(data_50) {
  metrics <- c("AverageDegreeRetention", "NaturalConnectivityRetention", "LCCFraction")
  comparisons <- lapply(metrics, function(metric) {
    large_values <- data_50[data_50$Group == "large", metric]
    small_values <- data_50[data_50$Group == "small", metric]
    test <- suppressWarnings(wilcox.test(large_values, small_values, exact = FALSE))
    data.frame(Metric = metric, ConF_median = median(large_values), USF_median = median(small_values),
      Difference_ConF_minus_USF = median(large_values) - median(small_values), Wilcoxon_W = unname(test$statistic),
      P_value = test$p.value, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, comparisons)
  result$BH_q_value <- p.adjust(result$P_value, method = "BH")
  result
}

plot_stability_curves <- function(summary_data) {
  average_degree <- data.frame(Group = summary_data$Group, RemovalPercent = summary_data$RemovalPercent,
    Metric = "Average degree", Mean = summary_data$AverageDegree_Mean, Lower95 = summary_data$AverageDegree_Lower95,
    Upper95 = summary_data$AverageDegree_Upper95)
  natural_connectivity <- data.frame(Group = summary_data$Group, RemovalPercent = summary_data$RemovalPercent,
    Metric = "Natural connectivity", Mean = summary_data$NaturalConnectivity_Mean, Lower95 = summary_data$NaturalConnectivity_Lower95,
    Upper95 = summary_data$NaturalConnectivity_Upper95)
  plot_data <- rbind(average_degree, natural_connectivity)
  plot_data$Metric <- factor(plot_data$Metric, levels = c("Average degree", "Natural connectivity"))
  plot_data$Group <- unname(group_display[as.character(plot_data$Group)])
  plot_data$Group <- factor(plot_data$Group, levels = c("USF", "ConF"))
  palette <- c(USF = "#009E73", ConF = "#D55E00")
  p <- ggplot(plot_data, aes(x = RemovalPercent, y = Mean, color = Group, fill = Group)) + geom_ribbon(aes(ymin = Lower95,
    ymax = Upper95), alpha = 0.16, linewidth = 0) + geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~Metric, ncol = 1, scales = "free_y") + scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), expand = c(0.01, 0)) + labs(x = "Randomly removed nodes (%)",
    y = NULL, color = "Group", fill = "Group", title = "Static network stability under random node removal",
    subtitle = paste0("Mean and empirical 95% interval; ", REMOVAL_REPS, " repetitions per removal proportion")) +
    theme_classic(base_size = 11) + theme(strip.background = element_rect(fill = "grey95", color = "grey70",
    linewidth = 0.4), strip.text = element_text(face = "bold"), legend.position = "top", plot.title = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines"))
  ggsave(file.path(plot_dir, "Static_Network_Node_Removal_Stability.png"), p, width = 7.2, height = 7.2,
    dpi = 600, bg = "white")
  ggsave(file.path(plot_dir, "Static_Network_Node_Removal_Stability.pdf"), p, width = 7.2, height = 7.2,
    device = grDevices::cairo_pdf)
  ggsave(file.path(plot_dir, "Static_Network_Node_Removal_Stability.svg"), p, width = 7.2, height = 7.2,
    device = svglite::svglite)
}

plot_50_percent_boxplot <- function(data_50) {
  plot_data <- rbind(data.frame(Group = data_50$Group, Replicate = data_50$Replicate, Metric = "Average degree retention",
    Value = data_50$AverageDegreeRetention), data.frame(Group = data_50$Group, Replicate = data_50$Replicate,
    Metric = "Natural connectivity retention", Value = data_50$NaturalConnectivityRetention), data.frame(Group = data_50$Group,
    Replicate = data_50$Replicate, Metric = "Largest-component node fraction", Value = data_50$LCCFraction))
  plot_data$Group <- unname(group_display[as.character(plot_data$Group)])
  plot_data$Group <- factor(plot_data$Group, levels = c("USF", "ConF"))
  plot_data$Metric <- factor(plot_data$Metric, levels = c("Average degree retention", "Natural connectivity retention",
    "Largest-component node fraction"))
  palette <- c(USF = "#009E73", ConF = "#D55E00")
  p <- ggplot(plot_data, aes(x = Group, y = Value, color = Group, fill = Group)) + geom_boxplot(width = 0.58,
    alpha = 0.22, outlier.shape = NA, linewidth = 0.55) + geom_jitter(width = 0.12, height = 0, size = 0.9,
    alpha = 0.45) + facet_wrap(~Metric, nrow = 1, scales = "free_y") + scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) + labs(x = NULL, y = "Residual value / baseline", title = "Network robustness after 50% random node removal",
    subtitle = paste0(REMOVAL_REPS, " random removals per group")) + theme_classic(base_size = 11) +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95", color = "grey70",
      linewidth = 0.4), strip.text = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold"))
  ggsave(file.path(plot_dir, "Static_Network_Robustness_50pct.png"), p, width = 10, height = 4.2, dpi = 600,
    bg = "white")
  ggsave(file.path(plot_dir, "Static_Network_Robustness_50pct.pdf"), p, width = 10, height = 4.2, device = grDevices::cairo_pdf)
  ggsave(file.path(plot_dir, "Static_Network_Robustness_50pct.svg"), p, width = 10, height = 4.2, device = svglite::svglite)
}

results <- lapply(target_groups, build_group_network)

names(results) <- target_groups

for (group_name in target_groups) {
  output_group <- unname(group_display[group_name])
  write.csv(relabel_group_for_output(results[[group_name]]$nodes), file.path(OUTPUT_DIR, paste0("Node_",
    output_group, ".csv")), row.names = FALSE, quote = TRUE, na = "")
  write.csv(relabel_group_for_output(results[[group_name]]$edges), file.path(OUTPUT_DIR, paste0("Edge_",
    output_group, ".csv")), row.names = FALSE, quote = TRUE, na = "")
}

node_wide <- do.call(rbind, lapply(results, `[[`, "nodes"))

node_long <- make_node_long_table(node_wide)

global_metrics <- do.call(rbind, lapply(results, `[[`, "global"))

filter_audit <- do.call(rbind, lapply(results, `[[`, "audit"))

write.csv(relabel_group_for_output(node_wide), file.path(OUTPUT_DIR, "Network_Node_Metrics_Wide.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(relabel_group_for_output(node_long), file.path(OUTPUT_DIR, "Network_Node_Metrics_Long.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(relabel_group_for_output(global_metrics), file.path(OUTPUT_DIR, "Network_Global_Metrics.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(relabel_group_for_output(filter_audit), file.path(OUTPUT_DIR, "Feature_Filtering_Audit.csv"),
  row.names = FALSE, quote = TRUE, na = "")

removal_raw <- rbind(simulate_node_removal(results$small, RANDOM_SEED + 1L), simulate_node_removal(results$large,
  RANDOM_SEED + 2L))

removal_summary <- summarize_removal(removal_raw)

baseline_average_degree <- setNames(global_metrics$Average_degree, global_metrics$Group)

baseline_natural_connectivity <- setNames(global_metrics$Natural_connectivity_unweighted, global_metrics$Group)

data_50 <- removal_raw[abs(removal_raw$RemovalFraction - 0.5) < 1e-09, , drop = FALSE]

data_50$AverageDegreeRetention <- data_50$AverageDegree/baseline_average_degree[data_50$Group]

data_50$NaturalConnectivityRetention <- data_50$NaturalConnectivity/baseline_natural_connectivity[data_50$Group]

comparison_50 <- compare_50_percent(data_50)

write.csv(relabel_group_for_output(removal_raw), file.path(OUTPUT_DIR, "Node_Removal_Stability_Raw.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(relabel_group_for_output(removal_summary), file.path(OUTPUT_DIR, "Node_Removal_Stability_Summary.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(relabel_group_for_output(data_50), file.path(OUTPUT_DIR, "Node_Removal_50pct_Robustness.csv"),
  row.names = FALSE, quote = TRUE, na = "")

write.csv(comparison_50, file.path(OUTPUT_DIR, "Node_Removal_50pct_Group_Comparison.csv"), row.names = FALSE,
  quote = TRUE, na = "")

parameter_table <- data.frame(Parameter = c("Network design", "Groups", "Samples per group", "ASV prevalence threshold",
  "Minimum detections per 60-sample group", "ASV selection within each group", "Common-ASV restriction",
  "Top-N abundance restriction", "Association", "Correlation threshold", "Multiple-testing correction",
  "Significance threshold", "Edge signs retained", "Community detection", "Node-removal proportions",
  "Node-removal trajectory design", "Random repetitions per proportion", "Random seed", "Natural connectivity adjacency"),
  Value = c("One independent static network per group; no sliding window", "USF, ConF", paste(group_display[filter_audit$Group],
    filter_audit$Samples, sep = "=", collapse = "; "), PREVALENCE_CUTOFF, paste(group_display[filter_audit$Group],
    filter_audit$Minimum_positive_samples, sep = "=", collapse = "; "), "Group-specific prevalence filtering",
    "No", "No", "Spearman rank correlation", paste0("|rho| >= ", COR_R_CUTOFF), "Benjamini-Hochberg across all ASV pairs within each group",
    paste0("BH q < ", COR_Q_CUTOFF), "Positive and negative", "Fast-greedy modularity using |rho| as edge weight",
    paste(sort(unique(removal_raw$RemovalPercent)), collapse = ", "), "Nested sequential removal along one fixed random node order per repetition",
    REMOVAL_REPS, RANDOM_SEED, "Unweighted binary adjacency, matching the paper-style stability metric"),
  stringsAsFactors = FALSE)

write.csv(parameter_table, file.path(OUTPUT_DIR, "Network_Analysis_Parameters.csv"), row.names = FALSE,
  quote = TRUE, na = "")

plot_stability_curves(removal_summary)

plot_50_percent_boxplot(data_50)

writeLines(c("Static ConF-versus-USF ASV co-occurrence network analysis", "", "Design: one static network from all samples in each group; no sliding windows.",
  paste0("Prevalence: ASV detected in at least ", PREVALENCE_CUTOFF * 100, "% of group samples."),
  paste0("Edges: Spearman |rho| >= ", COR_R_CUTOFF, " and BH q < ", COR_Q_CUTOFF, "."), "Positive and negative correlations were retained.",
  "ASVs were selected independently within each group; common ASVs were not required.", "Node-removal stability used 100 nested random trajectories, recorded every 1% removal by default.",
  "Global network metrics are one descriptive value per static network and are not biological replicates.",
  "Node-level metric tables and node-removal replicates can be used for boxplots, with dependence limitations stated in interpretation."),
  con = file.path(OUTPUT_DIR, "README_Method.txt"))

cat("\nAnalysis completed: ", format(Sys.time()), "\n", sep = "")

cat("Results written to ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
