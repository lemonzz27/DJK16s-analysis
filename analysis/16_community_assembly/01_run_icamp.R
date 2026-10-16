#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

user_r_library <- Sys.getenv("ICAMP_R_LIBRARY", unset = "")
if (nzchar(user_r_library) && dir.exists(user_r_library)) {
  .libPaths(unique(c(user_r_library, .libPaths())))
}

required_packages <- c("iCAMP", "ape")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nCurrent library paths:\n", paste(.libPaths(), collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(iCAMP)
  library(ape)
})

base_dir <- Sys.getenv("ASV_EXPORT_DIR", unset = file.path("data", "qiime2", "06_export"))
table_file <- Sys.getenv("ASV_TABLE_FILE", unset = file.path(base_dir, "ASV_table.tsv"))
tree_file <- Sys.getenv(
  "TREE_FILE",
  unset = file.path(base_dir, "phylogeny", "rooted_tree", "tree.nwk")
)

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)

if (length(file_arg) > 0) {
  script_file <- sub("^--file=", "", file_arg[1])
  code_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
} else {
  code_dir <- normalizePath(getwd(), mustWork = TRUE)
}

output_dir <- Sys.getenv(
  "ICAMP_OUTPUT_DIR",
  unset = file.path("results", "icamp_bin48_rand1000")
)
pd_dir <- file.path(output_dir, "pd_big")

if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Output directory is not empty. Choose a new ICAMP_OUTPUT_DIR: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pd_dir, recursive = TRUE, showWarnings = FALSE)

nworker <- 16L
rand_time <- 1000L
bin_size_limit <- 48L
memory_gb <- 50
random_seed <- 123L

cat("iCAMP analysis started:", format(Sys.time()), "\n")
cat("R version:", R.version.string, "\n")
cat("iCAMP version:", as.character(packageVersion("iCAMP")), "\n")
cat("R libraries:\n", paste(.libPaths(), collapse = "\n"), "\n")
cat("ASV table:", table_file, "\n")
cat("Tree:", tree_file, "\n")
cat("Code directory:", code_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Workers:", nworker, "\n")
cat("Randomizations:", rand_time, "\n")
cat("bin.size.limit:", bin_size_limit, "\n")

if (!file.exists(table_file)) {
  stop("ASV table does not exist: ", table_file)
}

if (!file.exists(tree_file)) {
  stop("Phylogenetic tree does not exist: ", tree_file)
}

cat("Reading the phylogenetic tree...\n")
tree <- ape::read.tree(tree_file)

if (is.null(tree$tip.label) || length(tree$tip.label) == 0) {
  stop("The phylogenetic tree has no tip labels.")
}

if (anyDuplicated(tree$tip.label)) {
  duplicate_tips <- unique(tree$tip.label[duplicated(tree$tip.label)])
  writeLines(duplicate_tips, file.path(output_dir, "duplicated_tree_tips.txt"))
  stop("Duplicated tree tip labels were detected.")
}

if (is.null(tree$edge.length)) {
  stop("The phylogenetic tree does not contain branch lengths.")
}

if (any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
  stop("The tree contains invalid or negative branch lengths.")
}

tree_tips_before <- length(tree$tip.label)
tree_rooted_before <- ape::is.rooted(tree)

cat("Tree tips:", tree_tips_before, "\n")
cat("Tree is rooted:", tree_rooted_before, "\n\n")

cat("Reading the ASV table...\n")
first_line <- readLines(table_file, n = 1, warn = FALSE, encoding = "UTF-8")

skip_lines <- if (
  length(first_line) > 0 &&
  grepl("^# Constructed from biom file", first_line)
) 1L else 0L

asv_raw <- read.delim(
  table_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  skip = skip_lines,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

if (nrow(asv_raw) == 0 || ncol(asv_raw) == 0) {
  stop("The ASV table is empty.")
}

if (anyDuplicated(rownames(asv_raw))) {
  duplicate_ids <- unique(rownames(asv_raw)[duplicated(rownames(asv_raw))])
  writeLines(duplicate_ids, file.path(output_dir, "duplicated_table_row_ids.txt"))
  stop("Duplicated IDs were detected in the first column of the ASV table.")
}

if (anyDuplicated(colnames(asv_raw))) {
  duplicate_ids <- unique(colnames(asv_raw)[duplicated(colnames(asv_raw))])
  writeLines(duplicate_ids, file.path(output_dir, "duplicated_table_column_ids.txt"))
  stop("Duplicated column names were detected in the ASV table.")
}

row_overlap <- sum(rownames(asv_raw) %in% tree$tip.label)
column_overlap <- sum(colnames(asv_raw) %in% tree$tip.label)

cat("Tree tips matching table rows:", row_overlap, "\n")
cat("Tree tips matching table columns:", column_overlap, "\n")

if (row_overlap == 0 && column_overlap == 0) {
  stop(
    "No ASV IDs match between the abundance table and the tree. ",
    "Confirm that the files came from the same QIIME2 run."
  )
}

asv_matrix <- as.matrix(asv_raw)
suppressWarnings(storage.mode(asv_matrix) <- "numeric")

if (anyNA(asv_matrix) || any(!is.finite(asv_matrix))) {
  stop(
    "The abundance table contains non-numeric or missing values. ",
    "Use ASV_table.tsv, not ASV_table_with_taxonomy.tsv."
  )
}

if (any(asv_matrix < 0)) {
  stop("The abundance table contains negative values.")
}

if (row_overlap >= column_overlap) {
  comm <- t(asv_matrix)
  cat("Detected ASVs in rows; the table was transposed.\n")
} else {
  comm <- asv_matrix
  cat("Detected ASVs in columns; no transpose was needed.\n")
}

storage.mode(comm) <- "numeric"
samples_before <- nrow(comm)
taxa_before <- ncol(comm)

zero_samples <- rownames(comm)[rowSums(comm) == 0]
zero_taxa <- colnames(comm)[colSums(comm) == 0]

writeLines(zero_samples, file.path(output_dir, "zero_samples_removed.txt"))
writeLines(zero_taxa, file.path(output_dir, "zero_ASVs_removed.txt"))

comm <- comm[rowSums(comm) > 0, , drop = FALSE]
comm <- comm[, colSums(comm) > 0, drop = FALSE]

if (nrow(comm) < 2) {
  stop("Fewer than two non-zero samples remain after filtering.")
}

common_taxa <- intersect(colnames(comm), tree$tip.label)
table_only_taxa <- setdiff(colnames(comm), tree$tip.label)
tree_only_taxa <- setdiff(tree$tip.label, colnames(comm))

writeLines(table_only_taxa, file.path(output_dir, "ASV_only_in_table.txt"))
writeLines(tree_only_taxa, file.path(output_dir, "ASV_only_in_tree.txt"))

if (length(common_taxa) < bin_size_limit) {
  stop(
    "Only ", length(common_taxa),
    " ASVs are shared by the table and tree, which is less than ",
    "bin.size.limit = ", bin_size_limit, "."
  )
}

tree <- ape::keep.tip(tree, common_taxa)
comm <- comm[, tree$tip.label, drop = FALSE]
comm <- comm[rowSums(comm) > 0, , drop = FALSE]

if (!identical(colnames(comm), tree$tip.label)) {
  stop("The ASV order does not match the tree tip order after filtering.")
}

if (nrow(comm) < 2) {
  stop("Fewer than two samples remain after matching the tree.")
}

cat("Samples retained:", nrow(comm), "\n")
cat("ASVs retained:", ncol(comm), "\n")
cat("ASVs found only in table:", length(table_only_taxa), "\n")
cat("Tips found only in tree:", length(tree_only_taxa), "\n")
cat("Minimum sample depth:", format(min(rowSums(comm)), scientific = FALSE), "\n")
cat("Maximum sample depth:", format(max(rowSums(comm)), scientific = FALSE), "\n\n")

write.table(
  data.frame(SampleID = rownames(comm), comm, check.names = FALSE),
  file.path(output_dir, "community_used_by_icamp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ape::write.tree(tree, file.path(output_dir, "tree_used_by_icamp.nwk"))

qc <- data.frame(
  Item = c(
    "Samples_before_filtering",
    "ASVs_before_filtering",
    "Zero_samples_removed",
    "Zero_ASVs_removed",
    "Samples_in_final_analysis",
    "ASVs_in_final_analysis",
    "ASVs_only_in_table",
    "Tips_only_in_tree",
    "Minimum_sample_depth",
    "Maximum_sample_depth",
    "Input_tree_is_rooted",
    "Randomizations",
    "CPU_workers",
    "Bin_size_limit"
  ),
  Value = c(
    samples_before,
    taxa_before,
    length(zero_samples),
    length(zero_taxa),
    nrow(comm),
    ncol(comm),
    length(table_only_taxa),
    length(tree_only_taxa),
    min(rowSums(comm)),
    max(rowSums(comm)),
    tree_rooted_before,
    rand_time,
    nworker,
    bin_size_limit
  ),
  stringsAsFactors = FALSE
)

write.csv(qc, file.path(output_dir, "input_QC.csv"), row.names = FALSE)

cat("Starting iCAMP null-model analysis...\n")

set.seed(random_seed)
analysis_start <- Sys.time()

icamp_out <- iCAMP::icamp.big(
  comm = comm,
  tree = tree,
  pd.wd = pd_dir,
  output.wd = output_dir,
  prefix = "DJK16s_bin48",
  rand = rand_time,
  ds = 0.2,
  phylo.metric = "bMPD",
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  sig.index = "SES.RC",
  ses.cut = 1.96,
  rc.cut = 0.95,
  bin.size.limit = bin_size_limit,
  nworker = nworker,
  memory.G = memory_gb,
  sp.check = TRUE,
  ignore.zero = TRUE,
  correct.special = TRUE,
  detail.save = TRUE,
  detail.null = FALSE,
  qp.save = TRUE,
  rtree.save = FALSE,
  taxo.metric = "bray",
  omit.option = "no",
  unit.sum = rowSums(comm)
)

analysis_end <- Sys.time()

cat("iCAMP calculation finished.\n")
cat("Elapsed time:", as.character(analysis_end - analysis_start), "\n")

saveRDS(
  icamp_out,
  file.path(output_dir, "DJK16s_bin48_icamp_result.rds"),
  compress = FALSE
)

writeLines(
  names(icamp_out),
  file.path(output_dir, "icamp_result_object_names.txt")
)

result_name <- grep("^bNRIiRCa$", names(icamp_out), value = TRUE)
if (length(result_name) == 0) {
  result_name <- grep("bNRI.*RCa", names(icamp_out), value = TRUE)
}

if (length(result_name) > 0) {
  result_name <- result_name[1]
  pairwise_result <- as.data.frame(icamp_out[[result_name]])

  write.csv(
    pairwise_result,
    file.path(output_dir, "DJK16s_bin48_pairwise_process_importance.csv"),
    row.names = FALSE
  )

  if (ncol(pairwise_result) >= 7) {
    process_values <- pairwise_result[, 3:ncol(pairwise_result), drop = FALSE]
    process_values[] <- lapply(process_values, as.numeric)
    mean_values <- colMeans(process_values, na.rm = TRUE)

    overall_summary <- data.frame(
      Process = names(mean_values),
      Mean_relative_importance = as.numeric(mean_values),
      Mean_percentage = as.numeric(mean_values) * 100,
      stringsAsFactors = FALSE
    )

    write.csv(
      overall_summary,
      file.path(output_dir, "DJK16s_bin48_overall_process_summary.csv"),
      row.names = FALSE
    )

    cat("\nMean relative importance across all sample pairs:\n")
    print(overall_summary)
  }
} else {
  warning(
    "No bNRI/RCa result object was found. See ",
    file.path(output_dir, "icamp_result_object_names.txt")
  )
}

parameters <- data.frame(
  Parameter = c(
    "table_file", "tree_file", "output_dir", "nworker", "rand",
    "bin.size.limit", "memory.G", "random_seed", "phylo.metric",
    "phylo.rand.scale", "taxa.rand.scale", "sig.index", "ses.cut",
    "rc.cut", "taxo.metric", "ds"
  ),
  Value = c(
    table_file, tree_file, output_dir, nworker, rand_time,
    bin_size_limit, memory_gb, random_seed, "bMPD", "within.bin",
    "across.all", "SES.RC", 1.96, 0.95, "bray", 0.2
  ),
  stringsAsFactors = FALSE
)

write.csv(
  parameters,
  file.path(output_dir, "icamp_analysis_parameters.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "R_sessionInfo.txt")
)

cat("iCAMP analysis completed successfully:\n")
cat(format(Sys.time()), "\n")
cat("Results:", output_dir, "\n")
