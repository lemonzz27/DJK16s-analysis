local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  options(stringsAsFactors = FALSE)

  root <- normalizePath(Sys.getenv("DJK16S_WORK_ROOT", unset = file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA")),
    winslash = "/", mustWork = TRUE)

  ranks <- c(d = "domain", p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")

  abundance <- read.csv(file.path(root, "tpm.csv"), row.names = 1L, check.names = FALSE)

  all_ids <- rownames(abundance)

  pathogen_ids <- rownames(read.csv(file.path(root, "致病菌", "tpm_raw_counts_before_rarefaction.csv"),
    row.names = 1L, check.names = FALSE))

  parse_one <- function(taxon_string) {
    out <- setNames(rep("unclassified", length(ranks)), unname(ranks))
    parts <- trimws(strsplit(ifelse(is.na(taxon_string), "", taxon_string), ";", fixed = TRUE)[[1L]])
    for (part in parts) {
      match_obj <- regexec("^([dpcofgs])__(.*)$", part)
      pieces <- regmatches(part, match_obj)[[1L]]
      if (length(pieces) != 3L)
        next
      rank_key <- pieces[2L]
      value <- trimws(pieces[3L])
      if (!nzchar(value))
        value <- "unclassified"
      out[[ranks[[rank_key]]]] <- value
    }
    out
  }

  prepare_variant <- function(source_file, label) {
    raw <- read.delim(source_file, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE,
      check.names = FALSE, fileEncoding = "UTF-8-BOM")
    names(raw)[1:3] <- c("mags", "Taxon", "Confidence")
    raw$mags <- trimws(as.character(raw$mags))
    if (anyDuplicated(raw$mags)) {
      stop(label, " contains duplicated Feature IDs.")
    }
    parsed <- t(vapply(raw$Taxon, parse_one, setNames(character(length(ranks)), unname(ranks))))
    parsed <- as.data.frame(parsed, stringsAsFactors = FALSE)
    tax <- cbind(data.frame(mags = raw$mags, stringsAsFactors = FALSE), parsed, data.frame(confidence = suppressWarnings(as.numeric(raw$Confidence)),
      stringsAsFactors = FALSE))
    tax <- tax[match(all_ids, tax$mags), , drop = FALSE]
    missing_ids <- all_ids[is.na(tax$mags)]
    if (length(missing_ids) > 0L) {
      filler <- data.frame(mags = missing_ids, domain = "unclassified", phylum = "unclassified", class = "unclassified",
        order = "unclassified", family = "unclassified", genus = "unclassified", species = "unclassified",
        confidence = NA_real_, stringsAsFactors = FALSE)
      tax[is.na(tax$mags), ] <- filler
    }
    out_dir <- file.path(root, "annotation_results", label)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(tax, file.path(out_dir, "tax.csv"), row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
    pathogen_tax <- tax[match(pathogen_ids, tax$mags), , drop = FALSE]
    write.csv(pathogen_tax, file.path(out_dir, "pathogen_tax.csv"), row.names = FALSE, quote = FALSE,
      fileEncoding = "UTF-8")
    coverage <- data.frame(annotation = label, source_file = normalizePath(source_file, winslash = "/"),
      source_rows = nrow(raw), ASV_rows = length(all_ids), matched_ASVs = sum(all_ids %in% raw$mags),
      unmatched_ASVs = sum(!all_ids %in% raw$mags), duplicated_source_IDs = anyDuplicated(raw$mags),
      mean_confidence = mean(tax$confidence, na.rm = TRUE), classified_domain = sum(tax$domain != "unclassified"),
      classified_phylum = sum(tax$phylum != "unclassified"), classified_class = sum(tax$class != "unclassified"),
      classified_order = sum(tax$order != "unclassified"), classified_family = sum(tax$family != "unclassified"),
      classified_genus = sum(tax$genus != "unclassified"), classified_species = sum(tax$species !=
        "unclassified"), stringsAsFactors = FALSE)
    write.csv(coverage, file.path(out_dir, "taxonomy_coverage.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    coverage
  }

  coverage <- prepare_variant(file.path(root, "taxonomysilva"), "SILVA")
})

local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  rm(list = ls())

  options(stringsAsFactors = FALSE, scipen = 999)

  suppressPackageStartupMessages(library(vegan))

  project_root <- normalizePath(Sys.getenv("DJK16S_WORK_ROOT", unset = file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "99_运行环境/串行运行工作区/SILVA")),
    winslash = "/", mustWork = TRUE)

  raw_backup_path <- file.path(project_root, "tpm_raw_counts_before_rarefaction.csv")

  raw_path <- if (file.exists(raw_backup_path)) {
    raw_backup_path
  } else {
    file.path(project_root, "tpm.csv")
  }

  pathogen_backup_path <- file.path(project_root, "致病菌", "tpm_raw_counts_before_rarefaction.csv")

  pathogen_raw_path <- if (file.exists(pathogen_backup_path)) {
    pathogen_backup_path
  } else {
    file.path(project_root, "致病菌", "tpm.csv")
  }

  raw <- read.csv(raw_path, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE)

  counts <- as.matrix(raw)

  storage.mode(counts) <- "numeric"

  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("The full ASV table contains invalid values.")
  }

  if (any(abs(counts - round(counts)) > 1e-08)) {
    stop("The full ASV table is not an integer raw-count matrix.")
  }

  original_depth <- colSums(counts)

  rarefaction_depth <- 21000L

  keep_samples <- names(original_depth)[original_depth >= rarefaction_depth]

  excluded_samples <- names(original_depth)[original_depth < rarefaction_depth]

  counts_for_rarefaction <- counts[, keep_samples, drop = FALSE]

  set.seed(123)

  rarefied <- t(vegan::rrarefy(t(round(counts_for_rarefaction)), sample = rarefaction_depth))

  stopifnot(all(colSums(rarefied) == rarefaction_depth))

  full_rarefied_path <- file.path(project_root, "tpm_rarefied.csv")

  write.csv(data.frame(mags = rownames(rarefied), rarefied, check.names = FALSE), full_rarefied_path, row.names = FALSE,
    quote = FALSE, fileEncoding = "UTF-8")

  pathogen_raw <- read.csv(pathogen_raw_path, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE)

  pathogen_ids <- rownames(pathogen_raw)

  if (!all(pathogen_ids %in% rownames(rarefied))) {
    stop("Some pathogen-candidate ASVs are absent from the full ASV table.")
  }

  pathogen_rarefied <- rarefied[pathogen_ids, , drop = FALSE]

  pathogen_rarefied_path <- file.path(project_root, "致病菌", "tpm_rarefied.csv")

  write.csv(data.frame(mags = rownames(pathogen_rarefied), pathogen_rarefied, check.names = FALSE), pathogen_rarefied_path,
    row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

  pathogen_backup <- pathogen_backup_path

  if (!file.exists(pathogen_backup)) {
    file.copy(pathogen_raw_path, pathogen_backup)
  }

  full_analysis_dirs <- c("a多样性", "a多样性年际变化", "CCA,rda", "decay", "lefse", "Mantel test",
    "network", "network年际变化", "PCOA", "PCOA1年际变化", "超微细菌", "韦恩图", "物种组成",
    "物种组成年际变化")

  pathogen_analysis_dirs <- c("a多样性", "a多样性年际变化", "CCA,rda", "decay", "lefse", "Mantel test",
    "network", "network年际变化", "PCOA", "PCOA1年际变化", "韦恩图", "物种组成", "物种组成年际变化")

  for (directory in full_analysis_dirs) {
    destination <- file.path(project_root, directory, "tpm.csv")
    if (file.exists(destination)) {
      file.copy(full_rarefied_path, destination, overwrite = TRUE)
    }
  }

  for (directory in pathogen_analysis_dirs) {
    destination <- file.path(project_root, "致病菌", directory, "tpm.csv")
    if (file.exists(destination)) {
      file.copy(pathogen_rarefied_path, destination, overwrite = TRUE)
    }
  }

  file.copy(pathogen_rarefied_path, file.path(project_root, "致病菌", "tpm.csv"), overwrite = TRUE)

  audit <- data.frame(sample = names(original_depth), original_reads = as.numeric(original_depth), retained_for_rarefaction = names(original_depth) %in%
    keep_samples, exclusion_reason = ifelse(names(original_depth) %in% keep_samples, "", "original reads below 21000"),
    rarefied_reads = ifelse(names(original_depth) %in% keep_samples, rarefaction_depth, NA_integer_),
    rarefaction_depth = rarefaction_depth, random_seed = 123L, pathogen_candidate_reads_after_full_table_rarefaction = c(as.numeric(colSums(pathogen_rarefied)),
      rep(NA_real_, length(excluded_samples)))[match(names(original_depth), c(colnames(pathogen_rarefied),
      excluded_samples))], stringsAsFactors = FALSE)

  write.csv(audit, file.path(project_root, "rarefaction_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  summary <- data.frame(item = c("input_table", "input_features", "samples", "minimum_original_depth",
    "maximum_original_depth", "original_depth_P05", "original_depth_P10", "original_depth_median", "rarefaction_depth",
    "random_seed", "retained_samples", "excluded_samples", "excluded_sample_ids", "pathogen_candidate_features",
    "full_table_md5", "rarefied_table_md5"), value = c(normalizePath(raw_path, winslash = "/"), nrow(counts),
    ncol(counts), min(original_depth), max(original_depth), unname(stats::quantile(original_depth, 0.05)),
    unname(stats::quantile(original_depth, 0.1)), unname(stats::median(original_depth)), rarefaction_depth,
    123L, length(keep_samples), length(excluded_samples), paste(excluded_samples, collapse = ";"), nrow(pathogen_rarefied),
    unname(tools::md5sum(raw_path)), unname(tools::md5sum(full_rarefied_path))), stringsAsFactors = FALSE)

  write.csv(summary, file.path(project_root, "rarefaction_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  png(file.path(project_root, "rarefaction_depth_distribution.png"), width = 2100, height = 1500, res = 220)

  hist(original_depth, breaks = 18, col = "#9ECAE1", border = "white", main = "ASV sequencing-depth distribution",
    xlab = "Reads per sample", ylab = "Number of samples")

  abline(v = rarefaction_depth, col = "#D7301F", lwd = 3, lty = 2)

  abline(v = median(original_depth), col = "#1B7837", lwd = 3, lty = 3)

  legend("topright", legend = c(paste0("Rarefaction depth = ", rarefaction_depth), paste0("Median = ",
    median(original_depth)), paste0("Excluded = ", length(excluded_samples), "/", length(original_depth))),
    col = c("#D7301F", "#1B7837", "black"), lty = c(2, 3, 0), lwd = c(3, 3, 0), pch = c(NA, NA, 16),
    bty = "n")

  dev.off()

  message("Rarefaction completed at ", rarefaction_depth, " reads/sample; retained ", length(keep_samples),
    "/", length(original_depth), " samples and distributed tables.")
})
