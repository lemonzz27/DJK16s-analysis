rm(list = ls())

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("DJK16S_WORK_ROOT", unset = Sys.getenv("DJK_PROJECT_ROOT", unset = ".")),
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

coverage <- rbind(prepare_variant(file.path(root, "taxonomysilva"), "SILVA"), prepare_variant(file.path(root,
  "taxonomygreen"), "Greengenes"))

write.csv(coverage, file.path(root, "annotation_results", "taxonomy_coverage_comparison.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

message("Prepared SILVA and Greengenes taxonomy variants.")
