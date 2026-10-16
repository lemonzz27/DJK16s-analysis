library(vegan)

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = "shannon", base = base)
  Simpson <- diversity(x, index = "simpson")
  Pielou <- Shannon/log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1)/rowSums(x)
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- picante::pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- "PD_whole_tree"
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}

otu <- read.csv("tpm.csv", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

otu <- t(otu)

otu <- as.data.frame(otu)

sample_depths <- rowSums(otu)

min_depth <- min(sample_depths)

set.seed(123)

otu_rarefied <- rrarefy(otu, sample = min_depth)

otu_rarefied <- as.data.frame(otu_rarefied)

alpha_all <- alpha(otu_rarefied, base = 2)

write.csv(alpha_all, "alpha1.csv", quote = FALSE)
