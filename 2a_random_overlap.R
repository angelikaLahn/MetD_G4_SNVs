# ==============================================================
# Random SNV overlap with G4 motifs
# ==============================================================
# Author: Angelika Lahnsteiner, Dr. 
# Date: 2025-01-01
# Purpose: Generate a randomized SNV dataset and overlap it  
#          with pqsfinder or G4Hunter datasets
# ==============================================================

# ==============================
# 1. Load libraries 
#  ==============================
suppressPackageStartupMessages({
library(dplyr)
library(purrr)
library(tidyr)
})

#  ==============================
# 2. Load datasets 
# ===============================

#load SNP dataset
snps <- read.delim("path/to/SNVs_hg38.bed", header = FALSE) # load the Park et al. SNV dataset
head(snps)
snps <- snps[c(1,2,3,10)]
colnames(snps) <- c("chr", "start", "end", "rs_id")
head(snps)

#load G4 SNPs pqsfinder
g4_snps <- read.delim("path/to/G4_SNVs.bed", header = FALSE) # load the G4-SNV datasets predicted by pqsfinder or G4Hunter
g4_snps <- g4_snps[c(1,2,3,10)]
colnames(g4_snps) <- c("chr", "start", "end", "rs_id")
head(g4_snps)

#load G4 motifs
g4_regions <- read.delim("path/to/G4hunter_hg38.bed", header = FALSE) # load genome-wide G4Hunter or pqsfinder predicted G4 motifs

g4_regions <- g4_regions[c(1:3)]
colnames(g4_regions) <- c("chr", "start", "end")


# ==============================
# 3. Analysis  
# ==============================

#count SNP occurence per chromosome
g4_snp_counts <- g4_snps %>%
  group_by(chr) %>%
  summarise(count = n())

print(g4_snp_counts)

# Count observed overlaps
system("bedtools intersect -u -a observed_snps.bed -b g4_regions.bed > observed_overlap.bed")
observed_overlap <- readLines("observed_overlap.bed")
observed_count <- length(observed_overlap) # this must result in the same counts as in the initial analysis 

#save data 
write.table(g4_snps[c("chr", "start", "end")], "observed_snps.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(g4_regions, "g4_regions.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


# Simulate overlaps 1000 times
set.seed(123)
n_iter <- 1000
overlap_counts <- numeric(n_iter)

# draw the same number of SNVs per chromosome as in the original SNV dataset.
# calculate p value via permutation test
for (i in 1:n_iter) {
  random_snps <- g4_snp_counts %>%
    group_by(chr) %>%
    reframe(sampled_snps = list(sample(snps$rs_id[snps$chr == chr], count, replace = FALSE))) %>%
    unnest(sampled_snps)
  
  random_snps_info <- snps %>%
    filter(rs_id %in% random_snps$sampled_snps) %>%
    select(chr, start, end)
  
  write.table(random_snps_info, "random_snps.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  system("bedtools intersect -u -a random_snps.bed -b g4_regions.bed > tmp_overlap.bed")
  
  overlap_counts[i] <- length(readLines("tmp_overlap.bed"))
}

# Summary statistics
mean_random <- mean(overlap_counts) 
mean_random 
# 64 pqsfinder
# 167.254 G4Hunter
sd_random <- sd(overlap_counts)
sd_random
#7.58143 pqsfinder
#12.51208 G4hunter

# Compute empirical p-value
p_value <- (sum(overlap_counts >= observed_count) + 1) / (n_iter + 1)
min_p <- 1 / (n_iter + 1)


if (p_value == 0) {
  cat("Empirical p-value: <", format(min_p, scientific = TRUE), "\n")
  cat("Exact empirical p-value (minimum possible with this n_iter):", min_p, "\n")
} else {
  cat("Empirical p-value:", format(p_value, scientific = TRUE), "\n")
  cat("Exact empirical p-value:", p_value, "\n")
}

# Apply Benjamini–Hochberg correction
p_adj <- p.adjust(p_value, method = "BH")



