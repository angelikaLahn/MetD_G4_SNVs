# ==============================================================
# Random SNV overlap with G4 motifs
# ==============================================================
# Author: Angelika Lahnsteiner, Dr. 
# Date: 2026-03-06
# Purpose: Generate a randomized SNV dataset and overlap it  
#          with pqsfinder or G4Hunter datasets
# ==============================================================

# ==============================
# 1. Load libraries 
#  ==============================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
})

#  ==============================
# 2. Load datasets 
# ===============================

#load SNP dataset
snps_file <- "pqsfinder_SNVs_binned.bed"        # SNVs with GC bin in column 5
genome_bins_file <- "genome_50bp_GC_binned.bed" # genome windows with GC bin in column 5
g4_file <- "pqsfinder_hg38.bed"                 # PQS/G4 regions

snps <- import(snps_file)
genome_bins <- import(genome_bins_file)
g4_regions <- import(g4_file)

# Add GC bin to GRanges metadata
snps$gc_bin <- snps$score
genome_bins$gc_bin <- genome_bins$score

length(snps)
table(seqnames(snps))

#  ==============================
# 3. Remove chromosomes not represented in genome bins
# ===============================
snps <- snps[seqnames(snps) %in% seqlevels(genome_bins)]

#  ==============================
# 4. Observed overlaps
# ===============================
observed_count <- sum(countOverlaps(snps, g4_regions) > 0)
cat("Observed overlaps:", observed_count, "\n")

# ===============================
# 5. Indexing
# ===============================
cat("Computing SNP and genome bins...\n")

# Unique chromosomes and bins
chroms <- intersect(unique(seqnames(snps)), unique(seqnames(genome_bins)))
bins <- sort(unique(snps$gc_bin))

# ---- SNPs ----
snps_by_chr_bin <- list()

for(chr in chroms){
  snps_chr <- snps[seqnames(snps) == chr]
  snps_by_chr_bin[[as.character(chr)]] <- split(snps_chr, snps_chr$gc_bin)
}

# ---- genome bins ----
genome_by_chr_bin <- list()

for(chr in chroms){
  bins_chr <- genome_bins[seqnames(genome_bins) == chr]
  genome_by_chr_bin[[as.character(chr)]] <- split(bins_chr, bins_chr$gc_bin)
}

cat("Computation done.\n")


# ===============================
# 6. Permutation function
# ===============================
permute <- function(snps_by_chr_bin, genome_by_chr_bin){
  
  shuffled <- GRanges()  # empty GRanges
  
  for(chr in names(snps_by_chr_bin)){
    
    snps_bins <- snps_by_chr_bin[[chr]]
    genome_bins_chr <- genome_by_chr_bin[[chr]]
    
    for(bin in names(snps_bins)){
      
      snps_bin <- snps_bins[[bin]]
      genome_bin <- genome_bins_chr[[bin]]
      
      if(length(snps_bin) == 0 || length(genome_bin) == 0)
        next
      
      sampled <- sample(genome_bin,
                        length(snps_bin),
                        replace = TRUE)
      
      new_snps <- resize(sampled, width = 1, fix = "center")
      
      shuffled <- c(shuffled, new_snps)
    }
  }
  
  shuffled
}
# ===============================
# 7. Run permutations
# ===============================
set.seed(123)
n_iter <- 1000
overlap_counts <- numeric(n_iter)

start_time <- Sys.time()

for(i in 1:n_iter){
  
  cat(sprintf("\n[%s] Iteration %d/%d\n",
              Sys.time(), i, n_iter))
  
  shuffled_snps <- permute(snps_by_chr_bin, genome_by_chr_bin)
  
  overlap_counts[i] <- sum(countOverlaps(shuffled_snps, g4_regions) > 0)
  
  elapsed <- Sys.time() - start_time
  eta <- elapsed / i * (n_iter - i)
  
  cat(sprintf("Overlaps: %d | elapsed: %.1f min | ETA: %.1f min\n",
              overlap_counts[i],
              as.numeric(elapsed, "mins"),
              as.numeric(eta, "mins")))
}

# ===============================
# 8. Compute statistics
# ===============================
mean_null <- mean(overlap_counts)
sd_null <- sd(overlap_counts)
fold_enrichment <- observed_count / mean_null
zscore <- (observed_count - mean_null) / sd_null
pval_two_sided <- mean(abs(overlap_counts - mean_null) >= abs(observed_count - mean_null))

cat("\n===== RESULTS =====\n")
cat("Mean null:", mean_null, "\n")
cat("SD null:", sd_null, "\n")
cat("Fold enrichment:", fold_enrichment, "\n")
cat("Z-score:", zscore, "\n")
cat("Two-sided p-value:", pval_two_sided, "\n")