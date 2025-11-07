# ==============================================================
# Find an overlap between G4-SNVs and alternative promoters
# ==============================================================
# Author: Angelika Lahnsteiner, Dr.
# Date: 2025-02-11
# ==============================================================

# ==============================
# 1. Load libraries 
# ==============================
suppressPackageStartupMessages({
library(regioneR)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg38")
library(annotatr)
library(dplyr)
})


# ==============================
# 2. Load data
# ==============================

# altProm file
altprom <- read.delim("path_to_your_file/altprom.txt", header = TRUE) #load Demircioglu altprom file

#load G4-SNPs files
df <- read.delim("path_to_your_file/G4_SNVs.txt", header = TRUE) #load pqsfinder or G4Hunter G4-SNV files

#prepare GRanges object for alternative promoters
gr_altProm_all <- GRanges(seqnames = altProm$seqnames, 
                          ranges = IRanges(start = altProm$start-250, end = altProm$start+250),
                          promoterID=altProm$promoterId)


#prepare GRanges object for SNVs
gr_SNPs <- GRanges(seqnames = df$chr,
                   ranges = IRanges(start = df$start_G4, end = df$end_G4),
                   rs_id=df$rs_id)



# ==============================
# 3. Find Overlap with G4s
# ==============================
#all regions
overlaps <- findOverlaps(gr_altProm_all,gr_SNPs)
# Extract the overlapping positions
altProm.G4.SNP.overlap <- as.data.frame(gr_altProm_all[queryHits(overlaps)]) 
# Extract corresponding SNPs (subjects)
snp_hits <- gr_SNPs[subjectHits(overlaps)]

altProm.G4.SNP.overlap$rs_id <- mcols(snp_hits)$rs_id

# Remove duplicates based on promoterID
altProm.G4.SNP.overlap <- altProm.G4.SNP.overlap %>%
  distinct(promoterID, .keep_all = T)
