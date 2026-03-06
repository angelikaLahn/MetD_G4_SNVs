###############################################################################
# Transcription Factor Motif Analysis in G4 SNP Regions
# Author: Angelika Lahnsteiner, Dr.
# Date: 2025-10-26
#
# Description:
#   - Load SNP-associated G4 regions (PQSfinder or G4Hunter)
#   - Extract both alleles and scan for TF motifs (JASPAR2024)
#   - Quantify, normalize, and compare TF motif occurrences
#   - Identify significantly differentially bound TFs (Fisher test)
#   - Prepare plots and tables
###############################################################################

# ==============================
# 1. Load libraries 
# ==============================
suppressPackageStartupMessages({
  library(JASPAR2024)
  library(TFBSTools)
  library(motifmatchr)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(tidyverse)
  library(scales)
  library(SummarizedExperiment)
})

# ==============================
# 2. Load datasets
# ==============================
df <- read.delim("path_to_your_file/G4_SNPs.bed", header = TRUE) #pqsfinder or G4Hunter predicted G4-SNVs

#run the code once with other and once with effect-alleles:
g4.snps <- df[c(2,5,6,10,1)] #other allele = non-effect
g4.snps <- df[c(2,5,6,9,1)] #effect allele

colnames(g4.snps)[4] <- "sequence"

# ==============================
# 3. Build PFM Matrix and get hits
# ==============================

# Create DNAStringSet
g4.snp.sequ <- DNAStringSet(g4.snps$sequence)

# add names to identify each entry
names(g4.snp.sequ) <- paste0(g4.snps$rs_id)

#JASPAR 2024
JASPAR2024 <- JASPAR2024()
db(JASPAR2024)
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

# 1. Pull metadata once
human_meta <- dbGetQuery(JASPARConnect, "
  SELECT M.ID, M.BASE_ID, M.VERSION, M.NAME, M.COLLECTION AS matrixClass
  FROM MATRIX M
  JOIN MATRIX_SPECIES MS ON M.ID = MS.ID
  WHERE MS.TAX_ID = 9606
")

buildPFM <- function(motif_id, conn, meta_df) {
  # a. Pull the long-format matrix data
  long_df <- dbGetQuery(conn, sprintf("
    SELECT row, col, val
    FROM MATRIX_DATA
    WHERE ID = %d
    ORDER BY row, col
  ", motif_id))
  
  # b. Pivot to wide so each 'col' becomes a column
  wide_df <- long_df %>%
    pivot_wider(names_from = col, values_from = val) %>%
    arrange(row)
  
  # c. Extract position columns and coerce to numeric matrix
  pos_cols <- setdiff(names(wide_df), "row")
  mat_vals <- as.numeric(unlist(wide_df[ , pos_cols, drop = FALSE]))
  pfm_mat <- matrix(mat_vals, nrow = 4, byrow = FALSE)
  
  # d. Assign correct rownames and colnames
  rownames(pfm_mat) <- c("A","C","G","T")
  colnames(pfm_mat) <- pos_cols
  
  # e. Fetch metadata for this motif
  md <- meta_df %>% filter(ID == motif_id)
  
  # f. Build and return the PFMatrix
  PFMatrix(
    ID            = paste0(md$BASE_ID, ".", md$VERSION),
    name          = md$NAME,
    matrixClass   = md$matrixClass,
    strand        = "+",
    bg            = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    profileMatrix = pfm_mat
  )
}

# 2. Build list of PFMatrix objects
pfm_list <- lapply(human_meta$ID, buildPFM, conn = JASPARConnect, meta_df = human_meta)

# 3. Combine into a PFMatrixList
all(sapply(pfm_list, inherits, "PFMatrix"))
pfm_list <- do.call(PFMatrixList, pfm_list)

# Quick check
length(pfm_list)
pfm_list[[1]]

# full PWM list
all_pwms <- toPWM(pfm_list)

hits <- matchMotifs(
  pwms         = all_pwms,
  subject      = g4.snp.sequ,
  genome       = hg38
)

motif_matrix <- assay(hits, "motifMatches")  # Logical matrix: rows = SNP sequences, columns = motifs
seqs_per_motif <- colSums(motif_matrix)    
motif_names <- colData(hits)$name
matched_motifs <- motif_names[seqs_per_motif > 0]
motif_match_summary <- data.frame(
  motif = motif_names,
  matches = seqs_per_motif
)

# Only motifs with at least one match
motif_match_summary <- subset(motif_match_summary, matches > 0)
head(motif_match_summary)

motif_match_summary <- motif_match_summary %>%
  group_by(motif) %>%
  summarise(matches = sum(matches, na.rm = TRUE)) %>%
  arrange(desc(matches))



# ==============================
# 3. Rename datasets
# ==============================

# prepare two dataframes from the results
df.other <- motif_match_summary
df.other$allele <- "Non-effect allele"

df.effect <- motif_match_summary
df.effect$allele <- "Effect allele"

#combine dataframes
df_combined <- rbind(df.effect, df.other)

#normalize 
df_combined <- df_combined %>%
  group_by(allele) %>%
  mutate(prop_sites = matches / sum(matches) * 100) %>%
  ungroup()


# ==============================
# 3. Prepare plots of all TFBS hits
# ==============================

# Select top motifs (e.g., top 20 across both alleles)
top_motifs <- df_combined %>%
  group_by(motif) %>%
  summarise(total = sum(matches)) %>%
  slice_max(total, n = 20) %>%
  pull(motif)

df_top <- df_combined %>%
  filter(motif %in% top_motifs)

df_summary <- df_top %>%
  group_by(motif, allele) %>%
  summarise(
    total = sum(matches),
    mean = mean(matches),
    sd = sd(matches),
    n = n(),
    .groups = "drop"
  )

# Plot
p1 <-  ggplot(df_top, aes(x = fct_reorder(motif, matches, .fun = sum), y = matches,fill = allele)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "grey20", width = 0.7) +
  coord_flip() +
  labs(
    title = "G4Hunter TF Counts — Effect vs. Non-Effect Allele",
    x = "Motif",
    y = "Total binding sites",
    fill = "Allele"
  ) +
  ylim(0,3000) +
  scale_fill_manual(values = c("Effect allele" = "violetred4", "Other allele" = "steelblue")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "grey50", fill = NA, size = 0.8),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) 

p1



# ==============================
# 4. Extract the TFBS hits with the greatest difference
# ==============================

# find the greatest difference
df.diff <- merge(df.effect, df.other, by = "motif", suffixes = c(".effect", ".other"))

df.diff <- df.diff %>%
  dplyr::mutate(
    diff = matches.effect - matches.other,           # signed difference
    abs_diff = abs(matches.effect - matches.other)   # absolute difference
  )

# By absolute difference
df.diff_sorted <- df.diff %>%
  dplyr::arrange(desc(abs_diff))

# Top 10 TFs with greatest difference
top_TFs <- head(df.diff_sorted, 10)
top_TFs


# Prepare colors for positive/negative differences
top_TFs$direction <- ifelse(top_TFs$diff > 0, "Effect > Non-effect", "Non-effect > Effect")


p3 <- ggplot(top_TFs, aes(
  x = fct_reorder(motif, abs_diff),
  y = diff,
  fill = direction
)) +
  geom_col(color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # zero line
  ylim(-200,100)+
  coord_flip() +
  scale_fill_manual(values = c("Effect > Non-effect" = "violetred4", "Non-effect > Effect" = "steelblue")) +
  labs(title = "Top TFs with Greatest Difference",
       x = "Transcription Factor", y = "Difference", fill = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) 

p3


# ==============================
# 5. Fishers exact test for significance
# ==============================

# Totals for each group
total_effect_TF  <- sum(df.diff$matches.effect)    # 64568
total_effect_TF  <- sum(df.diff$matches.other)   # 65260


# Fisher's exact test for each pairwise comparison
df.diff$effect_vs_other <- mapply(function(a, b) {
  mat <- matrix(
    c(a, total_effect_TF - a,
      b, total_effect_TF - b),
    nrow = 2
  )
  fisher.test(mat)$p.value
}, df.diff$matches.effect, df.diff$matches.other)


#  adjust p-values for multiple testing (FDR)
df.diff$effect_vs_other_FDR <- p.adjust(df.diff$effect_vs_other, method = "fdr")
# View results
df.diff

# Function to convert p-value to stars
p_to_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Add significance stars for FDR-adjusted p-values instead
df.diff$sign <- sapply(df.diff$effect_vs_other_FDR, p_to_stars)

# View final table
df.diff

write.table(df.diff,"path_to_your_folder/TFBScounts.csv",sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
