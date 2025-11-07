# ==============================================================
# Canonical G4 Motif Detection and Stem/Loop Subsutitution Frequencies
# ==============================================================
# Author: Angelika Lahnsteiner, Dr.
# Date: 2025-10-26
# Purpose:
#   Identify canonical G4 motifs from reference genome regions,
#   detect whether the SNVs (non-effect=other vs effect sequence)
#   occurs in stem1–4 or loop1–3, and summarize distributions.
#   Create plots and tables for all analyses.
#   Analysis performed with R version 4.5.0 (2025-04-11)
# ==============================================================

# ==============================
# 1. Load libraries
# ==============================
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Biostrings)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ==============================
# 2. Set working directory
# ==============================
setwd("path/to/working_directory")

# ==============================
# 3. Load BED file
# ==============================
file_path <- "path/to/pqsfinder/pqsfinder_G4_SNVs.bed" # load pqsfinder or G4Hunter G4-SNVs
data <- read.csv(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- data[c(1,2,5,6,7,8,9,10,18)]  # keep essential columns
bed_data <- data
colnames(bed_data)
cat("Input rows:", nrow(bed_data), "\n")

# ==============================
# 4. Reverse-complement sequences if C-rich (C > G)
# ==============================
rc_if_C_rich <- function(seq_vec) {
  c_count <- str_count(seq_vec, "C")
  g_count <- str_count(seq_vec, "G")
  to_rc <- c_count > g_count
  seq_vec[to_rc] <- as.character(reverseComplement(DNAStringSet(seq_vec[to_rc])))
  return(seq_vec)
}

bed_data <- bed_data %>%
  mutate(
    otherSequence_comp  = rc_if_C_rich(other_sequence),
    effectSequence_comp = rc_if_C_rich(effect_sequence)
  )

bed_data <- bed_data[c(1,2,3,4,10,11)]
colnames(bed_data)[5] <- "otherSequence" # non_effect_sequence
colnames(bed_data)[6] <- "effectSequence"
colnames(bed_data)

# ==============================
# 5. Define canonical motif pattern
# ==============================
canonical_pattern <- paste0(
  "(G{3,})",            # stem1
  "([ATCG]{1,12})",     # loop1
  "(G{3,})",            # stem2
  "([ATCG]{1,12})",     # loop2
  "(G{3,})",            # stem3
  "([ATCG]{1,12})",     # loop3
  "(G{3,})"             # stem4
)

# ==============================
# 6. Identify which sequence contains canonical motif
# ==============================
bed_data <- bed_data %>%
  mutate(
    has_motif_other = str_detect(otherSequence, canonical_pattern),
    has_motif_effect = str_detect(effectSequence, canonical_pattern),
    motif_source = case_when(
      has_motif_other & has_motif_effect ~ "both",
      has_motif_other & !has_motif_effect ~ "only_otherSequence",
      !has_motif_other & has_motif_effect ~ "only_effectSequence",
      TRUE ~ NA_character_
    )
  )

cat("Motif origin summary:\n")
motif_table <- table(bed_data$motif_source, useNA = "ifany")
print(motif_table)

# Helper to safely extract counts
get_count <- function(x, name) ifelse(name %in% names(x), x[name], 0)
n_both    <- get_count(motif_table, "both")
n_ref     <- get_count(motif_table, "only_otherSequence")
n_var     <- get_count(motif_table, "only_effectSequence")
n_nomotif <- get_count(motif_table, NA) + get_count(motif_table, "<NA>")
total_all <- sum(motif_table)
total_detected  <- n_both + n_ref + n_var
fraction_detect <- round(total_detected / total_all * 100, 2)

cat("\n--- Summary counts ---\n")
cat("  Total entries:", total_all, "\n")
cat("  With canonical motif (any allele):", total_detected, "\n")
cat("     • Found in both alleles:", n_both, "\n")
cat("     • OtherSequence only (motif-lost):", n_ref, "\n")
cat("     • EffectSequence only (motif-gained):", n_var, "\n")
cat("  Fraction canonical motif:", fraction_detect, "%\n")

# Keep only canonical motifs
G4s_canonical <- bed_data %>% filter(!is.na(motif_source))

# ==============================
# 7. Classify mutation region
# ==============================
get_region_positions <- function(seq, match_start, groups) {
  groups <- as.numeric(groups)
  groups <- groups[!is.na(groups)]
  starts <- match_start + cumsum(c(0, head(groups, -1)))
  ends <- starts + groups - 1
  data.frame(start = starts, end = ends)
}


classify_mutation_region <- function(other_seq, effect_seq, motif_source) {
  if (is.na(motif_source)) {
    return(data.frame(region_type = NA, region_number = NA, region_pos_within = NA))
  }
  
  if (motif_source == "only_effectSequence") {
    ref_seq <- effect_seq
    alt_seq <- other_seq
  } else {
    ref_seq <- other_seq
    alt_seq <- effect_seq
  }
  
  match_info <- stringr::str_locate(ref_seq, canonical_pattern)
  if (any(is.na(match_info))) {
    return(data.frame(region_type = NA, region_number = NA, region_pos_within = NA))
  }
  
  match_start <- match_info[1]
  match_end   <- match_info[2]
  motif_seq <- substr(ref_seq, match_start, match_end)
  m <- stringr::str_match(motif_seq, canonical_pattern)
  if (all(is.na(m))) return(data.frame(region_type = NA, region_number = NA, region_pos_within = NA))
  
  groups <- nchar(m[2:8])
  region_pos <- get_region_positions(ref_seq, match_start, groups)
  stem_idxs <- c(1, 3, 5, 7)
  loop_idxs <- c(2, 4, 6)
  
  ref_chars <- strsplit(ref_seq, "")[[1]]
  alt_chars <- strsplit(alt_seq, "")[[1]]
  if (length(ref_chars) != length(alt_chars))
    return(data.frame(region_type = NA, region_number = NA, region_pos_within = NA))
  
  diffs <- which(ref_chars != alt_chars)
  if (length(diffs) == 0)
    return(data.frame(region_type = NA, region_number = NA, region_pos_within = NA))
  
  for (pos in diffs) {
    for (i in seq_along(stem_idxs)) {
      idx <- stem_idxs[i]
      if (pos >= region_pos$start[idx] && pos <= region_pos$end[idx]) {
        region_pos_within <- pos - region_pos$start[idx] + 1
        return(data.frame(region_type = "Stem", region_number = i, region_pos_within = region_pos_within))
      }
    }
    for (i in seq_along(loop_idxs)) {
      idx <- loop_idxs[i]
      if (pos >= region_pos$start[idx] && pos <= region_pos$end[idx]) {
        region_pos_within <- pos - region_pos$start[idx] + 1
        return(data.frame(region_type = "Loop", region_number = i, region_pos_within = region_pos_within))
      }
    }
  }
  
  if (any(diffs < match_start | diffs > match_end)) {
    return(data.frame(region_type = "Outside", region_number = NA, region_pos_within = NA))
  }
  
  data.frame(region_type = NA, region_number = NA, region_pos_within = NA)
}

# ==============================
# 8. Classify mutations
# ==============================
cat("Classifying mutations ...\n")

mutation_annotations <- mapply(
  function(other, effect, source) {
    classify_mutation_region(other, effect, source)
  },
  G4s_canonical$otherSequence,
  G4s_canonical$effectSequence,
  G4s_canonical$motif_source,
  SIMPLIFY = FALSE
)

mutation_df <- do.call(rbind, mutation_annotations)
G4s_canonical <- cbind(G4s_canonical, mutation_df)

cat("Mutation classification complete.\n")
print(table(G4s_canonical$region_type, useNA = "ifany"))


# some cases are located outside of our defined canonical motif e.g CGGGGCAGGGCACAAGGAGGGAGAGGGGAG compared to GGGGGCAGGGCACAAGGAGGGAGAGGGGAG
# although a canonical motif is detected by our algorithm, the SNP is located "outside". G4hunter and pqsfinder would predict such motifs as buldged structures
# which we are not including here.

# remove outside SNVs
G4s_canonical <- G4s_canonical[!(G4s_canonical$region_type=="Outside"),]

cat("\n--- Region type summary ---\n")
print(table(G4s_canonical$region_type, useNA = "ifany"))

write.table(G4s_canonical, "G4s_canonical_final.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



# ==============================
# 9. Extract stem/loop substitutions
# ==============================

extract_stem_loop_regions <- function(seq) {
  m <- stringr::str_match(seq, canonical_pattern)
  if (all(is.na(m))) {
    return(data.frame(stem1 = NA, loop1 = NA,
                      stem2 = NA, loop2 = NA,
                      stem3 = NA, loop3 = NA,
                      stem4 = NA))
  }
  data.frame(
    stem1 = m[2], loop1 = m[3],
    stem2 = m[4], loop2 = m[5],
    stem3 = m[6], loop3 = m[7],
    stem4 = m[8],
    stringsAsFactors = FALSE
  )
}

# Extract stems/loops from the sequence containing the canonical motif
stem_loop_data <- do.call(rbind, lapply(seq_len(nrow(G4s_canonical)), function(i) {
  seq_to_use <- switch(G4s_canonical$motif_source[i],
                       "only_otherSequence"  = G4s_canonical$otherSequence[i],
                       "only_effectSequence" = G4s_canonical$effectSequence[i],
                       "both"                = G4s_canonical$otherSequence[i],
                       NA)
  if (is.na(seq_to_use)) return(NULL)
  extract_stem_loop_regions(seq_to_use)
}))

# Combine
G4s_canonical <- cbind(G4s_canonical, stem_loop_data)


# ==============================
# 10. Calculate substitution frequencies and Poisson significance
# ==============================
length_cols_stem <- c("stem1", "stem2", "stem3", "stem4")
length_cols_loop <- c("loop1", "loop2", "loop3")


sum_nchar <- function(x) sum(nchar(x), na.rm = TRUE)

total_stem_bp <- sum(sapply(G4s_canonical[, length_cols_stem], nchar), na.rm = TRUE)
total_loop_bp <- sum(sapply(G4s_canonical[, length_cols_loop], nchar), na.rm = TRUE)

n_stem_mut <- sum(G4s_canonical$region_type == "Stem", na.rm = TRUE)
n_loop_mut <- sum(G4s_canonical$region_type == "Loop", na.rm = TRUE)

stem_rate_perkb <- (n_stem_mut / total_stem_bp) * 1000
loop_rate_perkb <- (n_loop_mut / total_loop_bp) * 1000

poisson_result <- poisson.test(c(n_stem_mut, n_loop_mut),
                               T = c(total_stem_bp, total_loop_bp),
                               alternative = "two.sided")

# Calculate mean lengths (in bp)
# Count number of individual stems and loops (non-NA entries)
n_stems <- sum(!is.na(unlist(G4s_canonical[, length_cols_stem])))
n_loops <- sum(!is.na(unlist(G4s_canonical[, length_cols_loop])))
mean_stem_length <- total_stem_bp / n_stems
mean_loop_length <- total_loop_bp / n_loops

summary_table <- data.frame(
  Region = c("Stem", "Loop"),
  Mutations = c(n_stem_mut, n_loop_mut),
  Analyzed_bp = c(total_stem_bp, total_loop_bp),
  Mutations_per_1kb = c(stem_rate_perkb, loop_rate_perkb),
  p_value = poisson_result$p.value,
  average_length = c(mean_stem_length, mean_loop_length)
)


# Plot: stems vs loops counts
pdf("G4_SNV_counts_stem_vs_loop.pdf", width = 9, height = 6)
ggplot(summary_table, aes(x = Region, y = Mutations, fill = Region)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = c("Stem" = "#6EA0C7", "Loop" = "#F4A261")) +
  labs(
    title = "Substitution Frequency Comparison: Stems vs Loops",
    x = "Region Type", y = "Counts"
  ) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()


# Plot: stems vs loops normalized
pdf("G4_SNV_frequency_stem_vs_loop.pdf", width = 9, height = 6)
ggplot(summary_table, aes(x = Region, y = Mutations_per_1kb, fill = Region)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = c("Stem" = "#6EA0C7", "Loop" = "#F4A261")) +
  labs(
    title = "Substitution Frequency Comparison: Stems vs Loops",
    subtitle = "Normalized per 1,000 bp analyzed",
    x = "Region Type", y = "Mutations per 1,000 bp"
  ) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()

# ==============================================================
# 11. Per-stem and per-loop substitution frequencies
# ==============================================================

stem_bp <- sapply(G4s_canonical[, length_cols_stem], sum_nchar)
loop_bp <- sapply(G4s_canonical[, length_cols_loop], sum_nchar)
stem_mut <- sapply(1:4, function(i)
  sum(G4s_canonical$region_type == "Stem" & G4s_canonical$region_number == i, na.rm = TRUE))
loop_mut <- sapply(1:3, function(i)
  sum(G4s_canonical$region_type== "Loop" & G4s_canonical$region_number == i, na.rm = TRUE))

per_region_df <- data.frame(
  Region = c(paste0("Stem", 1:4), paste0("Loop", 1:3)),
  Mutations = c(stem_mut, loop_mut),
  Analyzed_bp = c(stem_bp, loop_bp)
) %>%
  mutate(
    Mut_per_bp = Mutations / Analyzed_bp,
    Mut_per_kb = Mut_per_bp * 1000,
    p_value = sapply(1:n(), function(i)
      poisson.test(
        Mutations[i], T = Analyzed_bp[i],
        r = (n_stem_mut + n_loop_mut) / (total_stem_bp + total_loop_bp)
      )$p.value),
    padj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    padj_label = formatC(padj, format = "e", digits = 3)
  )

write.table(per_region_df, "per_region_mutation_rates.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot per-region mutation frequency
ggsave("G4_mutation_frequency_per_region3.pdf",
       ggplot(per_region_df, aes(x = Region, y = Mut_per_kb, fill = grepl("Stem", Region))) +
         geom_col(color = "black") +
         geom_text(aes(label = sprintf("%.2f", Mut_per_kb)), vjust = -0.6, size = 4.2) +
        #geom_text(aes(label = padj_label), vjust = -2.0, size = 5.2, fontface = "bold") +
         scale_fill_manual(values = c("TRUE" = "#6EA0C7", "FALSE" = "#F4A261"),
                           labels = c("Loops", "Stems")) +
         labs(title = "Substitution Frequency per G4 Subregion",
              x = "G4 Subregion", y = "Substitutions per 1,000 bp") +
         theme_minimal(base_size = 14) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 9, height = 6, device = cairo_pdf)


# ==============================
# 12. Position-wise mutation frequency within loops and stems
# ==============================

# Count mutations per position
pos_data <- G4s_canonical %>%
  filter(region_type %in% c("Stem", "Loop"),
         !is.na(region_pos_within),
         !is.na(region_number)) %>%
  mutate(region_label = paste0(tolower(region_type), region_number))

pos_counts <- pos_data %>%
  group_by(region_label, region_pos_within) %>%
  summarise(n_mut = n(), .groups = "drop")


# Compute how many sequences include each position
region_lengths_long <- G4s_canonical %>%
  select(stem1:stem4) %>%                     # assumes columns exist
  mutate(across(everything(), nchar)) %>%     # convert sequences to lengths
  pivot_longer(cols = everything(),
               names_to = "region_label",
               values_to = "len") %>%
  filter(!is.na(len), len > 0) %>%
  mutate(region_label = tolower(region_label))

# simple function: for each region, count sequences covering each position
pos_total <- region_lengths_long %>%
  group_by(region_label) %>%
  do({
    max_len <- max(.$len)
    tibble(
      region_label = unique(.$region_label),
      region_pos_within = 1:max_len,
      total_sequences = sapply(1:max_len, function(p) sum(.$len >= p))
    )
  }) %>%
  ungroup()


# Merge & calculate frequency
pos_freq <- pos_counts %>%
  left_join(pos_total, by = c("region_label", "region_pos_within")) %>%
  mutate(freq = ifelse(total_sequences > 0, n_mut / total_sequences, NA_real_),
         region_type = ifelse(grepl("^loop", region_label), "Loop", "Stem"))


write.table(pos_freq, "pos_freq_pqs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# Split for plots
pos_loops <- pos_freq %>% filter(region_type == "Loop")
pos_stems <- pos_freq %>% filter(region_type == "Stem")

#cutoff at stem pos 4, because there are to less with more than 4 nts per stem
pos_stems <- pos_stems %>%
  dplyr::filter(region_pos_within <= 4)

# Plot 1: Loops (Loop1–3)
pdf("G4_positionwise_mutation_counts_loops.pdf", width = 9, height = 5)

ggplot(pos_loops, aes(x = region_pos_within, y = n_mut, fill = region_type)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ region_label, ncol = 3, drop = FALSE) +
  scale_x_continuous(
    limits = c(1, 12),
    breaks = seq(1, 12, 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c("Loop" = "#F4A261")) +
  labs(
    title = "Position-wise Substitutions in G4 Loops",
    subtitle = "Counts of substitutions per position within each loop (1–12 bp)",
    x = "Position within loop",
    y = "Counts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

dev.off()
cat("\nPlot saved as 'G4_positionwise_mutation_counts_loops.pdf'\n")


pdf("G4_positionwise_mutation_frequency_loops.pdf", width = 9, height = 5)
ggplot(pos_loops, aes(x = region_pos_within, y = freq, fill = region_type)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ region_label, ncol = 3, drop = FALSE) +
  scale_x_continuous(breaks = 1:12)+
  scale_fill_manual(values = c("Loop" = "#F4A261")) +
  labs(
    title = "Position-wise Substitution Frequency in G4 Loops",
    subtitle = "SNV count per position normalized by number of loops containing that position",
    x = "Position within loop",
    y = "Mutation frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
dev.off()
cat("\nPlot saved as 'G4_positionwise_mutation_frequency_loops.pdf'\n")


# Plot 2: Stems (Stem1–4)
pdf("G4_positionwise_mutation_counts_stems.pdf", width = 9, height = 5)

ggplot(pos_stems, aes(x = region_pos_within, y = n_mut, fill = region_type)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ region_label, ncol = 4, drop = FALSE) +
  scale_x_continuous(
    breaks = 1:4,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c("Stem" = "#6EA0C7")) +
  labs(
    title = "Position-wise Substitutions in G4 Stems",
    subtitle = "Counts of substitutions per position within each stem",
    x = "Position within stem",
    y = "Counts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

dev.off()
cat("\nPlot saved as 'G4_positionwise_mutation_counts_stems.pdf'\n")



pdf("G4_positionwise_mutation_frequency_stems.pdf", width = 10, height = 5)

ggplot(pos_stems, aes(x = region_pos_within, y = freq, fill = region_type)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ region_label, ncol = 4, drop = FALSE) +
  scale_x_continuous(breaks = 1:4
  ) +
  scale_fill_manual(values = c("Stem" = "#6EA0C7")) +
  labs(
    title = "Position-wise Substitution Frequency in G4 Stems",
    subtitle = "Normalized by the number of stems that contain each position",
    x = "Position within stem",
    y = "Mutation frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

dev.off()
cat("\n Plot saved as 'G4_positionwise_mutation_frequency_stems.pdf'\n")


## Loop length distribution ----
# Compute lengths for each region
length_data <- G4s_canonical %>%
  mutate(
    stem1_len = nchar(stem1),
    loop1_len = nchar(loop1),
    stem2_len = nchar(stem2),
    loop2_len = nchar(loop2),
    stem3_len = nchar(stem3),
    loop3_len = nchar(loop3),
    stem4_len = nchar(stem4)
  )

# Convert to long format
lengths_long <- length_data %>%
  select(stem1_len, stem2_len, stem3_len, stem4_len,
         loop1_len, loop2_len, loop3_len) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Region",
    values_to = "Length"
  ) %>%
  mutate(
    Region = factor(Region,
                    levels = c("stem1_len", "stem2_len", "stem3_len", "stem4_len",
                               "loop1_len", "loop2_len", "loop3_len")),
    Type = ifelse(grepl("stem", Region), "Stem", "Loop"),
    Region = gsub("_len", "", Region)
  )


# Quick summary
summary_lengths <- lengths_long %>%
  group_by(Region, Type) %>%
  summarise(
    n = n(),
    Median = median(Length, na.rm = TRUE),
    Mean = mean(Length, na.rm = TRUE),
    SD = sd(Length, na.rm = TRUE),
    Min = min(Length, na.rm = TRUE),
    Max = max(Length, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Summary of region lengths ---\n")
print(summary_lengths)

pdf("G4_stem_loop_length_distribution.pdf", width = 9, height = 6)
ggplot(lengths_long, aes(x = Region, y = Length, fill = Type)) +
  geom_boxplot(outlier.alpha = 0.4, color = "black", width = 0.6) +
  scale_fill_manual(values = c("Stem" = "#6EA0C7", "Loop" = "#F4A261")) +
  labs(
    title = "Distribution of Canonical G4 Stem and Loop Lengths",
    subtitle = "Each box shows the nucleotide length distribution per subregion",
    x = "G4 Subregion",
    y = "Length (bp)",
    fill = "Region Type"
  ) +
  ylim(0,12)+
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1)
  )

dev.off()

cat("\nPlot saved as 'G4_stem_loop_length_distribution.pdf'\n")

# Save numeric summary
write.table(summary_lengths, "G4_stem_loop_length_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSummary table saved as 'G4_stem_loop_length_summary.tsv'\n")