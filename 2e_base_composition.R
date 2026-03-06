# ==============================================================
# Canonical G4 motif detection and GC content analysis
# ==============================================================
# Author: Angelika Lahnsteiner, Dr.
# Date: 26.10.2025
# Purpose: Identify canonical G4s (G>=3N1–12), 
#          calculate per-loop and per-base GC content,
#          and visualize GC distribution and significance.
# ==============================================================

# ==============================
# 1. Load libraries 
# ==============================
suppressPackageStartupMessages({
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Biostrings)
})

# ==============================
# 2. Load data
# ==============================

# Load non-SNV-G4 datasets

# Define file path
file_path <- "/usr/local/CAME/epi2/Angelika_Lahnsteiner/G4_SNP/Stem_loop_overlap/pqsfinder/Base_composition/pqs.hg38.woSNVregions.txt"
# Read the comma-separated file
bed_data <- read.csv(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
# View first few rows
head(bed_data)

## ALTERNATIVELY
# load G4Hunter or pqsfinder G4-SNV dataset and perform the same analysis
# Load the canonical motifs from the substitution analysis 
file_path <- "/usr/local/CAME/epi2/Angelika_Lahnsteiner/G4_SNP/Stem_loop_overlap/pqsfinder/substitution_frequency/NEW_9bp/G4s_canonical_final_short.tsv"
G4s_canonical_SNV_G4s <- read.csv(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(G4s_canonical_SNV_G4s)
bed_data <- G4s_canonical_SNV_G4s[c(2:5)] # extract chr, start, end, sequence columns


# ==============================
# 3. Reverse-complement sequences if C-rich (C > G)
# ==============================
# Prepare just the G-motifs, reverse complement in the case of C-rich regions
# Count C and G in all sequences at once
rc_if_C_rich <- function(seq_vec) {
  c_count <- str_count(seq_vec, "C")
  g_count <- str_count(seq_vec, "G")
  to_rc <- c_count > g_count
  seq_vec[to_rc] <- as.character(reverseComplement(DNAStringSet(seq_vec[to_rc])))
  return(seq_vec)
}

bed_data <- bed_data %>%
  mutate(
    sequence_comp  = rc_if_C_rich(sequence),
  )


bed_data <- bed_data[c(1,2,3,5)]
#bed_data <- bed_data[c(2,3,4,6)]
colnames(bed_data)[4] <- "sequence"
colnames(bed_data)

# View first few rows
head(bed_data)

# ==============================
# 4. Define canonical G4 pattern (strict canonical)
# ==============================
# Four G-runs (≥3 Gs each), three loops (1–12 nt), 
# loops may include Gs or GG but not GGG tracts.
canonical_pattern <- paste0(
  "(G{3,})",                # stem1
  "([ATCG]{1,12})",           # loop1
  "(G{3,})",                # stem2
  "([ATCG]{1,12})",           # loop2
  "(G{3,})",                # stem3
  "([ATCG]{1,12})",           # loop3
  "(G{3,})"                # stem4
)


# ==============================
# 5. Filter sequences matching canonical G4 motif
# ==============================
G4s_canonical <- bed_data[str_detect(bed_data$sequence, canonical_pattern), ]
head(G4s_canonical)
nrow(G4s_canonical) 
#for later downsampling
non_SNV_G4s<-G4s_canonical[c(1:4)]

# ==============================
# 6. Extract GC metrics 
# ==============================

#define GC content calculation
gc_content <- function(sequence) {
  sequence <- toupper(sequence)
  if (nchar(sequence) == 0) return(NA)
  gc <- str_count(sequence, "[GC]")
  return(gc / nchar(sequence))
}

analyze_G4_GC_canonical <- function(sequence) {
  sequence <- toupper(sequence)
  m <- str_match(sequence, canonical_pattern)
  
  if (all(is.na(m))) {
    return(data.frame(
      n_stems = NA, n_loops = NA,
      mean_stem_GC = NA, mean_loop_GC = NA,
      stem1_GC = NA, stem2_GC = NA, stem3_GC = NA, stem4_GC = NA,
      loop1_GC = NA, loop2_GC = NA, loop3_GC = NA,
      loop1_GC_pos = I(list(NA)), loop2_GC_pos = I(list(NA)), loop3_GC_pos = I(list(NA))
    ))
  }

  # Extract stems and loops
  stems <- m[c(2,4,6,8)]
  loops <- m[c(3,5,7)]

  # Compute GC content
  stem_gc <- sapply(stems, gc_content)
  loop_gc <- sapply(loops, gc_content)

  # Compute per-base GC (binary vector for each loop)
  loop_gc_pos <- lapply(loops, function(l) {
    nts <- strsplit(l, "")[[1]]
    as.numeric(nts %in% c("G", "C"))
  })

  # Build output
  out <- data.frame(
    n_stems = 4,
    n_loops = 3,
    mean_stem_GC = mean(stem_gc),
    mean_loop_GC = mean(loop_gc)
  )

  for (i in 1:4) out[[paste0("stem", i, "_GC")]] <- stem_gc[i]
  for (i in 1:3) out[[paste0("loop", i, "_GC")]] <- loop_gc[i]
  for (i in 1:3) out[[paste0("loop", i, "_GC_pos")]] <- I(list(loop_gc_pos[[i]]))

  return(out)
}

# Apply analysis to all canonical motifs
G4s_canonical_GC_list <- lapply(G4s_canonical$sequence, analyze_G4_GC_canonical)
G4s_canonical_GC_df <- bind_rows(G4s_canonical_GC_list)
G4s_canonical_GC <- bind_cols(G4s_canonical, G4s_canonical_GC_df)

# Map GC content across individual loop positions
loop_positions <- G4s_canonical_GC %>%
  select(chr, start, end, loop1_GC_pos, loop2_GC_pos, loop3_GC_pos) %>%
  pivot_longer(cols = starts_with("loop"),
               names_to = "loop_id",
               values_to = "GC_pos_list") %>%
  unnest(GC_pos_list) %>%
  group_by(chr, start, end, loop_id) %>%
  mutate(position_in_loop = row_number(),
         GC = GC_pos_list) %>%
  ungroup() %>%
  mutate(loop_id = recode(loop_id,
                          "loop1_GC_pos" = "loop1",
                          "loop2_GC_pos" = "loop2",
                          "loop3_GC_pos" = "loop3"))


# Summarize GC content per loop position
gc_content <- loop_positions %>%
  group_by(loop_id, position_in_loop) %>%
  summarise(
    GC_content = mean(GC, na.rm = TRUE),
    n_sequences = n(),
    .groups = "drop"
  )


# Plot GC content for each loop
for(loop_name in c("loop1", "loop2", "loop3")){
  loop_data <- gc_content %>% filter(loop_id == loop_name)
  p <- ggplot(loop_data, aes(x = factor(position_in_loop), y = GC_content)) +
    geom_col(
      fill = case_when(
        loop_name == "loop1" ~ "darkblue",
        loop_name == "loop2" ~ "steelblue",
        TRUE ~ "skyblue"
      ),
      color = "black"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(
      title = paste0("GC Content per Position in ", toupper(loop_name)),
      x = paste0("Position in ", loop_name),
      y = "GC content"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Save as PDF
  ggsave(
    filename = paste0("G4_", loop_name, "_GC_ratio.pdf"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    device = cairo_pdf
  )
}


# ==============================
# 7. Identify and visualize top canonical G4 motifs (unique motifs)
# ==============================

# Extract canonical motif cores 
G4s_canonical <- G4s_canonical %>%
  mutate(motif_core = stringr::str_extract(sequence,
                                           "((G{3,}[ATCG]{1,12}){3,}G{3,})"))

n_distinct(G4s_canonical$motif_core)


# Count unique motifs
top_motifs <- G4s_canonical%>%
  filter(!is.na(motif_core)) %>%
  count(motif_core, sort = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(freq_percent = 100 * n / sum(n))


# Plot (horizontal, with right padding)
p_top_motifs <- ggplot(top_motifs,
                       aes(x = reorder(motif_core, n),
                           y = n,
                           fill = n)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", round(freq_percent, 1), "%)")),
            hjust = -0.1, size = 4, family = "sans") +
  scale_fill_gradient(low = "#BFD7ED", high = "#2F6690") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Top 20 Canonical G4 Motifs",
    subtitle = "Frequency among detected canonical G4 sequences",
    x = NULL,
    y = "Counts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray30"),
    axis.text.y = element_text(family = "monospace", size = 12),
    axis.text.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 80, 10, 10)  # ← increased right margin
  )


ggsave("Top20_Canonical_G4_Motifs_G4-SNVs.pdf",
       plot = p_top_motifs,
       width = 8, height = 8,
       device = cairo_pdf)


# ==============================
# 8. Downsampling
# ==============================
# since the number of non-SNV-G4s exceeds by far the number of SNV-G4s, we randomly downsample the dataset 
# using 1000 iterations 

# Parameters
n_snv <- 1096      # number of canonical motifs for SNV set
B <- 1000          # number of random resamples

# Draw 1096 motifs 100 times (resampling)
gc_sampled <- map_dfr(seq_len(B), function(b) {
  sampled <- G4s_canonical_GC %>%
    sample_n(size = n_snv, replace = FALSE)
  sampled$iteration <- b
  sampled
})

# Convert to per-base positional table
loop_positions <- gc_sampled %>%
  select(iteration, chr, start, end,
         loop1_GC_pos, loop2_GC_pos, loop3_GC_pos) %>%
  pivot_longer(
    cols = starts_with("loop"),
    names_to = "loop_id",
    values_to = "GC_pos_list"
  ) %>%
  unnest(GC_pos_list) %>%
  group_by(iteration, chr, start, end, loop_id) %>%
  mutate(position_in_loop = row_number(),
         GC = GC_pos_list) %>%
  ungroup() %>%
  mutate(loop_id = recode(loop_id,
                          "loop1_GC_pos" = "loop1",
                          "loop2_GC_pos" = "loop2",
                          "loop3_GC_pos" = "loop3"))


# Calculate mean GC per loop position within each iteration
gc_iteration_means <- loop_positions %>%
  group_by(loop_id, position_in_loop, iteration) %>%
  summarise(GC_mean = mean(GC, na.rm = TRUE), .groups = "drop")


# Compute global mean and 95% CI across iterations
gc_summary <- gc_iteration_means %>%
  group_by(loop_id, position_in_loop) %>%
  summarise(
    mean_GC = mean(GC_mean, na.rm = TRUE),
    lower_CI = quantile(GC_mean, 0.025, na.rm = TRUE),
    upper_CI = quantile(GC_mean, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


# dotplot with 95% confidence intervals
p <- ggplot(gc_summary, aes(x = position_in_loop, y = mean_GC, color = loop_id)) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(ymin = lower_CI, ymax = upper_CI, color = loop_id),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  theme_bw() +
  labs(
    x = "Position in loop",
    y = "Mean GC content",
    title = "GC fraction per loop position (1000× random resampling)",
    color = "Loop"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

ggsave("G4s_random_sampled_1000x.pdf",
       plot = p,
       width = 9, height = 5,
       device = cairo_pdf)

write.table(
  gc_summary,
  file = "gc_summary_downsampled.tsv", #add _downsampled for second dataset or _G4_SNVs 
  sep = "\t",
  row.names = FALSE,
  col.names= TRUE,
  quote = FALSE
)


# ==============================
# 9. Combined dotplot with G4-SNV GC content
# ==============================

# gc_content calculated for G4-SNVs

# first test statistics
gc_test <- gc_iteration_means %>%
  left_join(gc_content, by = c("loop_id", "position_in_loop"))

gc_pvals <- gc_test %>%
  group_by(loop_id, position_in_loop) %>%
  summarise(
    obs_GC = unique(GC_content),
    mean_bg = mean(GC_mean),
    p_emp = mean(abs(GC_mean - mean_bg) >= abs(obs_GC - mean_bg)),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_emp, method = "BH"))  # FDR correction

# Visulatization
gc_content_sig <- gc_content %>%
  left_join(gc_pvals, by = c("loop_id", "position_in_loop"))

#with grey background:
sig_regions <- gc_pvals %>%
  filter(p_adj < 0.05) %>%
  mutate(ymin = -Inf, ymax = Inf)  # full-height grey band

p4 <- ggplot() +
  # grey highlight for significant loop positions
  geom_rect(
    data = sig_regions,
    aes(
      xmin = position_in_loop - 0.5,
      xmax = position_in_loop + 0.5,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey80",
    alpha = 0.5
  ) +
  
  # bootstrap background (CI ribbon + mean line)
  geom_ribbon(
    data = gc_summary,
    aes(x = position_in_loop, ymin = lower_CI, ymax = upper_CI, fill = loop_id),
    alpha = 0.2
  ) +
  geom_line(
    data = gc_summary,
    aes(x = position_in_loop, y = mean_GC, color = loop_id),
    size = 1.2
  ) +
  
  # observed SNV-G4 points
  geom_point(
    data = gc_content,
    aes(x = position_in_loop, y = GC_content),
    color = "black",
    size = 2.5
  ) +
  ylim(0.2,0.8) +
  facet_wrap(~ loop_id, scales = "free_x") +
  theme_bw(base_size = 16) +  # increase overall font size
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +  # x-axis labels at 2,4,6,8,10,12
  labs(
    x = "Position in loop",
    y = "GC fraction",
    title = "SNV-G4 GC content vs. random background"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )


ggsave("G4s_random_sampled_1000x_combined_statistics_grey.pdf",
       plot = p4,
       width = 9, height = 5,
       device = cairo_pdf)

write.table(
  gc_pvals,
  file = "gc_statistics.tsv", 
  sep = "\t",
  row.names = FALSE,
  col.names= TRUE,
  quote = FALSE
)
