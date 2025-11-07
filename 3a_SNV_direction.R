# ==============================================================================
# G4-SNV Mutation Spectrum and Strand Bias Analysis
# ==============================================================================
# Author: Dr. Angelika Lahnsteiner
# Date: 2025-10-26
# Purpose:
#   Compare substitution direction of G4-associated SNPs on plus and minus strands,
#   assess differences between risk and protective associations,
#   analyze control SNVs, and evaluate G4-SNV orientation relative to the strand.
# ==============================================================================

# ------------------------------
# Load libraries
# ------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# ------------------------------
# Load data
# ------------------------------
# Load G4-SNVs (choose PQSfinder or G4Hunter)
df <- read.delim("path_to_your_file/G4_SNPs.bed", header = TRUE) #pqsfinder or G4Hunter predicted G4 motifs overlapping with SNVs

df.plus  <- df[df$G4_strand == "+",] # select G4s located on the plus strand
df.minus <- df[df$G4_strand == "-",] # selest G4s located on the minus

# ==============================================================================
# 1. Adjust base changes for minus strand (complementary substitution)
# ==============================================================================
complement_base <- function(base) {
  base <- toupper(base)
  comps <- c(A = "T", T = "A", C = "G", G = "C")
  comps[base]
}

df.minus$base_change <- sapply(df.minus$base_change, function(x) {
  parts <- unlist(strsplit(x, ">"))
  if (length(parts) == 2) {
    paste0(complement_base(parts[1]), ">", complement_base(parts[2]))
  } else {
    NA
  }
})

# ==============================================================================
# 2. Compute base-change frequencies for G4 (+) strand
# ==============================================================================
df.risk <- df.plus %>% filter(beta_value > 0) # dataframe contains all risk SNVs
change.summary.risk <- as.data.frame(table(df.risk$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "risk", Occurence = "G4")

df.prot <- df.plus %>% filter(beta_value < 0) # dataframe contains all protective SNVs
change.summary.prot <- as.data.frame(table(df.prot$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "protective", Occurence = "G4")

G4.change.plus <- rbind(change.summary.risk, change.summary.prot) %>%
  mutate(Strand = "+")

# ==============================================================================
# 3. Compute base-change frequencies for G4 (−) strand
# ==============================================================================
df.risk <- df.minus %>% filter(beta_value > 0) # dataframe contains all risk SNVs
change.summary.risk <- as.data.frame(table(df.risk$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "risk", Occurence = "G4")

df.prot <- df.minus %>% filter(beta_value < 0) # dataframe contains all protective SNVs
change.summary.prot <- as.data.frame(table(df.prot$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "protective", Occurence = "G4")

G4.change.minus <- rbind(change.summary.risk, change.summary.prot) %>%
  mutate(Strand = "-")

# ==============================================================================
# 4. Control dataset: Non-G4 SNVs matched for GC content
# ==============================================================================
all.SNPs <- read.delim("path_to_your_file/sign_SNPs_hg38.bed") # Park MetD-associated SNV dataset

control <- all.SNPs %>%
  filter(!(rs_id %in% df$rs_id)) %>%
  filter(GC_content_300bp > 0.43379)

control.risk <- control %>% filter(beta > 0)
control.prot <- control %>% filter(beta < 0)

change.summary.risk <- as.data.frame(table(control.risk$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "risk", Occurence = "control")

change.summary.prot <- as.data.frame(table(control.prot$base_change)) %>%
  setNames(c("Change", "Count")) %>%
  mutate(Ratio = round(Count / sum(Count), 4),
         Association = "protective", Occurence = "control")

control.change <- rbind(change.summary.risk, change.summary.prot) %>%
  mutate(Strand = "*")

# ==============================================================================
# 5. Merge datasets and prepare for plotting
# ==============================================================================
change.summary <- rbind(G4.change.plus, G4.change.minus, control.change) %>%
  mutate(
    FillGroup = case_when(
      Occurence == "control" ~ "control",
      Association == "risk" ~ "G4_pos",
      TRUE ~ "G4_neg"
    ),
    Occurence   = factor(Occurence, levels = c("G4", "control")),
    Association = factor(Association, levels = c("risk", "protective")),
    Proportion  = Ratio * 100
  )

# ==============================================================================
# 6. Plot mutation spectrum by risk/protective association
# ==============================================================================
p1 <- ggplot(change.summary, aes(x = Change, y = Proportion, fill = FillGroup)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~Association, scales = "free_y") +
  scale_fill_manual(values = c("G4_pos" = "darkred", "G4_neg" = "darkblue", "control" = "gray")) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y  = element_text(face = "bold"),
    strip.text   = element_text(face = "bold"),
    axis.title   = element_text(face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "grey", fill = NA, size = 1)
  ) +
  labs(y = "Proportion [%]", title = "Mutation spectrum in G4-SNVs and controls")

ggsave("G4_mutation_spectrum_by_association.png", p1, width = 8, height = 5, dpi = 300)

# ==============================================================================
# 7. Plot mutation spectrum including strand differentiation
# ==============================================================================
df.strand <- change.summary %>%
  mutate(Group = case_when(
    Occurence == "G4" & Strand == "+" ~ "G4+",
    Occurence == "G4" & Strand == "-" ~ "G4-",
    Occurence == "control" ~ "control"
  ),
  Proportion = Ratio * 100) %>%
  mutate(Group = factor(Group, levels = c("G4+", "G4-", "control")))

write.table(df.strand, "~/Desktop/2025-10-23-change_summary_G4Hunter.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

p2 <- ggplot(df.strand, aes(x = Change, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~Association) +
  scale_fill_manual(values = c("G4+" = "darkred", "G4-" = "darkblue", "control" = "gray")) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    legend.position = "bottom",
    panel.border = element_rect(color = "grey", fill = NA, size = 1)
  ) +
  labs(y = "Proportion [%]", title = "Mutation spectrum by strand and association")

ggsave("G4_mutation_spectrum_by_strand.png", p2, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 8. Compute strand difference (Δratio = G4+ − G4−)
# ==============================================================================
df.diff <- rbind(G4.change.plus, G4.change.minus) %>%
  group_by(Association, Change, Strand) %>%
  summarise(Ratio = mean(Ratio), .groups = "drop") %>%
  pivot_wider(names_from = Strand, values_from = Ratio, values_fill = 0) %>%
  mutate(diff_ratio = `+` - `-`, Proportion = diff_ratio * 100)

p3 <- ggplot(df.diff, aes(x = Change, y = Proportion)) +
  geom_bar(stat = "identity", fill = "darkred") +
  facet_wrap(~Association) +
  theme_minimal(base_size = 13) +
  labs(y = "Δ Proportion [%]", title = "Difference in +/− strand mutation frequencies") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1)
  )

ggsave("G4_mutation_delta_strand.png", p3, width = 8, height = 5, dpi = 300)

# ==============================================================================
# 9. Template vs. non-template strand analysis
# ==============================================================================
gr.snps <- GRanges(seqnames = df$chr,
                   ranges = IRanges(start = df$start_G4, end = df$end_G4),
                   rs_id = df$rs_id)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annotated <- annotatePeak(gr.snps, TxDb = txdb, tssRegion = c(-1000, 100), level = "gene")
df_annotated <- as.data.frame(annotated) %>%
  mutate(GeneStrand = ifelse(geneStrand == 1, "+", "-"),
         Annotation = gsub("\\s*\\(.*\\)", "", annotation))

G4.SNV <- cbind(df, df_annotated[, c("annotation", "geneId", "distanceToTSS", "geneStrand")]) %>%
  mutate(Association = ifelse(beta_value > 0, "risk", "protective"),
         TemplateStatus = ifelse(G4_strand == GeneStrand, "non-template", "template"))

# ==============================================================================
# 10. Plot template vs. non-template proportions
# ==============================================================================
df.template <- G4.SNV %>%
  group_by(Association, TemplateStatus) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Association) %>%
  mutate(Proportion = 100 * Count / sum(Count))

p.template <- ggplot(df.template, aes(x = TemplateStatus, y = Proportion, fill = TemplateStatus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Association) +
  scale_fill_manual(values = c("template" = "black", "non-template" = "gray")) +
  theme_minimal(base_size = 13) +
  labs(y = "Proportion [%]", title = "G4-SNV orientation relative to gene strand") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1)
  )

ggsave("G4_template_vs_nontemplate.png", p.template, width = 7, height = 5, dpi = 300)

# ==============================================================================
# 11. SNV distance to TSS (density distribution)
# ==============================================================================
p.tss <- ggplot(G4.SNV, aes(x = distanceToTSS / 1000, fill = TemplateStatus)) +
  geom_density(data = subset(G4.SNV, TemplateStatus == "template"),
               aes(y = ..density..), alpha = 0.6, fill = "firebrick") +
  geom_density(data = subset(G4.SNV, TemplateStatus == "non-template"),
               aes(y = -..density..), alpha = 0.6, fill = "steelblue") +
  labs(x = "Distance to TSS (kb)", y = "Density",
       title = "G4-SNV distance to TSS by strand orientation") +
  xlim(-1000, 1000) + ylim(-0.015, 0.015) +
  theme_minimal(base_size = 14)

ggsave("G4_distance_to_TSS_density.png", p.tss, width = 8, height = 5, dpi = 300)

cat("\n G4-SNV mutation and strand bias analysis completed.\n")
cat("Plots saved:\n",
    " - G4_mutation_spectrum_by_association.png\n",
    " - G4_mutation_spectrum_by_strand.png\n",
    " - G4_mutation_delta_strand.png\n",
    " - G4_template_vs_nontemplate.png\n",
    " - G4_distance_to_TSS_density.png\n", sep = "")
