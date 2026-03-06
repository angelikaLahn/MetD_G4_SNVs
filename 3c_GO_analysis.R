###############################################################################
# GO Enrichment of Genes Near G4-SNPs 
# Author: Angelika Lahnsteiner, Dr.
# Date: 2025-10-26
# Description:
#   - Load G4-SNP tables (PQSfinder or G4Hunter)
#   - Genomic annotation (promoters/genes) with annotatr
#   - Map SNVs to genes
#   - Perform GO enrichment (clusterProfiler)
#   - prepare a dotplot
###############################################################################

# ==============================
# 1. Load libraries 
# ==============================
suppressPackageStartupMessages({
  library(DOSE)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(IRanges)
  library(annotatr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(ggnewscale)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(readr)
})

# ==============================
# 2. Load data
# ==============================
df <- read.delim("path_to_your_file/G4_SNPs.bed", header = TRUE)

# ==============================
# 3. Genomic annotation (annotatr)
# ==============================
gr.G4.SNPs <- GRanges(seqnames = df$chr,
                      ranges   = IRanges(start = df$start, end = df$end),
                      rs_id    = df$rs_id)

annotations <- c("hg38_basicgenes",
                 "hg38_genes_intergenic",
                 "hg38_genes_intronexonboundaries")

annots <- build_annotations(genome = "hg38", annotations = annotations)
gr.annot <- annotate_regions(regions = gr.G4.SNPs, annotations = annots)
G4.SNPs.annot <- as.data.frame(gr.annot)

# Keep one annotation per rs_id if duplicates; trim annot.id to gene name/id core
G4.SNPs.annot$annot.id <- sub(":.*", "", G4.SNPs.annot$annot.id)
G4.SNPs.annot <- G4.SNPs.annot %>%
  distinct(rs_id, annot.symbol, .keep_all = TRUE)

# Extract genes 
genes <- G4.SNPs.annot %>% pull(annot.symbol) %>% unique() %>% na.omit()

# ==============================
# 4. GO enrichment (Molecular Function)
# ==============================
# Helper to run enrichGO safely (returns empty tibble if nothing enriched)
run_enrichGO_safe <- function(genes_vec) {
  if (length(genes.vec <- unique(genes_vec)) == 0) return(NULL)
  res <- tryCatch(
    enrichGO(gene          = genes.vec,
             keyType       = "SYMBOL",
             OrgDb         = org.Hs.eg.db,
             ont           = "MF",
             pAdjustMethod = "BH",
             qvalueCutoff  = 0.05,
             readable      = TRUE),
    error = function(e) NULL
  )
  res
}

ego <- run_enrichGO_safe(genes)

# ==============================
# 5. Tidy enrichment tables and compute RichFactor/logP
# ==============================

# RichFactor = Count / setSize, otherwise Count / size from BgRatio
compute_richfactor <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(tibble())
  df2 <- res_df@result
  if (!"RichFactor" %in% names(df2)) {
    if ("setSize" %in% names(df2)) {
      df2$RichFactor <- df2$Count / df2$setSize
    } else if ("BgRatio" %in% names(df2)) {
      denom <- as.numeric(sub(".*/", "", df2$BgRatio))
      df2$RichFactor <- df2$Count / denom
    } else {
      df2$RichFactor <- df2$Count
    }
  }
  df2$logP <- -log10(df2$p.adjust + 1e-300)
  as_tibble(df2)
}

tbl <- compute_richfactor(ego) 

# Take top 15 per group by Count 
top <- tbl %>% arrange(desc(Count)) %>% slice_head(n = 15)


# ==============================
# 5. Dotplot of top enriched MF terms
# ==============================
p <- ggplot(top, aes(x = RichFactor,
                               y = reorder(Description, RichFactor),
                               size = Count, color = logP)) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
    colors = c("#dceeff", "#91bfdb", "#4575b4", "#313695"),
    name   = expression(-log[10]("BH adj.p"))
  ) +
  theme_bw(base_size = 13) +
  labs(x = "Rich Factor", y = "GO term") +
  theme(axis.text.y = element_text(size = 10, face = "italic"),
        plot.title  = element_text(face = "bold", hjust = 0.5))


ggsave("GO_MF_dotplots.png", p, width = 8, height = 9, dpi = 300)

