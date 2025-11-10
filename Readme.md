# G4–SNV Overlap Analysis

Custom scripts for detecting and analyzing overlaps between **MetD-associated single-nucleotide variants (SNVs)** and **predicted G-quadruplex (G4) motifs** using **pqsfinder** [1] and **G4Hunter** [2] predictions.

## Overview

This repository contains scripts used to identify, classify, and analyze overlaps between SNVs associated with metabolic disorders (MetDs) identified by Park et al. [3] and predicted G4 motifs in the human genome (**GRCh38**).  
The pipeline supports both **pqsfinder**- and **G4Hunter**-predicted G4s and enables downstream analyses including canonical G4 motif identification, SNV classification by motif region (stem or loop), GC content evaluation, and overlap with transcriptional features.


## 1. Prediction of G4 Motifs

### pqsfinder
- G4 annotations and stability scores were obtained from the [pqsfinder genomes database](https://pqsfinder.fi.muni.cz/genomes) [1]
- Minimum stability score threshold: **47** (default parameters used).  

### G4Hunter
- G4 motifs predicted using the **`G4HunterDetect`** function from the `G4SNVHunter` R package [4].  
- Parameters:
  - Window size: **25**
  - Absolute G4Hunter score threshold: **1.5**
- Reference genome: **GRCh38** (downloaded from UCSC).


## 2. Detection of G4–SNV Overlaps

Overlaps between **MetD-associated SNVs** and **G4 motifs** were computed using **`bedtools intersect`**.

To assess whether overlaps occur by chance:
1. Chromosome-wise random SNV sets were generated from the full SNV dataset.
2. The randomization was repeated **1,000 times** to establish a null distribution.
3. P-values were computed using a **permutation test**.


## 3. Canonical G4 Motif Identification

Canonical G4 motifs were defined using the following regular expression pattern:
(G{3,})([ATCG]{1,12})(G{3,})([ATCG]{1,12})(G{3,})([ATCG]{1,12})(G{3,})

Only motifs with **≥3 consecutive guanines per stem** and **1–12 bp loops** were considered canonical.  
All analyses were performed in **R v4.5.0 (2025-04-11)**.

## 4. SNV Classification
Each SNV was classified based on the presence of a canonical G4 motif in the **effect** and **non-effect** alleles:

- Present in both alleles → **No change**
- Present only in the non-effect allele → **Motif loss**
- Present only in the effect allele → **Motif gain**

SNVs were further annotated by location:
- **Stem**
- **Loop**


## 5. Additional Analyses

### GC Content Analysis
Mean GC fraction and 95% confidence intervals of G4 loop regions were computed using **1,000 iterations**.

### Strand Bias and Substitution Direction
Mutation spectra and strand asymmetry were computed separately for G4s on the **plus** and **minus** strands, distinguishing between **risk-associated** and **protective** variants.

### Transcription Factor Overlaps
G4 sequences were intersected with **JASPAR (2024)** [5] transcription factor binding motifs using R.

### Gene Ontology (GO) Enrichment
Molecular Function (MF) enrichment was performed using the `clusterProfiler` package  
(`OrgDb = org.Hs.eg.db`, `ontology = "MF"`, FDR < 0.05).

### Alternative Promoter Overlaps
G4–SNVs were intersected with ±250 bp windows around alternative promoters  
(from *Demircioglu et al.*, 2020) [6].  

## 6. Dependencies

- **R ≥ 4.5.0**
  - `pqsfinder`
  - `G4SNVHunter`
  - `ChIPseeker`
  - `GenomicRanges`
  - `annotatr`
  - `clusterProfiler`
  - `org.Hs.eg.db`
- **bedtools ≥ 2.30.0**
- Reference genome: **GRCh38**



## Sources:
[1] Hon, J., et al., pqsfinder: an exhaustive and imperfection-tolerant search tool for potential quadruplex-forming sequences in R. Bioinformatics, 2017. 33(21): p. 3373-3379.

[2] Bedrat, A., L. Lacroix, and J.L. Mergny, Re-evaluation of G-quadruplex propensity with G4Hunter. Nucleic Acids Res, 2016. 44(4): p. 1746-59.

[3] Park, S., et al., Multivariate genomic analysis of 5 million people elucidates the genetic architecture of shared components of the metabolic syndrome. Nat Genet, 2024. 56(11): p. 2380-2391.

[4] Zhang, R., G4SNVHunter: Evaluating SNV-Induced Disruption of G-Quadruplex Structures, R package version 1.1.4

[5] Rauluseviciute, I., et al., JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res, 2024. 52(D1): p. D174-D182.

[6] Demircioglu, D., et al., A Pan-cancer Transcriptome Analysis Reveals Pervasive Regulation through Alternative Promoters. Cell, 2019. 178(6): p. 1465-1477 e17
