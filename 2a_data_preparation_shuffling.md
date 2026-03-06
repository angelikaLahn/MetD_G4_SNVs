---
title: "File preparation for random shuffling"
author: "angelikalahnsteiner"
---


#To prepare the files for random sampling, the following steps are conducted:

###############################################################################################
#PREPARE THE GENOME FILE hg38 (UCSC)
###############################################################################################


#1. segment the genome in 50bp windows
```bash
bedtools makewindows \
-g hg38_ucsc.fa.fai \
-w 50 \
> genome_WS50bp.bed
```

#2. remove sex chromosomes since they are not included in the MetD-associated SNV file
```bash
grep -v -E '^chrX|^chrY' genome_WS50bp.bed \
> genome_WS50bp_autosomes.bed
```

#3. download gap file from https://42basepairs.com/browse/s3/giab/technical/ucsc_bed_files?file=hg38.gap.bed.gz&preview=
#save the downloaded gap file as bed file 

```bash
awk 'BEGIN{OFS="\t"} {print $2, $3, $4}' gap.txt > gaps.bed
```

#4. remove gaps from the genome file
```bash
bedtools intersect \
-v \
-a genome_WS50bp_autosomes.bed \
-b gaps.bed \
> genome_WS50bp_autosomes_noGaps.bed
```

#check sequences for each entry --> no NNNNs in the bins included
```bash
bedtools getfasta \
-fi hg38_ucsc.fa \
-bed genome_WS50bp_autosomes_noGaps.bed \
-tab \
> genome_WS50bp_autosomes_sequences.txt
```

#5. Calculate nucleotide composition
```bash
bedtools nuc \
-fi hg38_ucsc.fa \
-bed genome_WS50bp_autosomes_noGaps.bed \
> genome_WS50bp_GC.txt
```

#6. save as bed file with positions and GC content
```bash
cut -f1,2,3,5 genome_WS50bp_GC.txt > genome_WS50bp_GC.bed
```

#remove header line
```bash
awk 'NR>1' genome_WS50bp_GC.bed > genome_WS50bp_GC_noH.bed
```

#7. bin the genome into 1-10 bins, with each 0.1 step increase
#bin 1: 0.0-0.1/bin 2: 0.1-0.2/../bin 10:0.9-1.0

```bash
awk 'BEGIN{OFS="\t"}{
    gc = $4
    bin = int(gc*10) + 1
    if (bin > 10) bin = 10
    print $1, $2, $3, $4, bin
}' genome_WS50bp_GC_noH.bed > genome_50bp_GC_binned.bed
```
#This file is now used as reference genome with the GC bin information for random shuffling.



###############################################################################################
# PREPARE THE SNV FILE
###############################################################################################
# this is done with the whole MetD-associated SNV datasets, as well as with the individual
# pqsfinder and G4Hunter G4-SNVs.
# just replace the file names.

#1. Create 50 bp windows centered around SNVs 
```bash
bedtools slop \
-i sign_SNPs_hg38.bed \
-g hg38_ucsc.genome \
-b 25 \
> snvs_WS50bp.bed
```

#2. save as bed file with chr and positions 
```bash
cut -f1-3,10 snvs_WS50bp.bed > snvs_WS50bp_short.bed
```


#4. Run bedtools nuc to get the nucleotide composition for each SNV and its +/-25bp windows
```bash
bedtools nuc \
-fi hg38_ucsc.fa \
-bed snvs_WS50bp_short.bed \
> snvs_WS50bp_nuc_comp.txt
```


#5. save as bed file with positions and GC content
```bash
cut -f1,2,3,4,6 snvs_WS50bp_nuc_comp.txt > snvs_WS50bp_nuc_comp.bed
```

#6. remove header line
```bash
awk 'NR>1' snvs_WS50bp_nuc_comp.bed > snvs_WS50bp_nuc_comp_noH.bed
```

#7. bin the snvs into 1-10 bins, with each 0.1 step increase
# bin 1: 0.0-0.1/bin 2: 0.1-0.2/../bin 10:0.9-1.0
```bash
awk 'BEGIN{OFS="\t"}{
    gc = $4
    bin = int(gc*10) + 1
    if (bin > 10) bin = 10
    print $1, $2, $3, $4, bin
}' snvs_WS50bp_nuc_comp_noH.bed > snvs_50bp_GC_binned.bed
```

## 8. save as bed file with positions and GC content
```bash
cut -f4,5 snvs_50bp_GC_binned.bed > sign_SNPs_hg38_binned.bed
```

