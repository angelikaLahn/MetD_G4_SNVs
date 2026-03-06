########### Calculation of pqsfinder and G4Hunter Scores ###########
#Install pqsfinder and G4SNVHunter package with Bioconductor before use

#written by Victoria Ellmer

library(ggplot2)
library(pqsfinder)
library(G4SNVHunter)

#Calculate the scores for both sequences containing the effect and the other allele

#load data first: 
MutatedPQS <- as.data.frame(read.table("MutatedG4s_pqsfinder.bed"))
names(MutatedPQS) <- MutatedPQS[1,]
MutatedPQS <- MutatedPQS[-1,]

MutatedHunter <- as.data.frame(read.table("MutatedG4s_G4Hunter.bed"))
names(MutatedHunter) <- MutatedHunter[1,]
MutatedHunter <- MutatedHunter[-1,]

#### Calculate pqsFinder Scores ####

pqsScoresMetD <- data.frame(other = rep(0,length(MutatedPQS$chr)),effect = rep(0,length(MutatedPQS$chr)),rsid=rep(0,length(MutatedPQS$chr)))


for (i in 1:length(MutatedPQS$chr)) {
  print(i)
  pqsScoresMetD$rsid[i] <- MutatedPQS$rsID[i]
  a <- max(score(pqsfinder(DNAString(MutatedPQS$otherSequence[i]),min_score =1))) #min_score = 1
  #print(a)
  if (is.na(a) ){ # if pqsfinder classifies the sequence to be no G4, score of 0 is given 
    pqsScoresMetD$other[i] <- 0
  }
  else {
    pqsScoresMetD$other[i] <- a
  }
  b <- max(score(pqsfinder(DNAString(MutatedPQS$effectSequence[i]), min_score =1)))
  if (is.na(b)) {
    pqsScoresMetD$effect[i] <- 0
  }
  else {
    pqsScoresMetD$effect[i] <- b
  }
}

#add additional info
pqsScoresMetD$diff <- pqsScoresMetD$effect - pqsScoresMetD$other
pqsScoresMetD$chr <- MutatedPQS$chr
pqsScoresMetD$start_G4 <- MutatedPQS$start_G4
pqsScoresMetD$end_G4 <- MutatedPQS$end_G4
pqsScoresMetD$beta_value <- MutatedPQS$beta_value
pqsScoresMetD$standard_error <- MutatedPQS$standard_error
pqsScoresMetD$frequency <- MutatedPQS$frequency
pqsScoresMetD$p_value <- MutatedPQS$`p-value`

#rearrange columns:
pqsScoresMetD <- pqsScoresMetD[c(3,5,6,7,1,2,4,8,9,10,11)]
#save data
write.table(pqsScoresMetD,"pqsScoresMetD_hg38_pval.txt",quote=F,row.names = F, sep = "\t")

#### calculate G4HunterScores ####

#Calculate the G4Hunter scores for overlaps detected with the G4Hunter
HunterScoresMetD <- data.frame(other = rep(0,length(MutatedHunter$chr)),effect = rep(0,length(MutatedHunter$chr)),rsid=rep(0,length(MutatedHunter$chr)))

for (i in 1:length(MutatedHunter$chr)){
  HunterScoresMetD$rsid[i] <- MutatedHunter$rsID[i]
  HunterScoresMetD$other[i] <- G4HunterScore(MutatedHunter$otherSequence[i])
  HunterScoresMetD$effect[i] <- G4HunterScore(MutatedHunter$effectSequence[i])
  
}
HunterScoresMetD$diff <- abs(HunterScoresMetD$effect) - abs(HunterScoresMetD$other)


HunterScoresMetD$chr <- MutatedHunter$chr
HunterScoresMetD$start_G4 <- MutatedHunter$start_G4
HunterScoresMetD$end_G4 <- MutatedHunter$end_G4
HunterScoresMetD$beta_value <- MutatedHunter$beta_value
HunterScoresMetD$standard_error <- MutatedHunter$standard_error
HunterScoresMetD$frequency <- MutatedHunter$frequency
HunterScoresMetD$p_value <- MutatedHunter$`p-value`


HunterScoresMetD <- HunterScoresMetD[c(3,5,6,7,1,2,4,8,9,10,11)]
#save data
write.table(HunterScoresMetD,"HunterScoresMetD_hg38_pval.txt",quote=F,row.names = F, sep = "\t")





