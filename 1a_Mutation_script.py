"""
                Mutation script
This script was used to mutate the identified G4 sequences which harbored a MetD risk SNP.
To mutate the G4 sequences we need:
        - a fasta file with all the SNPs within the annotated G4 sequences (bedtools intersect SNP.bed G4.bed; use bedtools getfasta after to retrieve a fasta file containing all the sequences)
        - a bed file with all the G4s harboring a SNP (bedtools intersect G4.bed SNP.bed)
Note that both files need to be sorted in a "natural order" (starting with chr1 lowest position; if the files are not sorted, use the bedtools sort function with additional genomic file (fa.fai) describing the order the files should be sorted in)

"""

# written by Victoria Ellmer, MSc.

#write a function to obtain the G4 sequence with the effect and with the non - effect allel
def mutateSequenceNew(chromG4,startG4,endG4,sequence,chromSNP,posSNP,otherSNP,effSNP):
    if chromG4 == chromSNP and startG4 <= posSNP and endG4 >= posSNP: #check if SNP is within the G4 sequence
        diff = posSNP - startG4 #lets get the position of our SNP respective to the G4
        sequence = sequence.upper() # we need to be sure that all characters are upper case! 
        if otherSNP == sequence[diff] or effSNP == sequence[diff]: #check if the base at the SNP position is either the effect or other allele
            otherSequence = sequence[:diff] + otherSNP + sequence[diff+1:]	
            effSequence = sequence[:diff] + effSNP + sequence[diff + 1:]
            return(otherSequence,effSequence)
        else:
            print("The given base in the G4 sequence is not the annotated as effect or non - effect allele.")
    else:
        print("The positions do not match. Are you sure you ordered the files beforehand?")
        

#write a function to obtain the G4 sequence with the effect and with the non - effect allel
#this function does not break when the base at the SNP position is not either the effect or non - effect allele 
def mutateSequenceModified(chromG4,startG4,endG4,sequence,chromSNP,posSNP,refSNP,snpSNP):
    if chromG4 == chromSNP and startG4 <= posSNP and endG4 >= posSNP:
        diff = posSNP - startG4 #lets get the position of our SNP respective to the G4
        sequence = sequence.upper() # we need to be sure that all characters are upper case! 
        if refSNP == sequence[diff] or snpSNP == sequence[diff]: 
            refSequence = sequence[:diff] + refSNP + sequence[diff+1:]	
            mutatedSequence = sequence[:diff] + snpSNP + sequence[diff + 1:]
            return(refSequence,mutatedSequence)
        else:
            print("The given base in the G4 sequence is not the annotated as effect or non - effect allele.")
            refSequence = sequence[:diff] + refSNP + sequence[diff+1:]	
            mutatedSequence = sequence[:diff] + snpSNP + sequence[diff + 1:]
            return(refSequence,mutatedSequence)
    else:
        print("The positions do not match. Are you sure you ordered the files beforehand?")


"""
Lets start with G4s predicted by pqsFinder:
    
"""
chrom_G4 = []
start_G4 = []
end_G4 = []
seq_G4 = []
with open("G4_overlap_pqsFinder.fa", "r+") as fasta:
    for line in fasta.readlines():
        if ">" in line:
            chro = line.split(":")[0][1:] #get chromosome
            info = line.split(":")[1]
            start = info.split("-")[0]
            end = info.split("-")[1]
            chrom_G4.append(chro)
            start_G4.append(int(start))
            end_G4.append(int(end))
        else:
            seq_G4.append(line.strip())

G4InfoMetD = {
    "chrom":chrom_G4,
    "start":start_G4,
    "end":end_G4,
    "seq":seq_G4
    }

#load data of SNPs overlapped with pqsFinder

chromSNP = []
posSNP = []
infoSNP = []
otherSNP = []
effSNP = []  
beta_value =[]
standard_error =[]
frequency = []
p_value = []

with open("SNP_overlap_pqsFinder_sorted.bed", "r+") as bed:
    for line in bed.readlines():
            chromSNP.append(line.split()[0]) #get chromosome
            posSNP.append(int(line.split()[1])) #get start coordinate
            infoSNP.append(line.split()[9])
            otherSNP.append(line.split()[4])
            effSNP.append(line.split()[3])
            beta_value.append(line.split()[5])
            standard_error.append(line.split()[6])
            frequency.append(line.split()[7])
            p_value.append(line.split()[8])

SNPInfoMetD = {
    "chrom":chromSNP,
    "pos":posSNP,
    "info":infoSNP,
    "other":otherSNP,
    "eff":effSNP,
    "beta":beta_value,
    "SE":standard_error,
    "frequency":frequency,
    "pValue":p_value
    }

"""
#Testing:        
        
mutateSequenceNew(G4InfoMetD["chrom"][0], G4InfoMetD["start"][0], G4InfoMetD["end"][0], G4InfoMetD["seq"][0],
                                SNPInfoMetD["chrom"][0], SNPInfoMetD["pos"][0], SNPInfoMetD["other"][0], SNPInfoMetD["eff"][0])

mutateSequenceNew(G4InfoMetD["chrom"][2], G4InfoMetD["start"][2], G4InfoMetD["end"][2], G4InfoMetD["seq"][2],
                                SNPInfoMetD["chrom"][2], SNPInfoMetD["pos"][2], SNPInfoMetD["other"][2], SNPInfoMetD["eff"][2])

"""
refG4MetD = []
mutG4MetD = []

for i in range(0,len(G4InfoMetD["chrom"])):
    print(i)
    mutant = mutateSequenceNew(G4InfoMetD["chrom"][i], G4InfoMetD["start"][i], G4InfoMetD["end"][i], G4InfoMetD["seq"][i],
                                    SNPInfoMetD["chrom"][i], SNPInfoMetD["pos"][i], SNPInfoMetD["other"][i], SNPInfoMetD["eff"][i])
    refG4MetD.append(mutant[0])
    mutG4MetD.append(mutant[1])


#save file: 
outfile = open("MutatedG4s_MetD_pqsfinder.bed","a")
#we make a BED file with the following settings: chromosome, start, stop, mutated sequence
s = "chr" + "\t" + "start_G4" + "\t" + "end_G4" + "\t" + "rsID" + "\t" + "otherSequence" + "\t" + "effectSequence" + "\t" + "ohter_allele" + "\t" + "effect_allele" + "\t" + "beta_value" + "\t" + "standard_error" + "\t" + "frequency" + "\t" + "p-value"
outfile.write(s + "\n")
for i in range(0,len(G4InfoMetD["chrom"])):
    s = G4InfoMetD["chrom"][i] + "\t" + str(G4InfoMetD["start"][i]) + "\t" + str(G4InfoMetD["end"][i]) + "\t"  + SNPInfoMetD["info"][i] +"\t" + refG4MetD[i].upper() + "\t" + mutG4MetD[i].upper()  + "\t" + SNPInfoMetD["other"][i].upper() + "\t" + SNPInfoMetD["eff"][i].upper() + "\t" + SNPInfoMetD["beta"][i] + "\t" + SNPInfoMetD["SE"][i] + "\t" + SNPInfoMetD["frequency"][i] + "\t" + SNPInfoMetD["pValue"][i]                
    outfile.write(s + "\n")
outfile.close()

"""
Mutate G4s predicted with G4Hunter:
    
"""
chrom_G4 = []
start_G4 = []
end_G4 = []
seq_G4 = []

with open("G4_overlap_G4Hunter.fa", "r+") as fasta:
    for line in fasta.readlines():
        if ">" in line:
            #s = str() #leerer string
            chro = line.split(":")[0][1:] #get chromosome
            info = line.split(":")[1]
            start = info.split("-")[0]
            end = info.split("-")[1]
            chrom_G4.append(chro)
            start_G4.append(int(start))
            end_G4.append(int(end))
        else:
            s = line.strip()
            seq_G4.append(s)
    
G4InfoMetD = {
    "chrom":chrom_G4,
    "start":start_G4,
    "end":end_G4,
    "seq":seq_G4
    }

#load SNP data

chromSNP = []
posSNP = []
infoSNP = []
otherSNP = []
effSNP = []  
beta_value =[]
standard_error =[]
frequency = []
p_value = []

with open("SNP_overlap_G4Hunter_sorted.bed", "r+") as fasta:
    for line in fasta.readlines():
            chromSNP.append(line.split()[0]) #get chromosome
            posSNP.append(int(line.split()[1])) #get start coordinate
            infoSNP.append(line.split()[9])
            otherSNP.append(line.split()[4])
            effSNP.append(line.split()[3])
            beta_value.append(line.split()[5])
            standard_error.append(line.split()[6])
            frequency.append(line.split()[7])
            p_value.append(line.split()[8])

SNPInfoMetD = {
    "chrom":chromSNP,
    "pos":posSNP,
    "info":infoSNP,
    "other":otherSNP,
    "eff":effSNP,
    "beta":beta_value,
    "SE":standard_error,
    "frequency":frequency,
    "pValue":p_value
    }

#Test for rs12123298:
mutateSequenceModified(G4InfoMetD["chrom"][5257], G4InfoMetD["start"][5257], G4InfoMetD["end"][5257], G4InfoMetD["seq"][5257],
                                SNPInfoMetD["chrom"][5257], SNPInfoMetD["pos"][5257], SNPInfoMetD["other"][5257], SNPInfoMetD["eff"][5257])


refG4MetD = []
mutG4MetD = []

for i in range(0,len(G4InfoMetD["chrom"])):
    print(i)
    mutant = mutateSequenceModified(G4InfoMetD["chrom"][i], G4InfoMetD["start"][i], G4InfoMetD["end"][i], G4InfoMetD["seq"][i],
                                    SNPInfoMetD["chrom"][i], SNPInfoMetD["pos"][i], SNPInfoMetD["other"][i], SNPInfoMetD["eff"][i])
    refG4MetD.append(mutant[0])
    mutG4MetD.append(mutant[1])

#save the file:
outfile = open("MutatedG4s_G4Hunter.bed","a")
s = "chr" + "\t" + "start_G4" + "\t" + "end_G4" + "\t" + "rsID" + "\t" + "otherSequence" + "\t" + "effectSequence" + "\t" + "ohter_allele" + "\t" + "effect_allele" + "\t" + "beta_value" + "\t" + "standard_error" + "\t" + "frequency" + "\t" + "p-value"
outfile.write(s + "\n")
for i in range(0,len(G4InfoMetD["chrom"])):
    s = G4InfoMetD["chrom"][i] + "\t" + str(G4InfoMetD["start"][i]) + "\t" + str(G4InfoMetD["end"][i]) + "\t"  + SNPInfoMetD["info"][i] +"\t" + refG4MetD[i].upper() + "\t" + mutG4MetD[i].upper()  + "\t" + SNPInfoMetD["other"][i].upper() + "\t" + SNPInfoMetD["eff"][i].upper() + "\t" + SNPInfoMetD["beta"][i] + "\t" + SNPInfoMetD["SE"][i] + "\t" + SNPInfoMetD["frequency"][i] + "\t" + SNPInfoMetD["pValue"][i]                
    outfile.write(s + "\n")
outfile.close()


