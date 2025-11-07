"""
                Mutation script for G4s harboring more than one SNP
For this Script we need bed files with sequences of G4s harboring more than one SNP. The coordinates of these G4s can be extracted using the bedtools intersect count option - with a count greater than 1 indicating more than one SNP
Also a bed file with the corresponding SNPs is needed. Make again sure, that both files are sorted in a "natural order".
Note that the corresponding SNPs were extracted from the MutatedG4s.bed files - which also harbor the coordinates of the corresponding G4s.
For the pqsFinder predicted G4s a max of two SNPs per G4 was found. For G4Hunter one sequence even harbored three SNPs.
"""

#written by Victoria Ellmer

def mutateSequenceNew(chromG4,startG4,endG4,sequence,chromSNP,posSNP,refSNP,snpSNP):
    if chromG4 == chromSNP and startG4 <= posSNP and endG4 >= posSNP:
        diff = posSNP - startG4 #lets get the position of our SNP respective to the G4
        sequence = sequence.upper() # we need to be sure that all characters are upper case! 
        if refSNP == sequence[diff] or snpSNP == sequence[diff]: 
            refSequence = sequence[:diff] + refSNP + sequence[diff+1:]	
            mutatedSequence = sequence[:diff] + snpSNP + sequence[diff + 1:]
            return(refSequence,mutatedSequence)
        else:
            print("The given base in the G4 sequence is not the annotated as effect or non - effect allele.")
    else:
        print("The positions do not match. Are you sure you ordered the files beforehand?")


### We want to know what double SNPs do in a G4 sequences ###
#Lets start with G4s predicted by pqsFinder:
chromG4 = []
startG4 = []
endG4 = []
seqG4 = []
with open("More_than_one_SNP_PQSFinder.bed", "r+") as fasta:
    for line in fasta.readlines():
        chro = line.split()[0] #get chromosome
        info = line.split()[3]
        start = line.split()[1]
        end = line.split()[2]
        seqG4.append(line.split()[4])
        chromG4.append(chro)
        startG4.append(int(start))
        endG4.append(int(end))
            

G4InfoMetD = {
    "chrom":chromG4,
    "start":startG4,
    "end":endG4,
    "seq":seqG4
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

with open("count_PQSFinder_SNPs.bed", "r+") as fasta:
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


chrom = 0 
start = 0 
stop = 0 
seq = str()
eff = str()
oth = str()
ID = str()
mutseq = str()
mut2seq = []
seq2 = []
        
#mutate the sequence two times
for i in range(0,len(chromG4)):
    #if two G4s have the same coordinates, the sequence is mutated again
    if chrom == chromG4[i] and start == startG4[i] and stop == endG4[i]:
        seq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[0]
        mutseq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], mutseq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[1]
        mut2seq.append(mutseq) #at this point the sequence has two mutations --> so we can save
        seq2.append(seq)
        print(i)
    else:
        chrom = chromG4[i]
        start = startG4[i]
        stop = endG4[i]
        seq = seqG4[i]
        eff = effSNP[i]
        oth = otherSNP[i]
        ID = infoSNP[i]
        seq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[0]
        mutseq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[1]
        #mut2seq.append(seq)

    
#save information 
chromG4new = []
startG4new = []
endG4new = []
rsID1 = []
rsID2 = []
posSNP1 = []
posSNP2 = []
seqG4new = []
other1 = []
effect1 = []
beta_value1 = []
other2 = []
effect2 = []
beta_value2 = []       
chrom = 0 
start = 0 
stop = 0 

for i in range(0,len(chromG4)):
    
    if chrom == chromG4[i] and start == startG4[i] and stop == endG4[i]:
        posSNP2.append(posSNP[i])
        effect2.append(effSNP[i])
        other2.append(otherSNP[i])
        beta_value1.append(beta_value[i])
        rsID2.append(infoSNP[i])
        print(i)
    else:
        chrom = chromG4[i]
        start = startG4[i]
        stop = endG4[i]
        seq = seqG4[i]
        eff1 = effSNP[i]
        oth1 = otherSNP[i]
        ID1 = infoSNP[i]
        pos1 = posSNP[i]
        
        chromG4new.append(chrom)
        startG4new.append(start)
        endG4new.append(stop)
        beta_value2.append(beta_value[i])
        effect1.append(eff1)
        other1.append(oth1)
        rsID1.append(ID1)
        seqG4new.append(seq2)
        posSNP1.append(pos1)

#save file:

outfile = open("Double_SNPs_pqsFinder.bed","a")
#we make a BED file with the following settings: chromosome, start, stop, mutated sequence
s = "chr" + "\t" + "start" + "\t" + "stop" + "\t" + "seq" + "\t" + "mutseq" + "\t" + "posSNP1" + "\t" + "posSNP2" + "\t" + "rsID1" + "\t" + "rsID2" + "\t" + "beta_value1" + "\t" + "beta_value2" + "\n"
outfile.write(s)
for i in range(0,len(chromG4new)):
    s = chromG4new[i] + "\t" + str(startG4new[i]) + "\t" + str(endG4new[i]) + "\t" + seq2[i] + "\t" + mut2seq[i] + "\t" + str(posSNP1[i]) + "\t" + str(posSNP2[i]) + "\t" + rsID1[i] + "\t" + rsID2[i] + "\t" + beta_value1[i] + "\t" + beta_value2[i]
    outfile.write(s + "\n")
outfile.close()
    
        
### Do the same for G4Hunter:
    

### We want to know what double SNPs do in a G4 sequences ###
#Lets start with G4s predicted by pqsFinder:
chromG4 = []
startG4 = []
endG4 = []
seqG4 = []
with open("More_than_one_SNP_G4Hunter.bed", "r+") as fasta:
    for line in fasta.readlines():
        chro = line.split()[0] #get chromosome
        info = line.split()[3]
        start = line.split()[1]
        end = line.split()[2]
        seqG4.append(line.split()[4])
        chromG4.append(chro)
        startG4.append(int(start))
        endG4.append(int(end))
            

G4InfoMetD = {
    "chrom":chromG4,
    "start":startG4,
    "end":endG4,
    "seq":seqG4
    }

#load data of SNPs overlapped with G4Hunter

chromSNP = []
posSNP = []
infoSNP = []
otherSNP = []
effSNP = []  
beta_value =[]
standard_error =[]
frequency = []
p_value = []

with open("count_G4Hunter_SNPs.bed", "r+") as fasta:
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


chrom = 0 
start = 0 
stop = 0 
seq = str()
eff = str()
oth = str()
ID = str()
mutseq = str()
mut2seq = []
seq2 = []

counter = [] #add this counter to know which G4 harbors how many SNPS
count = 0
#first append endG4, because else we will run out of indices later on *
endG4.append("")
#mutate the sequence two or three times (depending on how many SNPs occur in the G4)
for i in range(0,len(chromG4)):
    count +=1
    
    if chrom == chromG4[i] and start == startG4[i] and stop == endG4[i]:
        seq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[0]
        mutseq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], mutseq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[1]
        if (endG4[i] != endG4[i+1]) or i == len(chromG4): #*
            mut2seq.append(mutseq) #at this point the sequence has all the mutations --> so we can save
            seq2.append(seq)
            print(i)
            counter.append(count)
            count = 0
    else:
       
        chrom = chromG4[i]
        start = startG4[i]
        stop = endG4[i]
        seq = seqG4[i]
        eff = effSNP[i]
        oth = otherSNP[i]
        ID = infoSNP[i]
        seq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[0]
        mutseq = mutateSequenceNew(chromG4[i], startG4[i], endG4[i], seq, chromG4[i], posSNP[i], otherSNP[i], effSNP[i])[1]
        #mut2seq.append(seq)

#save information 
chromG4new = []
startG4new = []
endG4new = []

seqG4new = []
       
chrom = 0 
start = 0 
stop = 0 

for i in range(0,len(chromG4)):
    
    if chrom == chromG4[i] and start == startG4[i] and stop == endG4[i]:
      
        print(i)
    else:
        chrom = chromG4[i]
        start = startG4[i]
        stop = endG4[i]
        seq = seqG4[i]
        
        chromG4new.append(chrom)
        startG4new.append(start)
        endG4new.append(stop)
        seqG4new.append(seq2)
        
doubleSNPs = []
a = 0
for k in counter:
    doubleSNPs.append(infoSNP[a:a+k]) #extract the SNP IDs for the corresponding G4s
    a+=k


#

#save file:

outfile = open("Double_SNPs_G4Hunter.bed","a")
#we make a BED file with the following settings: chromosome, start, stop, mutated sequence
s = "chr" + "\t" + "start" + "\t" + "stop" + "\t" + "seq" + "\t" + "mutseq" + "\t" + "SNP1" + "\t" + "SNP2" + "\t" + "SNP3" + "\n"
outfile.write(s)
for i in range(0,len(chromG4new)):
    s = chromG4new[i] + "\t" + str(startG4new[i]) + "\t" + str(endG4new[i]) + "\t" + seq2[i] + "\t" + mut2seq[i] + "\t" 
    #we now need to add the SNPs if we have any 
    if counter[i] == 2:
        s = s + doubleSNPs[i][0] + "\t" + doubleSNPs[i][1] + "\t" + "NA"
    if counter[i] == 3:
        s = s + doubleSNPs[i][0] + "\t" + doubleSNPs[i][1] + "\t" + doubleSNPs[i][2]
    outfile.write(s + "\n")
outfile.close()





