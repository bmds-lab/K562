
##########################################
##                                      ##
##  File: prepareListOfftargetSites.py  ##
##  Author: Dimitri Perrin              ##
##                                      ##
##########################################

#
# Purpose: identify all offtarget sites in the whole genome
#

#
# Inputs:
#           - one file with all the chromosome sequences, one chromosome per line
#

#
# Output:
#           - one file with all the sites
#


from time import localtime, strftime
import re
import string
from sys import argv
from ast import literal_eval



# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq



if len(argv)!=4:
    print "\nUsage: "+argv[0]+" <dir> <ref_genome> <variants>"
    quit()

dir = argv[1]
sequenceFile = dir+argv[2]
variants = dir+argv[3]
offtargetSites = dir+"offtargetSites.txt"



# The input .fa file has one whole chromosome per line. This gives the order in which they were stored. It is obtained from prepareRefInput.py (which generated that .fa file)
chr_offset = {'chr10': 0, 'chr11': 1, 'chr12': 2, 'chr13': 3, 'chr14': 4, 'chr15': 5, 'chr16': 6, 'chr17': 7, 'chr18': 8, 'chr19': 9, 'chr1': 10, 'chr20': 11, 'chr21': 12, 'chr22': 13, 'chr2': 14, 'chr3': 15, 'chr4': 16, 'chr5': 17, 'chr6': 18, 'chr7': 19, 'chr8': 20, 'chr9': 21, 'chrX': 22, 'chrY': 23}


# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"


print strftime("%H:%M:%S", localtime())+":\tLoading the reference genome sequences"

chr_sequences = []
inFile = open(sequenceFile,'r')
for line in inFile:
    chr_sequences.append(line.rstrip())
inFile.close()

outFile = open(offtargetSites,'w')

print strftime("%H:%M:%S", localtime())+":\tExtracting from reference genome"



# For every chromosome
for chr in range(0,len(chr_sequences)):
    
    line_chr = chr_sequences[chr]
    
    print strftime("%H:%M:%S", localtime())+":\tLine "+str(chr)
    
    # we parse the line and look for forward sequences
    print strftime("%H:%M:%S", localtime())+":\t\tForward-parsing the chromosome."
    match_chr = re.findall(pattern_forward_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing and saving the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each sequence
        for i in range(0,len(match_chr)):
            outFile.write(match_chr[i][0:20]+"\n")
    else:
        print "\t\t\tWe did not detect any possible off-target sites."

    # we parse the line and look for reverse sequences
    print strftime("%H:%M:%S", localtime())+":\t\tReverse-parsing the chromosome."
    match_chr = re.findall(pattern_reverse_offsite,line_chr)
    print strftime("%H:%M:%S", localtime())+":\t\tProcessing the parsed sequences."
    if match_chr:
        print "\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites."
        # we save each reverse-complement sequence
        for i in range(0,len(match_chr)):
            # we reverse-complement the sequence and count the number of mismatches between the rc and the target
            outFile.write(rc(match_chr[i])[0:20]+"\n")

inFile.close()



print strftime("%H:%M:%S", localtime())+":\tExtracting from phased variants"

nb = 0

inFile = open(variants,'r')

for line in inFile:
    
    if line[0] == "#": # comments are ignored
        continue
    tempArray = line.rstrip().split("\t")

    filter = tempArray[6]
    if filter != "PASS": # in our case, it should always be "PASS"
        continue # but no point doing anything to this one if not passing

    #chr = tempArray[0][3:] # ***FOR K562 only*** ignoring the first three characters: chr
    chr = tempArray[0] 
    pos = literal_eval(tempArray[1])-1 # with -1 because VCF POS is 1-based
    ID = tempArray[2] # in our case, always equal to "."
    ref = tempArray[3]
    alt = tempArray[4]
    qual = tempArray[5]
    filter = tempArray[6]
    info = tempArray[7]
    format = tempArray[8]
    additionalField = tempArray[9]



    if chr not in chr_offset:
        print chr+" not in chr_offset\n"+line
	quit()

    ref_seq = chr_sequences[chr_offset[chr]][pos-22:pos+len(ref)+23]
    
    # we parse the ref_seq and look for forward sequences
    match_chr = re.findall(pattern_forward_offsite,ref_seq)
    if match_chr:
        # we save each sequence
        for i in range(0,len(match_chr)):
            outFile.write(match_chr[i][0:20]+"\n")
        nb += len(match_chr)

    # we parse the line and look for reverse sequences
    match_chr = re.findall(pattern_reverse_offsite,ref_seq)
    if match_chr:
        # we save each reverse-complement sequence
        for i in range(0,len(match_chr)):
            # we reverse-complement the sequence and count the number of mismatches between the rc and the target
            outFile.write(rc(match_chr[i])[0:20]+"\n")
        nb += len(match_chr)


    tempArray = alt.split(",")
    for el in tempArray:
        alt_seq = chr_sequences[chr_offset[chr]][pos-22:pos]+el+chr_sequences[chr_offset[chr]][pos+len(ref):pos+len(ref)+23]

        # we parse the ref_seq and look for forward sequences
        match_chr = re.findall(pattern_forward_offsite,alt_seq)
        if match_chr:
            # we save each sequence
            for i in range(0,len(match_chr)):
                outFile.write(match_chr[i][0:20]+"\n")
            nb += len(match_chr)

        # we parse the line and look for reverse sequences
        match_chr = re.findall(pattern_reverse_offsite,alt_seq)
        if match_chr:
            # we save each reverse-complement sequence
            for i in range(0,len(match_chr)):
                # we reverse-complement the sequence and count the number of mismatches between the rc and the target
                outFile.write(rc(match_chr[i])[0:20]+"\n")
            nb += len(match_chr)

inFile.close()
print "\t\t\tWe detected "+str(nb)+" possible off-target sites."


outFile.close()

print "\n"+strftime("%H:%M:%S", localtime())+":\tDone."
