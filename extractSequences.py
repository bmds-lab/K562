
################################
##                            ##
## File: extractSequences.py  ##
## Author: Dimitri Perrin     ##
##                            ##
################################


# Purpose:
#  - reads the vcf file
#  - extracts sequences that may be useful in terms of allele-specific CRISPR sites
#
# What makes a sequence useful?
#  - several bases are changed between REF and ALT
#  - or there is an indel that shiftes ALT out of frame from REF

# Assumptions:
#  - variants are more than one CRISPR target away from each other
#  - this is true for the K562 data (min distance: 1868)

# Inputs:
#
#
# Outputs:
#


from time import localtime, strftime
from sys import argv
from ast import literal_eval


if len(argv)!=5:
    print "\nUsage: "+argv[0]+" <dir> <vcf_file> <ref_genome.fa> <output.txt>"
    quit()

dir_ = argv[1]
vcf_file = dir_+argv[2]
ref_genome = dir_+argv[3]
out_ = dir_+argv[4]

# obtained from prepareRefInput.py
chr_offset = {'chr10': 0, 'chr11': 1, 'chr12': 2, 'chr13': 3, 'chr14': 4, 'chr15': 5, 'chr16': 6, 'chr17': 7, 'chr18': 8, 'chr19': 9, 'chr1': 10, 'chr20': 11, 'chr21': 12, 'chr22': 13, 'chr2': 14, 'chr3': 15, 'chr4': 16, 'chr5': 17, 'chr6': 18, 'chr7': 19, 'chr8': 20, 'chr9': 21, 'chrX': 22, 'chrY': 23}


print strftime("%H:%M:%S", localtime())+": Loading the reference genome."

chr_sequences = []
inFile = open(ref_genome,'r')
for line in inFile:
    chr_sequences.append(line.rstrip())
inFile.close()


print strftime("%H:%M:%S", localtime())+": Starting to read the VCF file."

nb = 0
max = 0

inFile = open(vcf_file,'r')
outFile = open(out_,'w')

for line in inFile:

    if line[0] == "#": # comments are ignored
        continue
    tempArray = line.rstrip().split("\t")

    filter = tempArray[6]
    if filter != "PASS": # in our case, it should always be "PASS"
        continue # but no point doing anything to this one if not passing

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


    if len(ref)==1 and len(alt)==1: # not interested in single point mutation
        continue
    
    if additionalField[1]!="|": # not interested if unphased
        continue
    
    if additionalField[0]==additionalField[2]: # not interested here either
        continue

    tempAlt = alt.split(",")


    # default case
    if len(tempAlt)==1:
        
        # we extract the 'ref' sequence
        ref_seq = chr_sequences[chr_offset[chr]][pos-22:pos+len(ref)+23]
        # we tag it with the first field of the phasing
        ref_seq+="_"+additionalField[0]

        # we extract the 'alt' sequence
        alt_seq = chr_sequences[chr_offset[chr]][pos-22:pos]+tempAlt[0]+chr_sequences[chr_offset[chr]][pos+len(ref):pos+len(ref)+23]
        # we tag it with the second field of the phasing
        alt_seq+="_"+additionalField[2]

        if max<len(ref):
            max=len(ref)
        if max<len(alt):
            max=len(alt)
                
        outFile.write(chr+"\t"+str(pos)+"\t"+ref_seq.upper()+"\t"+alt_seq.upper()+"\n")



    # other case: two 'alt' alleles, so we focus on those and ignore 'ref' (because we have already checked that in this case, the phasing is always 1|2 or 2|1)
    else:
        
        # we extract the first 'alt' sequence
        alt_seq1 = chr_sequences[chr_offset[chr]][pos-22:pos]+tempAlt[0]+chr_sequences[chr_offset[chr]][pos+len(ref):pos+len(ref)+23]
        # we tag it with the first field of the phasing
        alt_seq1+="_"+additionalField[0]


        # we extract the second 'alt' sequence
        alt_seq2 = chr_sequences[chr_offset[chr]][pos-22:pos]+tempAlt[1]+chr_sequences[chr_offset[chr]][pos+len(ref):pos+len(ref)+23]
        # we tag it with the second field of the phasing
        alt_seq2+="_"+additionalField[2]

        if max<len(tempAlt[0]):
            max=len(tempAlt[0])
        if max<len(tempAlt[1]):
            max=len(tempAlt[1])

        outFile.write(chr+"\t"+str(pos)+"\t"+alt_seq1.upper()+"\t"+alt_seq2.upper()+"\n")


    nb+=1


inFile.close()
outFile.close()

print strftime("%H:%M:%S", localtime())+": Done. "+str(nb)+" sequences worth investigating. "+"Maximum allele length: "+str(max)



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
