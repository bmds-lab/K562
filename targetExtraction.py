
################################
##                            ##
## File: targetExtraction.py  ##
## Author: Dimitri Perrin     ##
##                            ##
################################


# Purpose:
#  - reads the allele-specific sequences
#  - extracts CRISPR sites
#


from subprocess import call
from time import localtime, strftime
from sys import argv
from ast import literal_eval
import re
import string
import ast


#############################
##   Auxiliary functions   ##
#############################


# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

# Function that replaces U with T in the sequence (to go back from RNA to DNA)
def transToDNA(rna):
    switch_UT = string.maketrans('U', 'T')
    dna = rna.translate(switch_UT)
    return dna

# Function that calculates the AT% of a given sequence
def AT_percentage(seq):
    total = 0.0
    length = float(len(seq))
    for c in seq:
        if c in "AT":
            total += 1
    return 100*total/length







if len(argv)!=4:
    print "\nUsage: "+argv[0]+" <dir> <sequenceFile> <bowtie_ref>"
    quit()

dir_ = argv[1]
name = argv[2]
sequenceFile = dir_+name
bowtie_ref = dir_+argv[3]




# Defining the patterns used to detect sequences
pattern_forward = r"(?=([ACG][ACGT]{19}[ACGT]GG))"
pattern_reverse = r"(?=(CC[ACGT][ACGT]{19}[TGC]))"


# Temporary files used in the method
tempTargetFile = dir_+"reads_"+name
alignmentFile = dir_+"alignedRead_"+name


# Parameters
nb_threads_Bowtie = "16"
nb_threads_C = "16"

target_limit = 1000000000 # JUST FOR TESTING PURPOSES
logFile = open("temp_log_targetExtraction_"+name,'w')

outputFile = dir_+"potentialTargets_"+name
out_RNAfold = dir_+"RNAfold_output_"+name
out_targetsToScore = dir_+"targetsToScore_"+name
in_targetScores = dir_+"targets_scored_"+name
C_program = "./findMismatches_threads"
offTargetSites = dir_+"offtargetSites.txt"

accepted_targets = dir_+"accepted_targets_"+name
rejected_targets = dir_+"rejected_targets_"+name

# Defining the patterns used for secondary structures
pattern_RNAstructure = r".{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)"
pattern_RNAenergy = r"\s\((.+)\)"


# Thresholds used when processing secondary structures
low_energy_threshold = -30
high_energy_threshold = -18


# Threshold when looking at the off-target sites
offtarget_threshold = 75


# guide RNA
guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"



###################################
##   Processing the input file   ##
###################################

output = strftime("%H:%M:%S", localtime())+":\tGetting ready to process "+sequenceFile
print output
logFile.write(output+"\n")

inFile = open(sequenceFile,'r')

possibleTargets=dict()
removedTargets=dict()

# For every line in the input file

for line in inFile:
    
    tempArray = line.rstrip().split("\t")
    chr = tempArray[0]
    pos = tempArray[1]
    temp_seq1 = tempArray[2]
    seq1 = temp_seq1[:-2] # the sequence itself stops two characters before the end
    allele_seq1 = temp_seq1[-1:] # the allele is the last character
    temp_seq2 = tempArray[3]
    seq2 = temp_seq2[:-2] # the sequence itself stops two characters before the end
    allele_seq2 = temp_seq2[-1:] # the allele is the last character
    
    # we parse the line and look for forward sequences in seq1
    for match_exon in re.finditer(pattern_forward,seq1):
        start = match_exon.start() # start of the target
        target23 = seq1[start:start+23] # target sequence
        match_pos = str(start-22) # position of the target relative to the start of the difference ref/alt
        if target23 in possibleTargets:
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq1+"_fw")
        else:
            possibleTargets[target23]=[]
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq1+"_fw")

    # we parse the line and look for reverse sequences in seq1
    for match_exon in re.finditer(pattern_reverse,seq1):
        start = match_exon.start() # start of the target
        target23 = rc(seq1[start:start+23]) # target sequence
        match_pos = str(start-22) # position of the target relative to the start of the difference ref/alt
        if target23 in possibleTargets:
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq1+"_rv")
        else:
            possibleTargets[target23]=[]
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq1+"_rv")

    # we parse the line and look for forward sequences in seq2
    for match_exon in re.finditer(pattern_forward,seq2):
        start = match_exon.start() # start of the target
        target23 = seq2[start:start+23] # target sequence
        match_pos = str(start-22) # position of the target relative to the start of the difference ref/alt
        if target23 in possibleTargets:
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq2+"_fw")
        else:
            possibleTargets[target23]=[]
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq2+"_fw")

    # we parse the line and look for reverse sequences in seq2
    for match_exon in re.finditer(pattern_reverse,seq2):
        start = match_exon.start() # start of the target
        target23 = rc(seq2[start:start+23]) # target sequence
        match_pos = str(start-22) # position of the target relative to the start of the difference ref/alt
        if target23 in possibleTargets:
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq2+"_rv")
        else:
            possibleTargets[target23]=[]
            possibleTargets[target23].append(chr+":"+pos+"_"+match_pos+"_"+allele_seq2+"_rv")


    if len(possibleTargets) > target_limit:
        break

inFile.close()
#print "\t\tSKIPPED THE IDENTIFICATION OF REVERSE SEQUENCES!"

output = "\n"+strftime("%H:%M:%S", localtime())+":\t%d potential targets have been identified." % (len(possibleTargets))
print output
logFile.write(output+"\n")


##############################################################
##   Removing targets that have multiple matches in exons   ##
##############################################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that have been observed more than once."
print output
logFile.write(output+"\n")

targetsToRemove=[]
for target23 in possibleTargets:
    # number of occurrences of the target
    total_occurrences = len(possibleTargets[target23])
    
    # number of occurrences of the reverse complement target
    reverse_target23 = rc(target23)
    reverse_also_exists = False
    if reverse_target23 in possibleTargets:
        total_occurrences += len(possibleTargets[reverse_target23])
        reverse_also_exists = True
    
    # we reject if the total is greater than 1
    if total_occurrences>1:
        targetsToRemove.append(target23)
        # we also reject the reverse complement if it exists
        if reverse_also_exists:
            targetsToRemove.append(reverse_target23)

for target23 in targetsToRemove:
    # if the target is not already removed (as reverse-complement of another one)...
    if target23 in possibleTargets:
        # ... then we remove it
        del possibleTargets[target23]
        removedTargets[target23] = "Multiple matches in exons"

output = "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))
print output
logFile.write(output+"\n")


############################################
##   Removing targets that contain TTTT   ##
############################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that contain TTTT."
print output
logFile.write(output+"\n")

targetsToRemove=[]
for target23 in possibleTargets:
    if "TTTT" in target23:
        targetsToRemove.append(target23)

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "TTTT"


output = "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))
print output
logFile.write(output+"\n")



#######################################################
##   Removing targets that have AT% < 20% or > 80%   ##
#######################################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tRemoving all targets that have extreme AT%."
print output
logFile.write(output+"\n")

targetsToRemove=[]
for target23 in possibleTargets:
    target = target23[0:20]
    ATpc = AT_percentage(target)
    if ATpc<20 or ATpc>80:
        targetsToRemove.append(target23)

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "AT%"


output = "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))
print output
logFile.write(output+"\n")



###############################################
##   Using Bowtie to find multiple matches   ##
###############################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tPreparing file for Bowtie analysis."
print output
logFile.write(output+"\n")

outFile = open(tempTargetFile,'w')

tempTargetDict_offset = dict()
for target23 in possibleTargets:
    similarTargets = [target23[0:20]+"AGG", target23[0:20]+"CGG", target23[0:20]+"GGG", target23[0:20]+"TGG", target23[0:20]+"AAG", target23[0:20]+"CAG", target23[0:20]+"GAG", target23[0:20]+"TAG"]
    for seq in similarTargets:
        outFile.write(seq+"\n")
        tempTargetDict_offset[seq] = target23
outFile.close()



output = "\n"+strftime("%H:%M:%S", localtime())+":\tFile ready. Calling Bowtie."
print output
logFile.write(output+"\n")

cmd = "bowtie2 -x "+bowtie_ref+" -p "+nb_threads_Bowtie+" --reorder --no-hd -t -r -U "+tempTargetFile+" -S "+alignmentFile
call([cmd],shell=True)



output = "\n"+strftime("%H:%M:%S", localtime())+":\tStarting to process the Bowtie results."
print output
logFile.write(output+"\n")

inFile = open(alignmentFile,'r')
bowtieLines = inFile.readlines()
inFile.close()

targetsToRemove=[]

i=0
pc_step = 0.1
nb_lines = len(bowtieLines)

while i<nb_lines:
    if i>pc_step*nb_lines:
        output = strftime("%H:%M:%S", localtime())+":\t\t"+str(pc_step*100)+"%"
        print output
	logFile.write(output+"\n")
	pc_step+=0.1

    nb_occurences = 0
    # we extract the read and use the dictionnary to find the corresponding target
    read = bowtieLines[i].rstrip().split("\t")[9]
    seq = ""
    if read in tempTargetDict_offset:
        seq = tempTargetDict_offset[read]
    elif rc(read) in tempTargetDict_offset:
        seq = tempTargetDict_offset[rc(read)]
    else:
        output = "Problem? "+read
	print output
	logFile.write(output+"\n")

    
    
    # we count how many of the eight reads for this target have a perfect alignment
    for j in range(i,i+8):
        if "XM:i:0" in bowtieLines[j]:
            nb_occurences += 1
            # we also check whether this perfect alignment also happens elsewhere
            if "XS:i:0"  in bowtieLines[j]:
                nb_occurences += 1

    # if that number is at least two, the target is removed
    if nb_occurences > 1:
        targetsToRemove.append(seq)
    
    # we continue with the next target
    i+=8


# we can remove the dictionnary
del tempTargetDict_offset

for target23 in targetsToRemove:
    # if the target is not already removed (as reverse-complement of another one)...
    if target23 in possibleTargets:
        # ... then we remove it
        del possibleTargets[target23]
        removedTargets[target23] = "Multiple matches in genome"

output = "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))
print output
logFile.write(output+"\n")


##########################################
##   Calculating secondary structures   ##
##########################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tCalculating secondary structures."
print output
logFile.write(output+"\n")

# WARNING: removing any existing version of the RNAfold output file.
call(["rm -f "+out_RNAfold],shell=True)

# Calling RNAfold on all targets

temp_counter = 0
for target23 in possibleTargets:
    target = "G"+target23[1:20]
    structure = target+guide
    cmd = "echo "+structure+" | RNAfold --noPS >> "+out_RNAfold
    call([cmd],shell=True)
    temp_counter+=1
    if (temp_counter%10000)==0:
        output = strftime("%H:%M:%S", localtime())+":\t\t"+str(temp_counter)+" targets processed."
	print output
	logFile.write(output+"\n")

#    if temp_counter>break_point:
#        print strftime("%H:%M:%S", localtime())+":\t\tReaching a break point!"
#        break

total_number_structures = temp_counter


#########################################
##   Processing secondary structures   ##
#########################################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tProcessing secondary structures."
print output
logFile.write(output+"\n")

inFile = open(out_RNAfold,'r')
RNA_structures = inFile.readlines()
inFile.close()

targetsToRemove=[]
i=0
for target23 in possibleTargets:
    L1 = RNA_structures[2*i].rstrip()
    L2 = RNA_structures[2*i+1].rstrip()
    target = L1[:20]
    if transToDNA(target) != target23[0:20] and transToDNA("C"+target[1:]) != target23[0:20] and transToDNA("A"+target[1:]) != target23[0:20]:
        output = "Error? "+target23+"\t"+target
        print output
	logFile.write(output+"\n")
	quit()
    #    print L1
    #    print target
    #    print L2
    match_structure = re.search(pattern_RNAstructure,L2)
    if match_structure:
        # The structure is correct, we only reject if the energy is too low
        energy = ast.literal_eval(match_structure.group(1))
        if energy < low_energy_threshold:
            targetsToRemove.append(transToDNA(target23))
    else:
        match_energy = re.search(pattern_RNAenergy,L2)
        if match_energy:
            # The structure is not correct, we only reject if the energy is not high enough
            energy = ast.literal_eval(match_energy.group(1))
            if energy <= high_energy_threshold:
                targetsToRemove.append(transToDNA(target23))
    i+=1
#    if i>break_point-1:
#        print strftime("%H:%M:%S", localtime())+":\t\tReaching a break point!"
#        break

for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "Secondary structure or energy"


output = "\t\t%d potential targets are selected for the next step." % (len(possibleTargets))
print output
logFile.write(output+"\n")



#########################
##   Scoring targets   ##
#########################

output = "\n"+strftime("%H:%M:%S", localtime())+":\tScoring targets."
print output
logFile.write(output+"\n")

i=0
targetsToRemove=[]

output = "\n"+strftime("%H:%M:%S", localtime())+":\tWriting targets to file."
print output
logFile.write(output+"\n")

outFile = open(out_targetsToScore,'w')

temp_counter = 0
for target23 in possibleTargets:
    #    print strftime("%H:%M:%S", localtime())+":\t\t"+target23
    target = target23[0:20]
    outFile.write(target+"\n")
    temp_counter += 1
outFile.close()

total_number_scores = temp_counter

output = "\n"+strftime("%H:%M:%S", localtime())+":\tCalling C program.\n"
print output
logFile.write(output+"\n")

cmd = C_program+" "+nb_threads_C+" "+out_targetsToScore+" "+offTargetSites+" "+in_targetScores+" "+str(offtarget_threshold)
call([cmd],shell=True)

output = "\n"+strftime("%H:%M:%S", localtime())+":\tReading the results from file, and processing."
print output
logFile.write(output+"\n")

inFile = open(in_targetScores,'r')
scores = dict()

for target23 in possibleTargets:
    line = inFile.readline().rstrip()
    t20 = line.split("\t")[0]
    score = ast.literal_eval(line.split("\t")[1])
    if t20 != target23[0:20]:
        output = "Problem? "+t20+" - "+target
        print output
	logFile.write(output+"\n")
	quit()
    #    print target23+"\t"+str(score)
    if score <offtarget_threshold:
        targetsToRemove.append(target23)
    else:
        scores[target23] = score

inFile.close()


for target23 in targetsToRemove:
    del possibleTargets[target23]
    removedTargets[target23] = "Off-target score"


output = "\t\t%d potential targets are selected as successful candidates." % (len(possibleTargets))
print output
logFile.write(output+"\n")



output = "\n"+strftime("%H:%M:%S", localtime())+":\tSaving the results."
print output
logFile.write(output+"\n")

outFile = open(accepted_targets,'w')


i=0
for target23 in possibleTargets:
    
    # For each target, we save...
    
    output_line = ""
    target = target23[0:20]
    
    # ... the sequence
    output_line += target+"\t"
    
    # ... the secondary and the energy
    while i<total_number_structures:
        L1 = RNA_structures[2*i].rstrip()
        L2 = RNA_structures[2*i+1].rstrip()
        target_RNA = L1[:20]
        if transToDNA(target_RNA) == target or transToDNA("C"+target_RNA[1:]) == target or transToDNA("A"+target_RNA[1:]) == target:
            structure = L2.split(" ")[0]
            energy = L2.split(" ")[1][1:-1]
            output_line += L1+"\t"+structure+"\t"+energy+"\t"
            break
        i+=1
    if i == total_number_structures:
        print "Error? "+target+" not found in "+out_RNAfold
        quit()
    
    # ... the off-target score
    output_line += str(scores[target23])+"\t"

    # ... the position
    output_line += possibleTargets[target23][0]+"\n"

    outFile.write(output_line)

outFile.close()

outFile = open(rejected_targets,'w')
for target23 in removedTargets:
    outFile.write(target23+"\t"+removedTargets[target23]+"\n")
outFile.close()

output = "\n"+strftime("%H:%M:%S", localtime())+":\tDone."
print output
logFile.write(output+"\n")

logFile.close()

#######################
##                   ##
##    End of File    ##
##                   ##
#######################
