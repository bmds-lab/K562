
###############################
##                           ##
## File: prepareRefInput.py  ##
## Author: Dimitri Perrin    ##
##                           ##
###############################


# Inputs:
#  - one .fa file containing the whole reference genome (whith line breaks within chromosome)
#
# Outputs:
#  - one .fa file with the whole genome (without those lines breaks)


from sys import argv


if len(argv)!=4:
    print "\nUsage: "+argv[0]+" <dir> <input.fa> <output.fa>"
    quit()

dir_ = argv[1]
in_ = dir_+argv[2]
out_ = dir_+argv[3]

inFile = open(in_,'r')
outFile = open(out_,'w')

nb=0
for line in inFile:
    if line[0]==">":
        print "'"+line[1:].rstrip().split(" ")[0]+"': "+str(nb)
        if nb>0:
            outFile.write("\n")
        nb+=1
    else:
        outFile.write(line.rstrip())

inFile.close()
outFile.close()



#######################
##                   ##
##    End of File    ##
##                   ##
#######################
