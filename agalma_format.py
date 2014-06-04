#! /usr/python
# AUGUST GUANG
# MAY 2014

# quick tool for concatenating sequences

import sys, getopt, string
import itertools
import pprint #tmp for pretty printing

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile):

    speciesDict = {}
    geneID = 0

    for f in iFile:
        print f
        with open(f, 'r') as fin:
            print geneID
            # read lines in chunks of n
            # format will be:
            
            #   >species_ID
            #   ACTGACTGNNNNNNN
            for key, group in itertools.groupby(fin, lambda line: line.startswith('>')):
                if key:
                    header = next(group).strip()
                else:
                    lines=''.join(group).strip()
                    geneName = header + '-' + str(geneID)
                    if header in speciesDict:
                        speciesDict[header].append((geneName, lines))
                    else:
                        speciesDict[header] = [(geneName, lines)]
        geneID = geneID + 1
    return speciesDict

def writeSpecies(oFile, speciesDict):
    f = open(oFile, 'w')
    for speciesKey in speciesDict:
        for seq in speciesDict[speciesKey]:
            f.write(seq[0] + '\n')
            f.write(seq[1] + '\n')
    f.close()

def main(argv):
    iFile = ''
    oFile = 'all.fa'
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:n:", ["ifile="])
    except getopt.error, msg:
        print 'simmod.py -i <sequenceFile> -o <outputFiles>'
        sys.exit(2)
    # process arguments
    iFile = []
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            iFile.append(arg)

    speciesDict = sepGenes(iFile)
    writeSpecies(oFile, speciesDict)
#    pprint.pprint(speciesDict)

if __name__ == "__main__":
    main(sys.argv[1:])
