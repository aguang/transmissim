#! /usr/python
# AUGUST GUANG
# MAY 2014: quick tool for concatenating sequences for homologize in agalma pipeline
# designed specifically for how I've set up indelseqgen
# JULY 2014: combined agalma_format and simmod.py to create tool for creating sequences in agalma pipeline format
# added species naming on original tree

import sys, getopt, string
import itertools
import pprint #tmp for pretty printing
import os

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile, numGenes):
    speciesDict = {}
    geneID = 0
    for i in range(1,numGenes+1):
        f = iFile + str(i) + '_sequences.seq'
        with open(f, 'r') as fin:
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

# writes all genes for each species with accompanying species ID
def writeSpecies(directory, speciesDict):
    i = 0
    for speciesKey in speciesDict:
        out = os.path.join(directory, speciesKey[1:] + '.fa')
        f = open(out, 'w')
        for seq in speciesDict[speciesKey]:
            f.write(seq[0] + '\n')
            f.write(seq[1] + '\n')
        f.close()

def main(argv):
    numGenes = 0
    inExt = ''
    directory = ''
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:d:n:", ["ifile="])
    except getopt.error, msg:
        print 'agalma_format.py -i <sequenceFile> -n <numGenes> -d <directory>'
        sys.exit(2)
    # process arguments
    iFile = []
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            # extension for sequence files
            inExt = arg
        if opt in ("-d", "--directory"):
            # directory for out files
            directory = arg
        if opt in ("-n", "--numGenes"):
            numGenes = int(arg)

    speciesDict = sepGenes(inExt, numGenes)
    writeSpecies(directory, speciesDict)
#    pprint.pprint(speciesDict)

if __name__ == "__main__":
    main(sys.argv[1:])
