#! /usr/python
# AUGUST GUANG
# NOVEMBER 2013

# quick tool for separating sequences

import sys, getopt, string
import itertools
import pprint #tmp for pretty printing

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile, oFile, geneNum):
    with open(iFile, 'r') as fin:
        # read lines in chunks of n
        # format will be:

        #   >species_ID
        #   ACTGACTGNNNNNNN
        geneCount = 0
        speciesCount = 1
        outStr = oFile + str(speciesCount) + '.fa'
        f = open(outStr, 'w')
        for key, group in itertools.groupby(fin, lambda line: line.startswith('>')):
            if key:
                header = next(group).strip()
            else:
                lines=''.join(group).strip()
                f.write(header + '\n')
                f.write(lines + '\n')
                geneCount = geneCount + 1
                if geneCount == geneNum:
                    geneCount = 0
                    f.close()
                    speciesCount = speciesCount + 1
                    outStr = oFile + str(speciesCount) + '.fa'
                    f = open(outStr, 'w')

def main(argv):
    iFile = ''
    oFile = ''
    numGenes = 0
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:n:", ["iFile=", "ofile=", "numGenes="])
    except getopt.error, msg:
        print 'simmod.py -i <sequenceFile> -o <outputFiles>'
        sys.exit(2)
    # process arguments
    for opt, arg in opts:
        if opt in ("-i", "--iFile"):
            iFile = arg
        elif opt in ("-o", "--output"):
            oFile = arg
        elif opt in ("-n", "--numGenes"):
            numGenes = arg

    numGenes = int(numGenes)

    speciesDict = sepGenes(iFile, oFile, numGenes)

if __name__ == "__main__":
    main(sys.argv[1:])
