#! /usr/python
# AUGUST GUANG
# JULY 2015: script to generate reads from gene sequences

import sys, getopt, string
import itertools
import os

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile, numGenes):
    speciesDict = {}
    geneID = 0
    for i in range(0,numGenes):
        f = iFile + str(i) + '_sequences.seq'
        try:
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
        except:
            print geneID, " gene family not found, skipping"
            geneID = geneID + 1
    return speciesDict

# writes all genes for each species with accompanying species ID into a reference.fa file along with
# a BED file for RNAseqreadsimulator
def writeSpecies(directory, speciesDict):
    for speciesKey in speciesDict:
        outf = os.path.join(directory, speciesKey[1:] + '.fa')
        bed = os.path.join(directory, speciesKey[1:] + '.bed')
        explv = os.path.join(directory, speciesKey[1:] + '_explv.txt')
        outr = os.path.join(directory, speciesKey[1:])
        fa = open(outf, 'w')
        fa.write(">chr1\n")
        outb = open(bed, 'w')
        pos = 0
        for seq in speciesDict[speciesKey]:
            # write BED
            l = len(seq[1])
            end = pos + l
            outb.write("chr1\t%i\t%i\t%s\t0\t+\t%i\t%i\t0\t1\t%i\t0\n" % (pos, end, seq[0][1:],pos,pos,l))
            pos = end

            # write fasta
            fa.write(seq[1])
        fa.close()
        outb.close()

        path = os.path.dirname(sys.argv[0])
        path = os.path.abspath(path)
        posbias = os.path.join(path, 'posbias.txt')
        readerr = os.path.join(path, 'readerror.txt')
        ret = os.system('python %s/genexplvprofile.py %s > %s' % (path, bed, explv))
        ret = os.system('python %s/gensimreads.py -e %s -b %s -l 100 -p 200,20 %s | python %s/getseqfrombed.py -b %s -f A -r 0.01 -l 100 - %s | python splitfasta.py -o %s' % (path, explv, posbias, bed, path, readerr, outf, outr))

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
