#!/usr/python
# AUGUST GUANG
# generate reads from gene sequences

import sys, getopt, string
import numpy as np
import random
import itertools
import os

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile, numGenes):
    speciesDict = {}
    geneID = 1
    for i in range(0,numGenes):
        #f = iFile + str(i) + '_sequences.seq'
        f = iFile
        try:
            with open(f, 'r') as fin:
                # read lines in chunks of n
                # format will be:
            
                #   >species_ID
                #   ACTGACTGNNNNNNN
                for key, group in itertools.groupby(fin, lambda line: line.startswith('>')):
                    if key:
                        species = next(group).strip()
#                        header = next(group).strip().split('_')
#                        species = header[0]
#                        geneCopy = header[1]
#                        individual = header[2]
                    else:
                        lines=''.join(group).strip()
#                        geneName = species + '-gene' + str(geneID) + '-copy' + geneCopy
                        geneName = species
                        print geneName
                 
                        if species in speciesDict:
                            speciesDict[species].append((geneName, lines))
                        else:
                            speciesDict[species] = [(geneName, lines)]
            geneID = geneID + 1
        except:
            print "error in creating species dictionary for gene family", geneID
    return speciesDict

# writes all genes for each species with accompanying species ID into a reference.fa file along with
# a BED file for RNAseqreadsimulator
# then writes as a fastq CASAVA file
def writeSpecies(directory, speciesDict):
    for speciesKey in speciesDict:
        # establish file paths
        outr = os.path.join(directory, speciesKey[1:])
        outf = outr + '.fa'
        bed = outr + '.bed'
        explv = outr + '_explv.txt'

        # write BED and fasta file
        fa = open(outf, 'w')
        fa.write(">chr1\n")
        outb = open(bed, 'w')
        pos = 0
        for seq in speciesDict[speciesKey]:
            # write BED
            l = len(seq[1])
            end = pos + l
            outb.write("chr1\t%i\t%i\t%s\t0\t+\t%i\t%i\t0\t1\t%i\t0\n" % (pos, end, seq[0][1:],pos,pos,l))
            pos = end + 1000

            # write fasta
            fa.write(seq[1])
        fa.close()
        outb.close()

        # run RNAseqreadsimulator
        path = os.path.dirname(sys.argv[0])
        path = os.path.abspath(path)
        posbias = os.path.join(path, 'posbias.txt')
        readerr = os.path.join(path, 'readerror.txt')
        ret = os.system('python %s/genexplvprofile.py %s > %s' % (path, bed, explv))
        ret = os.system('python %s/gensimreads.py -e %s -l 100 -n 1000000 -p 200,20 %s | python %s/getseqfrombed.py -l 100 - %s | python splitfasta.py -o %s' % (path, explv, bed, path, outf, outr))

        # turn fasta files into fastq files with highest possible Phred quality scores
        # and CASAVA headers
        lane = random.randint(1,3)
        tile = random.randint(1, 2000)
        make_headers(outr, 1, lane, tile)
        make_headers(outr, 2, lane, tile)

# make CASAVA headers
def make_headers(read, paired_num, lane, tile):
    fasta = read + '_' + str(paired_num) + '.fa'
    fastq = read + '_' + str(paired_num) + '.fq'
    qual_range = '0?@ABCDEFGH'
    with open(fasta, 'r') as fa:
        fq = open(fastq, 'w')
        for key, group in itertools.groupby(fa, lambda line: line.startswith('>')):
            if key:
                header = next(group).strip()
            else:
                lines = ''.join(group).strip()
                lines = lines.translate(None, '\n')
                num_bases = len(lines)

                splits = header.rsplit('_',4)
                read_num = splits[2]
                gene_num = splits[0].split('-')[1]

                fq.write('@SIM-12345:1:ABCDEFGHI:%i:%i:%s:%s %i:N:0:ATCACG\n' % (lane, tile, gene_num, read_num, paired_num))
                fq.write(lines+'\n')
                fq.write('+\n')
                qual_indices = np.random.binomial(10, 0.7, num_bases)
                qualities = ''.join(map(lambda x:qual_range[x], qual_indices))
                fq.write(qualities + '\n')        

def main(argv):
    numGenes = 0
    inExt = ''
    directory = ''
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:d:n:", ["ifile="])
    except getopt.error, msg:
        print 'rawreads.py -i <sequenceFile> -n <numGenes> -d <directory>'
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

if __name__ == "__main__":
    main(sys.argv[1:])
