#! /usr/python
# AUGUST GUANG
# MAY 2014: quick tool for concatenating sequences for homologize in agalma pipeline
# designed specifically for how I've set up indelseqgen
# JULY 2014: combined agalma_format and simmod.py to create tool for creating sequences in agalma pipeline format
# added species naming on original tree
# OCTOBER 2015: modified to work with SimPhy naming convention

import sys, getopt, string
import itertools
import pprint #tmp for pretty printing
import os

# separates genes for each species and places them into a species dictionary
def sepGenes(iFile, numGenes):
    speciesDict = {}
    geneID = 1
    for i in range(1,numGenes):
        f = iFile + str(i) + '_sequences.seq'
        try:
            with open(f, 'r') as fin:
                # read lines in chunks of n
                # format will be:
            
                #   >species_ID-geneID
                #   ACTGACTGNNNNNNN
                duplicates = {}
                for key, group in itertools.groupby(fin, lambda line: line.startswith('>')):
                    if key:
                        header = next(group).strip().split('_')
                        species = header[0]
                        geneCopy = header[1]
                        individual = header[2]
                    else:
                        lines=''.join(group).strip()
                        geneName = species + '-gene' + str(geneID) + '-copy' + geneCopy

                        if species in speciesDict:
                            speciesDict[species].append((geneName, lines))
                        else:
                            speciesDict[species] = [(geneName, lines)]
            geneID = geneID + 1
        except:
            print geneID, " gene family not found, skipping"
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

def main(sequences_directory, num_genes, out_directory):
    species_dictionary = sepGenes(sequences_directory, num_genes)
    writeSpecies(out_directory, species_dictionary)
