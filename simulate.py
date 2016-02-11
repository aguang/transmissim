#!/usr/bin/env python
import os
import sys, string
import random as rand
import math
from ete2 import Tree
from genesequence import simulate_sequences
import rawreads    
import agalma_format

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params
    if params:
        options = [line.strip().split('//')[0].strip() for line in open(params)]
        
        analysis = options[0]
        input_directory = options[1]
        intermediates = options[2]
        num_genes = int(options[3])
        tree_directory = options[5]
        gene_sequences = options[7]
        rsrs_path = options[9]
        reads_out = options[10]
        read_len = options[11]
        read_err = options[12]
        explv_flag = options[13]

        if analysis == "all" or "ST->GS":
            os.system("mkdir %s" % (tree_directory))
            os.system("mkdir %s" % (gene_sequences))
            os.system("mkdir %s/seqs/" % (gene_sequences))

            os.system("simphy -RS 1 -RL F:%i -SB l:-14,1 -SL F:15 -SP F:1000000000 -SG F:1 -V 0 -OM 1 -O %s -OD 1 -OP 1 -OC 1 -ON 1 -CS 22" % (num_genes, tree_directory))
            for i in range(1,num_genes+1):
                j = str(i).zfill(int(math.log10(num_genes))+1)
                gene_tree = Tree("%s/1/g_trees%s.trees" % (tree_directory, j))
                simulate_sequences(gene_tree, "long_roots/root_seqs.txt%i" % (i), "%s/seqs/tree%s" % (gene_sequences, i))

            if intermediates:
                gs_in = "%s/seqs/tree" % (gene_sequences)
                agalma_format.main(gs_in, num_genes, gene_sequences)
                print "gene sequences are in ", gene_sequences

        if analysis == "all" or "GS->RR":
            os.system("mkdir %s" % (reads_out))
            if analysis == "GS->RR":
                gene_sequences = input_directory
            rawreads.main("%s/seqs/tree" % (gene_sequences), read_err, read_len, reads_out, num_genes, rsrs_path, explv_flag)
            print "raw reads are in ", reads_out
