#!/usr/bin/env python
import os
import sys, string
import random as rand
from ete2 import Tree
from genesequence import simulate_sequences
import rawreads    

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
        intermediates = options[1]
        gene_sequences = options[3]
        rsrs_path = options[5]
        reads_out = options[6]
        num_genes = int(options[7])
        read_err = options[8]
        explv_flag = options[9]

        if analysis == "all" or "GT->GS":
            os.system("simphy -RS 1 -RL F:200 -SB l:-14,1 -SL F:15 -SP F:1000000000 -SG F:1 -V 0 -OM 1 -O /users/bguang/scratch/homology-testing/tests/stepanalysis/trees -OD 1 -OP 1 -OC 1 -ON 1 -CS 22")
            for i in range(1,201):
                j = str(i).zfill(3)
                gene_tree = Tree("/users/bguang/scratch/homology-testing/tests/stepanalysis/trees/1/g_trees%s.trees" %(j))
                simulate_sequences(gene_tree, "long_roots/root_seqs.txt%i" % (i), "/users/bguang/scratch/homology-testing/tests/stepanalysis/seqs/tree%s" % (j))
            gene_sequences = "/users/bguang/scratch/homology-testing/tests/stepanalysis/seqs/tree"
            print "gene trees are in /users/bguang/scratch/homology-testing/tests/stepanalysis/trees"
            print "gene sequences are in ", gene_sequences
        if analysis == "all" or "GT->GS":
            rawreads.main(gene_sequences, read_err, reads_out, num_genes, rsrs_path, explv_flag)
            print "raw reads are in ", reads_out
