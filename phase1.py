#!/usr/bin/env python

import os
import sys, string
import random
import dendropy

# simulate sequences from specified gene tree
def simulate_sequences(species_tree, root_file, species_tree_file, num_trees):
    print "Sequence Evolution Based on Species Tree"
    """Simulate sequences from set of gene trees using iSG"""
    species_tree.is_rooted=True
    # write tree_file in [:root_file] "Label" (tree) format
    for i in range(num_trees):
        infile = species_tree_file + "iSG" + str(i+1)
        with open(infile, 'w') as fout:
            fout.write('[:' + root_file + str(i+1) + '] ')
            fout.write('" Label ' + '" ')
#            i.write(fout, 'newick', suppress_rooting=True)
            # suppress edge lengths is temporary hack to deal with fact that iSG does not support having a final branch length at the end
            # which is the format dendropy outputs, i.e. (a: 0.1, b: 0.1):0.0;

            s = species_tree.as_string('newick')
            s = s[5:] # get rid of rooting
            last_parens = s.rfind(')')
            fout.write(s[:last_parens+1])
            fout.write(';\n')
            fout.close()

        outfile = infile + '_sequences'

        # simulate sequences with iSG
        ret = os.system('./indel-seq-gen -m HKY --codon_rate 0.2,0.05,0.75 --outfile %s < %s' % (outfile, infile))
#        ret = os.system('./indel-seq-gen -m HKY --outfile %s < %s' % (outfile, infile))
        if ret != 0:
            raise OSError, 'indel-seq-gen failed with code %s ' % ret

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='phase1 of homology pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree to simulate gene tree evolution on')
    parser.add_argument('-n', '--num_trees', default=10, help='number of trees to generate')
    parser.add_argument('-m', '--sequence_model_type', default='codon', help='model to use for sequence simulation, choices are codon, amino acid, nucleotide')
    parser.add_argument('-r', '--root_seq', help='root sequence file, should have all the root sequences you will intend to use, format should be id|sequence')

    opts = parser.parse_args()

    species_tree_file = opts.species_tree
    species_tree = dendropy.Tree.get_from_path(species_tree_file, schema="newick", rooted=True)
    num_trees = opts.num_trees
    root_seq = opts.root_seq

    simulate_sequences(species_tree, root_seq, species_tree_file, num_trees)
