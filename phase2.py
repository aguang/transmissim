#!/usr/bin/env python

import sys, string
import HyPhy
import random
import dendropy
from dendropy import treesim

# generate gene tree from shared species tree
def generate_gene_tree(species_tree, birth_rate, death_rate, num_periods):
    print(species_tree.as_string('newick'))
    for i in range(num_periods):
        tree = treesim.birth_death(birth_rate=1.0, # temp values for rates
                                   death_rate=0.9,
                                   tree=species_tree,
                                   max_time=1,
                                   assign_taxa=False)
#        tree.randomly_assign_taxa(create_required_taxa=True)
        print(tree.as_string('newick'))
        return tree

# simulate sequences from specified gene tree
def simulate_sequences(gene_trees, root_file, speices_tree_file):
    """Simulate sequences from set of gene trees using iSG"""
    c = 0
    for i in gene_trees:
        # write tree_file in <:root_file> "Label" (tree) format
        out = species_tree_file + '_gene' + str(c) + '_tree'
        with open(out, 'w') as fout:
            fout.write('<:' + root_file + str(c) + '> ')
            fout.write('" Label ' + str(c) + '" ')
            print(i.as_string('newick'))
            fout.write(i.as_string('newick'))
            fout.close()
        c = c + 1

def write_root(root_file, num_trees):
    """Write root sequences in iSG format by pulling from given file with root sequences"""
    # open file
    with open(root_file, 'r') as fin:
        for i in range(num_trees):
            line = fin.readline().split('|')
            j = line[0]
            seq = line[1]
            out = root_file + j
            with open(out, 'w') as fout:
                l = len(seq)
                fout.write(str(l))
                fout.write('\n')
                fout.write('0'*l+'\n')
                fout.write(seq)
                fout.close()

# assess homolog clustering
#def homolog_assessment():

# assess reconciliation
#def reconciliation_assessment():

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='phase2 of homology pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree to simulate gene tree evolution on')
    parser.add_argument('-n', '--num_gene_trees', default=10, type=int, help='number of gene trees to simulate from species tree'+ds)
    
    parser.add_argument('-b', '--birth_rate', default=random.gauss(0, 1), help='birth rate, if none are specified will be drawn from normal distribution')
    parser.add_argument('-d', '--death_rate', default=random.gauss(0, 1), help='death rate, if none are specified will be drawn from normal distribution')

    parser.add_argument('-m', '--sequence_model_type', default='codon', help='model to use for sequence simulation, choices are codon, amino acid, nucleotide')

    parser.add_argument('-r', '--root_seq', help='root sequence file, should have all the root sequences you will intend to use, format should be id|sequence')

    parser.add_argument('-hc', '--homolog_clustering_flag', default=True, help='output homology clustering aggregate statistics')
    parser.add_argument('-rc', '--reconciliation_flag', default=True, help='output tree reconciliation aggregate statistics')

    opts = parser.parse_args()

    species_tree_file = opts.species_tree
    species_tree = dendropy.Tree.get_from_path(species_tree_file, schema="newick")
    num_trees = opts.num_gene_trees
    
    birth_rate = opts.birth_rate
    death_rate = opts.death_rate
    hc_flag = opts.homolog_clustering_flag
    rc_flag = opts.reconciliation_flag
    root_seq = opts.root_seq

    # write root sequences first
#    write_root(root_seq, num_trees)
    tree = [generate_gene_tree(species_tree, birth_rate, death_rate, 2)]
    simulate_sequences(tree, root_seq, species_tree_file)
