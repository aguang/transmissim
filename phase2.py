#!/usr/bin/env python
import os
import sys, string
import phase1
import random
from ete2 import Tree

# generate gene tree from shared species tree
def filter_tree(tree):
    tree.resolve_polytomy(recursive=True)
    D = tree.search_nodes(name="D")
    for death in D:
        P = death.up
        if P.is_root():
            death.detach()
            children = P.children
            if children:
                child = P.children[0] #assumes bifurcating tree
                child.dist = P.dist + child.dist
                P.delete()
            else:
                return 0
        else:
            G = P.up
            G.dist=G.dist+P.dist
            death.detach()
            P.delete()
    return tree

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='phase2 of homology pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree to simulate gene tree evolution on')
    parser.add_argument('-n', '--num_gene_trees', default=10, type=int, help='number of gene trees to simulate from species tree'+ds)
    
    parser.add_argument('-b', '--birth_rate', default=random.gauss(0, 0.1), help='birth rate, if none are specified will be drawn from normal distribution')
    parser.add_argument('-d', '--death_rate', default=random.gauss(0, 0.1), help='death rate, if none are specified will be drawn from normal distribution')

    parser.add_argument('-w', '--wr_flag', default=0, help='flag for whether root sequences need to be written or not')
    parser.add_argument('-r', '--root_seq', help='root sequence file, should have all the root sequences you will intend to use, format should be id|sequence')
    parser.add_argument('-dir', '--out_dir', help='out directory')
    opts = parser.parse_args()

    species_tree_file = opts.species_tree
    num_trees = int(opts.num_gene_trees)
    
    birth_rate = float(opts.birth_rate)
    death_rate = float(opts.death_rate)

    wr_flag = int(opts.wr_flag)
    root_seq = opts.root_seq
    out_dir = opts.out_dir

    if wr_flag == 1:
        phase1.write_root(root_seq, num_trees)

#    out = out_dir + species_tree_file
#    generate_gene_tree(species_tree_file, birth_rate, death_rate, num_trees, out)
    for i in range(0,num_trees):
        gene_tree_file = species_tree_file + str(i)
        gene_tree = Tree(gene_tree_file, format=1)
        gene_tree = filter_tree(gene_tree)
        if gene_tree == 0:
            continue
        out = out_dir + "tree" + str(i)
        r = root_seq + str(i+1)
        phase1.simulate_sequences(gene_tree, r, gene_tree_file, out)
