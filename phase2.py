#!/usr/bin/env python

import os
import sys, string
import random
import dendropy
from dendropy import treesim

# generate gene tree from shared species tree
def generate_gene_tree(species_tree, birth_rate, death_rate, num_trees):
    print "Step 2: Generating gene trees based on species tree"
    gene_trees = []
    for i in range(num_trees):
        tree = treesim.birth_death(birth_rate=0.01, # temp values for rates
                                   death_rate=0.02,
                                   tree=species_tree,
                                   max_time=1,
                                   assign_taxa=False)
        gene_trees.append(tree)
    return gene_trees

# simulate sequences from specified gene tree
def simulate_sequences(gene_trees, root_file, species_tree_file):
    print "Step 3: Generating sequences based on gene trees"
    """Simulate sequences from set of gene trees using iSG"""
    c = 1
    for i in gene_trees:
        i.is_rooted=True
        # write tree_file in [:root_file] "Label" (tree) format
        infile = species_tree_file + '_gene' + str(c) + '_tree'
        with open(infile, 'w') as fout:
            fout.write('[:' + root_file + str(c) + '] ')
            fout.write('" Label ' + str(c) + '" ')
#            i.write(fout, 'newick', suppress_rooting=True)
            # suppress edge lengths is temporary hack to deal with fact that iSG does not support having a final branch length at the end
            # which is the format dendropy outputs, i.e. (a: 0.1, b: 0.1):0.0;

            s = i.as_string('newick')
#            print s
            s = s[5:] # get rid of rooting
            last_parens = s.rfind(')')
            print s[:last_parens+1]
            fout.write(s[:last_parens+1])
#            fout.write(s)
            fout.write(';\n')
#            print s[5:-6]
#            fout.write(s[5:-6])
#            fout.write(';')
            fout.close()
        c = c + 1

        outfile = infile + '_sequences'

        # simulate sequences with iSG
        print outfile
        print infile
#        ret = os.system('./indel-seq-gen -m HKY --codon_rate 0.2,0.05,0.75 --outfile %s < %s' % (outfile, infile))
        ret = os.system('./indel-seq-gen -m HKY --outfile %s < %s' % (outfile, infile))
        if ret != 0:
            raise OSError, 'indel-seq-gen failed with code %s ' % ret
        

def write_root(root_file, num_trees):
    """Write root sequences in iSG format by pulling from given file with root sequences"""
    print "Step 1: Writing root sequences in indelSeqGen format"
    with open(root_file, 'r') as fin:
        count = 1;
        while count <= num_trees:
            line = fin.readline().split('|')
            seq = line[1].upper()
            l = len(seq) - 1
            if l > 600: # only take sequences that will create at least 200 amino acids = ~180 kmers each
                out = root_file + str(count)
                with open(out, 'w') as fout:
                    fout.write(str(l))
                    fout.write('\n')
                    fout.write('0'*l+'\n')
                    fout.write(seq)
                    fout.close()
                count = count + 1

# assess homolog clustering
def homolog_assessment():
    """Assesses homology clustering using fablast. Future plans include tblastx+mcl and blastp+mcl"""
#    print "Step 4: Assess homology clustering with fablast"
#    ret = os.system('fablast -k 15 -i')

# assess reconciliation
#def reconciliation_assessment():

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='phase2 of homology pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree to simulate gene tree evolution on')
    parser.add_argument('-n', '--num_gene_trees', default=10, type=int, help='number of gene trees to simulate from species tree'+ds)
    
    parser.add_argument('-b', '--birth_rate', default=random.gauss(0, 0.1), help='birth rate, if none are specified will be drawn from normal distribution')
    parser.add_argument('-d', '--death_rate', default=random.gauss(0, 0.1), help='death rate, if none are specified will be drawn from normal distribution')

    parser.add_argument('-m', '--sequence_model_type', default='codon', help='model to use for sequence simulation, choices are codon, amino acid, nucleotide')

    parser.add_argument('-r', '--root_seq', help='root sequence file, should have all the root sequences you will intend to use, format should be id|sequence')

    parser.add_argument('-hc', '--homolog_clustering_flag', default=True, help='output homology clustering aggregate statistics')
    parser.add_argument('-rc', '--reconciliation_flag', default=True, help='output tree reconciliation aggregate statistics')

    opts = parser.parse_args()

    species_tree_file = opts.species_tree
    species_tree = dendropy.Tree.get_from_path(species_tree_file, schema="newick", rooted=True)
    num_trees = opts.num_gene_trees
    
    birth_rate = opts.birth_rate
    death_rate = opts.death_rate
    hc_flag = opts.homolog_clustering_flag
    rc_flag = opts.reconciliation_flag
    root_seq = opts.root_seq

    # write root sequences first
    write_root(root_seq, num_trees)

    gene_trees = generate_gene_tree(species_tree, birth_rate, death_rate, num_trees)
    simulate_sequences(gene_trees, root_seq, species_tree_file)
