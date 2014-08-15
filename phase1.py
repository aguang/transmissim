#!/usr/bin/env python
import os
import sys, string
import random as rand
import dendropy

# simulate sequences from specified gene tree
def simulate_sequences(species_tree, root_file, species_tree_file, num_genefams,out_dir):
    print "Step 2: Sequence Evolution Based on Species Tree"

    """Simulate sequences from set of gene trees using iSG. Rates of evolution drawn from random distributions."""
    species_tree.is_rooted=True
    # write tree_file in [:root_file] "Label" (tree) format
    for i in range(num_genefams):
        infile = out_dir + species_tree_file + str(i+1)
        with open(infile, 'w') as fout:
            fout.write('[:' + root_file + str(i+1) + '] ')
            fout.write('" Label ' + '" ')
            s = species_tree.as_string('newick')
            s = s[5:] # get rid of rooting
            last_parens = s.rfind(')')
            fout.write(s[:last_parens+1])
            fout.write(';\n')
            fout.close()
        outfile = infile + '_sequences'
        paramfile = infile + '_params'

        # set up iSG simulation parameters
        options = ''
        matrix = ''
        # rate substitution parameters
        # -m (string) HKY, F84 for substitution matrix
        # -a (float) gamma rate heterogeneity OR
        # -c (list) codon position heterogeneity OR
        # -g (int) numbe categories for discrete gamma-dist rate heterogeneity (b/w 2 and 32)
        rate_type = rand.randint(0, 2)
        options = {
            0: '-a ' + str(rand.uniform(0.01, 2)),
            1: '-c ' + ','.join([str(rand.uniform(0.1,0.5)),str(rand.uniform(0.01,0.1)),str(rand.uniform(0.5,0.99))]),
            2: '-a ' + str(rand.uniform(0.01, 2)) + ' -g ' + str(rand.randint(2,32))
            }[rate_type]
        matrix = {
            0: 'HKY',
            1: 'F84'
            }[rand.randint(0,1)]
        with open(paramfile, 'a') as p:
            p.write("codon\n")
            p.write("rate type: " + options + "\n")
            p.write("matrix: " + matrix + "\n")

        # invariance parameter
        # -i (float) proportion of invariable sites
        inv_rate = '-i ' + str(rand.random())
        with open(paramfile, 'a') as p:
            p.write("invariance rate: " + inv_rate + "\n")
            p.write("---------------\n")

        ret = os.system('./indel-seq-gen -m %s %s --outfile %s < %s' % (matrix, options, outfile, infile))
        if ret != 0:
            raise OSError, 'indel-seq-gen failed with code %s ' % ret

def write_root(root_file, num_genefams):
    """Write root sequences in iSG format by pulling from given file with root sequences"""
    print "Step 1: Writing root sequences in indelSeqGen format"
    with open(root_file, 'r') as fin:
        count = 1;
        while count <= num_genefams:
            line = fin.readline()
            if not line:
                break
            line = line.split('|')
            seq = line[1].upper()
            l = len(seq) - 1
#            if l > 200: # only take amino acid sequences with at least 200
#            if l > 600: # only take sequences with at least 600
            out = root_file + str(count)
            with open(out, 'w') as fout:
                fout.write(str(l))
                fout.write('\n')
                fout.write('0'*l+'\n')
                fout.write(seq)
                fout.close()
            count = count + 1

if __name__ == "__main__":

    print "Warning: This script assumes your root sequence file corresponds to your sequence type choice, i.e. nucleotides or amino acids"

    import argparse
    parser = argparse.ArgumentParser(description='phase1 of homology simulation pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-sp', '--species_tree', help='shared species tree')
    parser.add_argument('-n', '--num_genefams', default=10, help='number of genes to generate')
    parser.add_argument('-r', '--root_seq', help='root sequence file, should have all the root sequences you will intend to use, format should be id|sequence')
    parser.add_argument('-w', '--wr_flag', default=0, help='flag for whether root sequences need to be written or not')
    parser.add_argument('-d', '--out_dir', help='out directory')

    opts = parser.parse_args()

    species_tree_file = opts.species_tree
    species_tree = dendropy.Tree.get_from_path(species_tree_file, schema="newick", rooted=True)
    num_genefams = int(opts.num_genefams)
    root_seq = opts.root_seq
    wr_flag = int(opts.wr_flag)
    out_dir = opts.out_dir

    if wr_flag == 1:
        write_root(root_seq, num_genefams)

    simulate_sequences(species_tree, root_seq, species_tree_file, num_genefams, out_dir)
