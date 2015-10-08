#!/usr/bin/env python
import os
import sys, string
import random as rand
from ete2 import Tree

# simulate sequences from specified tree
def simulate_sequences(tree, root_file, out):
    """Simulate sequences from set of gene trees using iSG. Rates of evolution drawn from random distributions."""
    # write tree_file in [:root_file] "Label" (tree) format
    infile = out
    with open(infile, 'w') as fout:
        fout.write('[:' + root_file + '] ')
        fout.write('" Label ' + '" ')
        fout.write("{9, .0012, indel-length-dist.txt}")
        s=tree.write(format=5)
        fout.write(s)
        fout.write('\n')
        fout.close()
    outfile = infile + '_sequences'

    # set up iSG simulation parameters
    options = ''
    matrix = ''
    # rate substitution parameters
    # -m (string) HKY for substitution matrix
    # -a (float) gamma rate heterogeneity OR
    # -c (list) codon position heterogeneity OR
    # -g (int) number categories for discrete gamma-dist rate heterogeneity (b/w 2 and 32)
    options = '-a 0.8473'
    matrix = 'HKY'
        
    # invariance parameter
    # -i (float) proportion of invariable sites
    inv_rate = '-i .2577'

    path = os.path.dirname(sys.argv[0])
    path = os.path.abspath(path)
    ret = os.system('%s/indel-seq-gen -m %s %s --outfile %s < %s' % (path, matrix, options, outfile, infile))
#    raise OSError, 'indel-seq-gen failed with code %s ' % ret

"""
def write_root(root_file, num_genefams):
    Write root sequences in iSG format by pulling from given file with root sequences
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
            if l > 600: # only take sequences with at least 600
                out = root_file + str(count)
                with open(out, 'w') as fout:
                    fout.write(str(l))
                    fout.write('\n')
                    fout.write('0'*l+'\n')
                    fout.write(seq)
                    fout.close()
                count = count + 1
"""

if __name__ == "__main__":

#    print "Warning: This script assumes your root sequence file corresponds to your sequence type choice, i.e. nucleotides or amino acids"

    import argparse
    parser = argparse.ArgumentParser(description='phase1 of homology simulation pipeline')

    ds = ' [%(default)s]'
    parser.add_argument('-t', '--tree', help='tree for simulation')
    parser.add_argument('-r', '--root_seq', help='root sequence file, should be in a format suitable for iSG, if not rerun script with -w 1 to turn the write_root flag on')
#    parser.add_argument('-w', '--wr_flag', default=0, help='flag for whether root sequences need to be written or not')
    parser.add_argument('-d', '--out', help='out directory')

    opts = parser.parse_args()

    tree_file = opts.tree
    tree = Tree(tree_file, format=1)
    root_seq = opts.root_seq
#    wr_flag = int(opts.wr_flag)
    out = opts.out

#    if wr_flag == 1:
#        write_root(root_seq)

    simulate_sequences(tree, root_seq, out)
