#!/usr/bin/env python
from rpy2.robjects.packages import importr
from transmission import binary_tree
from ete3 import Tree
import genesequence
import pyvolve
from itertools import groupby
import os
from Bio import SeqIO
from collections import defaultdict

def sequence(full_tree, root_file):
    # gene sequences
    tree = pyvolve.read_tree(tree=full_tree, scale_tree = 0.0001)
    model = pyvolve.Model("nucleotide")
    root = ''
    with open(root_file, 'r') as f:
        root=f.read()
    my_partition = pyvolve.Partition(models = model, root_sequence=root)
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree=tree)
    my_evolver()

def transmission(R0, w, n_hosts, duration, rate_import_case):
    outbreaker = importr('outbreaker')
    base = importr('base')
    w = base.rep(0.8, 350)
    success = 0
    test = 0
    while success == 0: # make sure simulation yields a transmission network
        test = outbreaker.simOutbreak(R0 = R0, infec_curve=w, n_hosts=n_hosts, duration=duration, rate_import_case=rate_import_case)
        if len(test[4]) > 1:
            success = 1
            
    full_tree = binary_tree(test)
    with open('simulated_tree.tre', 'w') as f:
        f.write(full_tree)

    return full_tree

def reads(art, sequencing_system, reads_out, read_length, coverage):
    # genomic reads
    os.system('%s -ss %s -i simulated_alignment.fasta -o %s -l %s -f %s' % (art, sequencing_system, reads_out, read_length, coverage))
    # split by taxon
    record_dict = defaultdict(list)
    for record in SeqIO.parse(reads_out+".fq", "fastq"):
        taxon = record.id.split('-')[0]
        record_dict[taxon].append(record)

    for k in record_dict.keys():
        with open(reads_out+k+".fa", 'w') as f:
            for record in record_dict[k]:
                SeqIO.write(record, f, "fasta")
                f.write(">%s'\n" % record.id)
                f.write(str(record.seq.reverse_complement())+"\n")

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params
    if params:
        options = [line.strip().split('//')[0].strip() for line in open(params)]
        
        analysis_start = options[0]
        analysis_end = options[1]

        R0 = int(options[2])
        w = options[3] # currently disabled
        n_hosts = int(options[4])
        duration = int(options[5])
        rate_import_case = int(options[6])

        # pyvolve options
        full_tree = options[8]
        root_file = options[9]

        # ART options
        art = options[11]
        reads_out = options[12]
        sequencing_system = options[13]
        read_length = options[14]
        coverage = options[15]

        if analysis_start == "all" or "TN":
            full_tree = transmission(R0, w, n_hosts, duration, rate_import_case)
            if analysis_end != "TT":
                sequence(full_tree, root_file)
                if analysis_end != "GS":
                    reads(art, sequencing_system,reads_out,read_length,coverage)

        if analysis_start == "TT":
            sequence(full_tree, root_file)
            if analysis_end != "GS":
                reads(art, sequencing_system, reads_out, read_length, coverage)

        if analysis_start == "GS":
            reads(art, sequencing_system, reads_out, read_length, coverage)