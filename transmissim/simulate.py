#!/usr/bin/env python
import readline
from rpy2.robjects.packages import importr
from transmission import binary_tree
from viraltree import viral
from ete3 import Tree
import pyvolve
from itertools import groupby
import os
from Bio import SeqIO
from collections import defaultdict

def sequence(full_tree, root_file):
    # gene sequences
    tree = pyvolve.read_tree(tree=full_tree, scale_tree = 0.001)
    model = pyvolve.Model("nucleotide")
    root = ''
    with open(root_file, 'r') as f:
        root=f.read()
    my_partition = pyvolve.Partition(models = model, root_sequence=root)
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree=tree)
    my_evolver()

def transmission(R0, w, n_hosts, duration, rate_import_case, out, simphy_path, seed, birth_rate, death_rate):
    outbreaker = importr('outbreaker')
    base = importr('base')
    w = base.rep(0.8, 350)
    success = 0
    test = 0
    while success == 0: # make sure simulation yields a transmission network
        test = outbreaker.simOutbreak(R0 = R0, infec_curve=w, n_hosts=n_hosts, duration=duration, rate_import_case=rate_import_case)
        if len(test[4]) > 1:
            success = 1
    ances = test[4]
    onset = test[2]

    vt = viral(onset, ances, duration, birth_rate, death_rate, simphy_path, seed, out)
    full_tree = binary_tree(test)
    with open('simulated_tree.tre', 'w') as f:
        f.write(full_tree)
    with open('simulated_viral.tre', 'w') as f:
        f.write(vt.write(format=3))

    return vt,full_tree

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
        seed = int(options[2])

        # transmission tree options
        R0 = int(options[3])
        w = options[4] # currently disabled
        n_hosts = int(options[5])
        duration = int(options[6])
        rate_import_case = int(options[7])

        # viral tree options
        simphy_path = options[9]
        birth_rate = options[10]
        print("birth rate: ", birth_rate)
        death_rate = options[11]
        print("death rate: ", death_rate)

        # pyvolve options
        full_tree = options[13]
        root_file = options[14]

        # ART options
        art = options[16]
        reads_out = options[17]
        sequencing_system = options[18]
        read_length = options[19]
        coverage = options[20]

        if analysis_start == "all" or "TN":
            vt,full_tree = transmission(R0, w, n_hosts, duration, rate_import_case, "./", simphy_path, seed, birth_rate, death_rate)
            if analysis_end != "TT":
                t = vt.write(format=5)
                sequence(t, root_file)
                if analysis_end != "GS":
                    reads(art, sequencing_system,reads_out,read_length,coverage)

        if analysis_start == "TT":
            sequence(full_tree, root_file)
            if analysis_end != "GS":
                reads(art, sequencing_system, reads_out, read_length, coverage)

        if analysis_start == "GS":
            reads(art, sequencing_system, reads_out, read_length, coverage)
