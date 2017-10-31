#!/usr/bin/env python
import readline
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from transmissim.transmission import binary_trees
from transmissim.viraltree import viral
#from transmission import binary_trees
#from viraltree import viral
import random
from ete3 import Tree
import pyvolve
from itertools import groupby
import os
from Bio import SeqIO
from collections import defaultdict
import random
import multiprocessing as mp

def sequence(full_tree, root_file):
    # gene sequences
    tree = pyvolve.read_tree(tree=full_tree, scale_tree = 0.001)
    print("tree is ok...")
    model = pyvolve.Model("nucleotide")
    root = ''
    with open(root_file, 'r') as f:
        root=f.read()
    my_partition = pyvolve.Partition(models = model, root_sequence=root)
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree=tree)
    my_evolver()

def outbreaker(cluster_R0, n_hosts, cluster_duration, rate_import_case):
    outbreaker = importr('outbreaker')
    base = importr('base')
    w = base.rep(0.8, 365)
    rseed = robjects.r['set.seed']
    rseed(seed)
    success = 0
    test = 0
    while success == 0: # make sure simulation yields a transmission network
        rseed(seed+1)
        test = outbreaker.simOutbreak(R0 = cluster_R0, infec_curve=w, n_hosts=n_hosts, duration=cluster_duration, rate_import_case=rate_import_case)
        if len(test[4]) > 1:
            success = 1
    return(test)

def transmission_tree(network, cluster_duration, ancestral_duration, out, simphy_path, seed, birth_rate, death_rate):
    ances = [i for i in network[4]]
    onset = [i for i in network[2]]
    vt = viral(onset, cluster_duration, ancestral_duration, ances, birth_rate, death_rate, seed, simphy_path, out)
    full_trees = binary_trees(network)
    i = 0
    with open('%s/clusters.txt' % (out), 'w') as g:
        for tree in full_trees:
            print(tree)
            taxa = Tree(tree).get_leaves()
            for j in taxa:
                g.write(j.name + ' ')
            g.write('\n')
            with open('%s/simulated_tree_%s.tre' % (out, i), 'w') as f:
               f.write(tree)
            i = i+1
    with open('%s/simulated_viral.tre' % (out), 'w') as f:
        f.write(vt.write(format=5))

    return vt

def reads(art, sequencing_system, reads_out, read_length, coverage):
    # genomic reads
    os.system('%s -ss %s -i simulated_alignment.fasta -o %s -l %s -f %s -m %s -s %s -p' % (art, sequencing_system, reads_out, read_length, coverage, mean_fragment_length, sd_fragment_length))
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
        print(analysis_start)
        analysis_end = options[1]
        seed = options[2]
        if not seed:
            seed = random.randint(1,4294967295)
        else: seed = int(seed)
        # transmission tree options
        cluster_R0 = float(options[7])
        cluster_death = float(options[8])
        n_hosts = int(options[9])
        cluster_duration = int(options[10])
        rate_import_case = float(options[11])

        ancestral_R0 = float(options[12])
        ancestral_death = float(options[13])
        ancestral_duration = float(options[14])

        tree_out = options[15]

        # viral tree options
        simphy_path = options[18]
        birth_rate = options[19]
        death_rate = options[20]

        # pyvolve options
        full_tree = options[21]
        root_file = options[25]

        # ART options
        art = options[28]
        reads_out = options[29]
        sequencing_system = options[30]
        read_length = options[31]
        coverage = options[32]
        mean_fragment_length = options[33]
        sd_fragment_length = options[34]

        print(analysis_start == "Pipeline: all")
        if analysis_start == "all":
            network = outbreaker(cluster_R0, n_hosts, cluster_duration, rate_import_case)
            vt = transmission_tree(network, cluster_duration, ancestral_duration, tree_out, simphy_path, seed, birth_rate, death_rate)
            if analysis_end != "TT":
                t = vt.write(format=5)
                sequence(t, root_file)
                print("sequence finished")
                if analysis_end != "GS":
                    reads(art, sequencing_system,reads_out,read_length,coverage)

        if analysis_start == "TT":
            print("Pipeline: Transmission Tree -> Reads")
            t=''
            with open(full_tree, 'r') as f:
                t=f.readline()
            sequence(t, root_file)
            if analysis_end != "GS":
                reads(art, sequencing_system, reads_out, read_length, coverage)

        if analysis_start == "GS":
            reads(art, sequencing_system, reads_out, read_length, coverage)
