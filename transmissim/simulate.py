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
import os, sys
from Bio import SeqIO
from collections import defaultdict
import random
import numpy

def sequence(full_tree, root_file, sequence_out):
    print("Sequence simulation: Pyvolve")
    tree = pyvolve.read_tree(tree=full_tree, scale_tree=0.0001)
    model = pyvolve.Model("nucleotide")
    root = ''
    with open(root_file, 'r') as f:
        root=f.read()
    my_partition = pyvolve.Partition(models = model, root_sequence=root)
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree=tree)
    my_evolver(seqfile=os.path.join(sequence_out,"simulated_alignment.fasta"),seqfmt="fasta")

def outbreaker(cluster_R0, n_hosts, cluster_duration, rate_import_case, seed, w):
    outbreaker = importr('outbreaker')
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

def transmission_tree(network, phylogeny_out):
    print("Transmission phylogeny simulation")
    full_trees = binary_trees(network)
    i = 0
    with open('%sclusters.txt' % (phylogeny_out), 'w') as g:
        for tree in full_trees:
            taxa = Tree(tree).get_leaves()
            for j in taxa:
                g.write(j.name + ' ')
            g.write('\n')
            with open('%s/simulated_tree_%s.tre' % (phylogeny_out, i), 'w') as f:
               f.write(tree)
            i = i+1
    # Trees have branch lengths in # of mutations TOTAL over time (not /site)
    return(full_trees)

def viral_tree(network, duration, birth_rate, death_rate, seed, simphy_path, out):
    print("Viral phylogeny simulation")
    ances = [i for i in network[4]]
    onset = [i for i in network[2]]
    vt = viral(onset, duration, ances, birth_rate, death_rate, seed, simphy_path, out)
    with open('%s/simulated_viral.tre' % (out), 'w') as f:
        f.write(vt.write(format=5))
    # Trees have branch lengths in # of generations
    return vt.write(format=5)

def reads(art, sequencing_system, sequence_out, read_length, coverage):
    reads_out = os.path.join(sequence_out, "reference")
    sim_align = os.path.join(sequence_out, "simulated_alignment.fasta")
    # genomic reads
    #os.system('%s -ss %s -i simulated_alignment.fasta -o %s -l %s -f %s -m %s -s %s' % (art, sequencing_system, reads_out, read_length, coverage, mean_fragment_length, sd_fragment_length))
    os.system('%s -ss %s -i %s -o %s -l %s -f %s' % (art, sequencing_system, sim_align, reads_out, read_length, coverage))
    # split by taxon
    record_dict = defaultdict(list)
    for record in SeqIO.parse(os.path.join(sequence_out, "reference.fq"), "fastq"):
        taxon = record.id.split('-')[0]
        record_dict[taxon].append(record)

    for k in record_dict.keys():
        with open(os.path.join(sequence_out, k+".fa"), 'w') as f:
            for record in record_dict[k]:
                SeqIO.write(record, f, "fasta")
                f.write(">%s'\n" % record.id)
                f.write(str(record.seq.reverse_complement())+"\n")

def network(sim_contact, R0, number_of_hosts, duration, rate_import_case, seed, infection_curve):
    network = 0
    if(sim_contact != 0):
        print("Contact network simulation not implemented yet.")
        sys.exit()
    else:
        base = importr('base')
        v=[value for i in infection_curve for value in numpy.repeat(float(i[0]),i[1]).tolist()]
        w = robjects.FloatVector(v)
        network = outbreaker(cluster_R0 = R0, n_hosts = number_of_hosts,
                cluster_duration = duration, rate_import_case = rate_import_case,
                seed=seed, w=w)
    return(network)

def net_phylo_seq(sim_contact, R0, number_of_hosts, duration, rate_import_case,
    seed, sim_transmission_tree, sim_viral, root_sequence, phylogeny_out, birth_rate,
    death_rate, sequencing_system, read_length, coverage, viral_tree_program,
    reads_program, sim_reads, sequence_out, infection_curve):

    nk = network(sim_contact, R0, number_of_hosts, duration, rate_import_case, seed, infection_curve)

    tt = False
    vt = False
    if(sim_transmission_tree):
        tt = transmission_tree(nk, phylogeny_out)
    if(sim_viral):
        vt = viral_tree(nk, duration, birth_rate, death_rate, seed, viral_tree_program, phylogeny_out)

    if(vt):
        sequence(vt, root_sequence, sequence_out)
    elif(tt):
        if(len(tt) > 1):
            print("Multiple transmission trees not implemented yet.")
            sys.exit()
        else:
            sequence(tt[0], root_sequence, sequence_out)

    if(sim_reads):
        reads(reads_program, sequencing_system, sequence_out, read_length, coverage)

def main(cfg):
    seed = random.randint(1,4294967295)
    if cfg['main']['seed']:
        print('Seed: ', cfg['main']['seed'])
        seed = int(cfg['main']['seed'])

    sim_network = cfg['modules']['network']
    sim_phylogeny = cfg['modules']['phylogeny']
    sim_sequence = cfg['modules']['sequence']

    # assert sim_network, sim_phylogeny, sim_sequence are booleans

    viral_tree_program = cfg['programs']['viraltreeprogram']
    reads_program = cfg['programs']['readsprogram']

    phylogeny_out = cfg['output']['phylogenyout']
    sequence_out = cfg['output']['sequenceout']

    sim_contact = cfg['network']['contact']
    sim_transmission_network = cfg['network']['transmission']
    R0 = cfg['network']['R0']
    number_of_hosts = cfg['network']['numberofhosts']
    duration = cfg['network']['duration']
    rate_import_case = cfg['network']['rateimportcase']
    infection_curve = cfg['network']['infection_curve']

    sim_transmission_tree = cfg['phylogeny']['transmission']
    sim_viral = cfg['phylogeny']['viral']
    birth_rate = cfg['phylogeny']['birthrate']
    death_rate = cfg['phylogeny']['birthrate']

    sim_reads = cfg['sequence']['reads']
    root_sequence = cfg['sequence']['rootsequence']
    sequencing_system = cfg['sequence']['sequencingsystem']
    read_length = cfg['sequence']['readlength']
    coverage = cfg['sequence']['coverage']

    if(sim_network and sim_phylogeny and sim_sequence):
        print("Module option: Simulate network, phylogeny, and sequences.")
        net_phylo_seq(sim_contact, R0, number_of_hosts, duration, rate_import_case,
            seed, sim_transmission_tree, sim_viral, root_sequence, phylogeny_out,
            birth_rate, death_rate, sequencing_system, read_length, coverage,
            viral_tree_program, reads_program, sim_reads, sequence_out, infection_curve)

    elif(sim_network and sim_phylogeny and not sim_sequence):
        # run just network & phylogeny
        print("Not implemented yet.")

    elif(sim_network and not sim_phylogeny and sim_sequence):
        print("Network to Sequence simulation not implemented yet.")

    elif(not sim_network and sim_phylogeny and sim_sequence):
        # just simulate transmission tree & sequences
        print("Not implemented yet.")

    elif(not sim_network and not sim_phylogeny and sim_sequence):
        # just simulate genome sequences
        print("Not implemented yet.")

    elif(not sim_network and sim_phylogeny and not sim_sequence):
        # just simulate phylogeny
        print("Not implemented yet.")

    elif(sim_network and not sim_phylogeny and not sim_sequence):
        # just simulate network
        nk = network(sim_contact, R0, number_of_hosts, duration, rate_import_case, seed)
        tree = binary_trees(nk)
        print("Tree: ", tree)

    else:
        print("All modules 0, nothing to simulate.")
