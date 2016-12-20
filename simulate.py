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

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params
    if params:
        options = [line.strip().split('//')[0].strip() for line in open(params)]
        
        analysis = options[0]
        R0 = int(options[1])
        w = options[2] # currently disabled
        n_hosts = int(options[3])
        duration = int(options[4])
        rate_import_case = int(options[5])

        # pyvolve options
        root_file = options[7]

        # ART options
        art = options[9]
        reads_out = options[10]
        sequencing_system = options[11]
        read_length = options[12]
        coverage = options[13]

        if analysis == "all" or "TN->R":

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

            # gene sequences
            tree = pyvolve.read_tree(tree=full_tree)
            pyvolve.print_tree(tree)
            model = pyvolve.Model("nucleotide")
            root = ''
            with open(root_file, 'r') as f:
                root=f.read()
            my_partition = pyvolve.Partition(models = model, root_sequence=root)
            my_evolver = pyvolve.Evolver(partitions = my_partition, tree=tree)
            my_evolver()

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
