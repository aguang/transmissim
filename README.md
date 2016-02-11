#Overview

This repository aims to host a set of Python scripts that will simulate transcriptomes with known orthology/paralogy relationships, gene tree topologies, and species tree topologies. It will incorporate a summary statistics script for assessing phylogenetics program performance at the homolog clustering level and at the final tree reconstruction level.

#Testing

To test, simply run:

	python simulate.py -p params.txt

#Pipeline Organization

##Species and Gene Tree Simulation

Species and gene tree simulation is done based off [SimPhy](https://github.com/adamallo/SimPhy).

##Sequence Simulation

Sequence simulation for each gene tree is further done based off [indel-seq-gen](http://bioinfolab.unl.edu/~cstrope/iSG/) with root sequence input from assembled transcriptomes from the gastropod project.

##Read Simulation

Read simulation for each gene sequence, sorted by taxa is done using [RNASeqReadSimulator](https://github.com/davidliwei/RNASeqReadSimulator).

##Homology Assessment

Homology clustering is performed using either fablast and mcl or blastp and mcl from the [agalma](https://bitbucket.org/caseywdunn/agalma) package. Homology assessment has yet to be determined.