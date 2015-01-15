#Overview

This repository aims to host a set of Python scripts that will simulate transcriptomes with known orthology/paralogy relationships, gene tree topologies, and species tree topologies. It will incorporate a summary statistics script for assessing phylogenetics program performance at the homolog clustering level and at the final tree reconstruction level.

The pipeline is broken up into 2 phases:
## Phase 1
All gene trees share the same species tree, i.e. all homologs are orthologs. Codon evolution is simulated along each gene tree with biologically relevant parameters.

## Phase 2
All gene trees are simulated from the same species tree using a birth-death process with biologically relevant birth-and-death parameters. (Alternatively, user-supplied parameters) Codon evolution is simulated along each gene tree with biologically relevant parameters.

#Pipeline Organization

##Gene Tree Simulation

Gene tree simulation is performed based off of a birth-death process from a shared species tree. The process is as follows:
1. For each tree we go down each branch and draw waiting times between events from an exponential distribution.
2. If the waiting time exceeds the branch length, no births or deaths happen along that branch and we move onto the next branch.
3. If the waiting time is less than the branch length, we draw from the uniform distribution on [0, 1] to determine whether the event is a birth or a death. We then repeat the waiting time process down the remaining length of the branch.
4. Relative birth rates are 1/3 to reflect that death rates tend to be an order of magnitude higher than birth rates. Lambda in the exponential distribution is further parametrized by a gamma distribution that reflects our prior belief in birth and death (i.e. event) rates. Reference: [Demuth 2009](http://www.ncbi.nlm.nih.gov/pubmed/19153999).

Currently we have a single shared species tree based off of the gastropoda species tree generated with agalma, working on gene tree simulation.

##Sequence Simulation

Sequence simulation for each gene tree is further done based off [indel-seq-gen](http://bioinfolab.unl.edu/~cstrope/iSG/) with root sequence input from assembled transcriptomes from the gastropod project.

##Homology Assessment

Homology clustering is performed using either fablast and mcl or blastp and mcl from the [agalma](https://bitbucket.org/caseywdunn/agalma) package. Homology assessment has yet to be determined.

##Phylogeny Assessment

Phylogeny assessment is based off of congruence to the original species tree, with metrics given on the quartets, triplets, Robinson-Foulds, and Felsenstein tree metric.

#Current Progress

Sequence Simulation, i.e. Phase 1 is online. Phylogeny Assessment is completed, but not online. Gene Tree Simulation with Sequence Simulation, i.e. Phase 2 is online.