#Overview

The purpose of this repository is hosting a project on a pipeline that will perform homology testing on gene trees simulated from species tree, and assess their performance both at the homolog clustering level and at the final tree reconstruction level.

#Pipeline Organization

##Gene Tree Simulation

Gene tree simulation is performed based off of a birth-death model from a shared species tree. Currently we have a single shared species tree based off of the gastropoda species tree generated with agalma, working on gene tree simulation.

##Sequence Simulation

Sequence simulation for each gene tree is further done based off indel-seq-gen with root sequence input from assembled transcriptomes from the gastropod project.

##Homology Assessment

Homology clustering is performed using either fablast and mcl or blastp and mcl from the [agalma](https://bitbucket.org/caseywdunn/agalma) package. Assessment is done using either homologize-stats.py on the mcl file produced or cluster_analysis2.py (not uploaded to repo yet) on the blast hits.tab file. Output is respectively a histogram of the number of species for a given gene in each cluster, aggregate statistics of the percentage of clusters with all species for only one gene, the perccentage of clusters with more than 1 gene in them, and the percentage of clusters with only 1 gene in them but that did not cluster all species with that gene together, and graphs of each cluster colored by gene, or just graphs of each cluster colored by gene.

##Phylogeny Assessment

Phylogeny assessment is based off of congruence to the original species tree, with metrics given on both the quartet method and the nearest neighbor joining method. This has not actually been done yet.

#Current Progress