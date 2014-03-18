#Overview

The purpose of this repository is hosting a project on a pipeline that will perform homology testing on gene trees simulated from species tree, and assess their performance both at the homolog clustering level and at the final tree reconstruction level.

#Pipeline Organization

##Gene Tree Simulation

Gene tree simulation is performed based off of a birth-death model from a shared species tree. Currently we have a single shared species tree based off of the gastropoda species tree generated with agalma, and our birth-death process tree simulation is done using DendroPy.

##Sequence Simulation

Sequence simulation for each gene tree is further done based off of either a codon model or an amino acid model using INDELible.

##Homology Assessment

##Phylogeny Assessment

Phylogeny assessment is based off of congruence to the original species tree, with metrics given on both the quartet method and the nearest neighbor joining method.

#Current Progress