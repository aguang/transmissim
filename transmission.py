#!/usr/bin/env python
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from collections import defaultdict

outbreaker = importr('outbreaker')
base = importr('base')
w = base.rep(0.8, 350)
test = outbreaker.simOutbreak(R0=2, infec_curve=w, n_hosts=200, duration=350)

def binary_tree(ob):
	ances = ob[4]
	ances = [y for y in ances if y is not robjects.NA_Integer]

	print(ances)
	sources = defaultdict(list)
	for v,k in ances:
		sources[k].append(v)
	print(sources)
    # group all nodes that share an ancestor
    # ances = list(test[4]) #test[4] = test$ances

    # last node to be infected by ancetor forms pair with ancestor
    # branch length (will deal with later)
    # all other nodes from sister branches to P0 in reverse order
    # if the node infected others then it will have its own subtree as sister branch instead

binary_tree(test)