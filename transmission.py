#!/usr/bin/env python
from rpy2.robjects.packages import importr
from rpy2.robjects import NA_Integer
from collections import defaultdict
from itertools import izip as zip, count

outbreaker = importr('outbreaker')
base = importr('base')
w = base.rep(0.8, 350)
test = outbreaker.simOutbreak(R0=2, infec_curve=w, n_hosts=200, duration=350)

def binary_tree(ob):
	ances = ob[4]
	sources=group_ancestors(ances)
    # branch length (will deal with later)
    # all other nodes from sister branches to P0 in reverse order
    # if the node infected others then it will have its own subtree as sister branch instead

# group all nodes that share an ancestor
def group_ancestors(ances):
	sources = defaultdict(list)
	for i,s in zip(count(), ances):
		if s is not NA_Integer:
			sources[s].append(i)
	return sources

# for each group, construct pair tree
def pair_infected(sources):
	pair_infected = {}
	for group, infected in sources.items():
		group = group - 1 # -1 because python indices are -1 from R
		l=len(infected)-1
		pairs=[]
		for c in range(l):
			pairs.append((c+1, str(infected[c]))) # c+1 rep index of pair, str represents node
		pairs.append((str(infected[l]),str(group)))
		pair_infected[group] = pairs

	# last node to be infected by ancestor forms pair with ancestor
	return pair_infected

# assign branch lengths
def branch_lengths(pairs, onset):
	return tree	

binary_tree(test)
