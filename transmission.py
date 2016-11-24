#!/usr/bin/env python
from rpy2.robjects.packages import importr
from rpy2.robjects import NA_Integer
from collections import defaultdict
from itertools import izip as zip, count, repeat

outbreaker = importr('outbreaker')
base = importr('base')
w = base.rep(0.8, 350)
test = outbreaker.simOutbreak(R0=2, infec_curve=w, n_hosts=200, duration=350)

def binary_tree(ob):
	ances = ob[4]
	sources=group_ancestors(ances)
	pair_infected = pair_infected(sources)

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
		pairs.append((str(group), str(infected[l])))
		pair_infected[group] = pairs

	# last node to be infected by ancestor forms pair with ancestor
	return pair_infected

# assign branch lengths
def assign_bl(pair_infected, nmut, duration):
	# maintained list of branch lengths at infect time
	node_t = list(repeat(duration, len(nmut)))
	bl = {}
	for k,v in pair_infected.items():
		t0 = node_t[k] # initial infection time for computing tips
		t = t0 # initial infection time for internal bl
		bl_pairs = []
		for pair in v:
			p1 = int(pair[1])
			tips = t0 - nmut[p1]
			internal = t - tips
			t = tips # update t to tips for subsequent internals
			node_t[p1] = tips
			bl_pairs.append((tips, internal))
		bl[k] = bl_pairs
	return bl
