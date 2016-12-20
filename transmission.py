#!/usr/bin/env python
from rpy2.robjects.packages import importr
from rpy2.robjects import NA_Integer
from collections import defaultdict
from itertools import izip as zip, count, repeat
from re import sub

def binary_tree(ob):
	ances = ob[4] # 4=ances, 7=nmut, 2=onset
	sources=group_ancestors(ances)
	pi = pair_infected(sources)
	ft = full_tree(pi, ob[7])
	return ft # do we want Newick string or dendropy tree?

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

# full tree
def full_tree(pi, nmut):
	newick_str = ',0:1'
	for k in pi.keys():
		k_str = '{p0}'
		for pair in pi[k]:
			p1 = pair[1] # should be character
			p_str = '({p0},' + p1 + ':1):' + str(nmut[int(p1)])
			# then replace '_' in newick_str with p_str
			k_str = k_str.format(p0 = p_str)
		k_str = k_str.format(p0 = str(k)+':1')
		newick_str = newick_str.replace(',%s:1' % k, ',%s' % k_str)
	newick_str += ';'
	newick_str = newick_str[1:]
	return newick_str
