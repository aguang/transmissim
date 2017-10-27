#!/usr/bin/env python
from rpy2.robjects.packages import importr
from rpy2.robjects import NA_Integer
from collections import defaultdict
from itertools import count, repeat
from re import sub
from dendropy.simulate import treesim
import dendropy
import random

def binary_trees(ob):
	ances = ob[4] # 4=ances, 7=nmut, 2=onset
	sources, NA_set=group_ancestors(ances)
	pi = pair_infected(sources)
	ft = full_trees(pi, ob[2], NA_set)
	return ft # do we want Newick string or dendropy tree?

# group all nodes that share an ancestor
def group_ancestors(ances):
	sources = defaultdict(list)
	NA_set = set([])
	for i,s in zip(count(), ances):
		if s is not NA_Integer:
			sources[s].append(i)
		else:
			NA_set.add(i)
			#sources[0].append(i)

	return sources, NA_set

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

# full trees
def full_trees(pi, onset, NA_set):
	trees = []
	nodes_in_trees = defaultdict(int)
	i = 0
	for cluster_origin in NA_set:
		trees.append(',{p0}:1'.format(p0=cluster_origin))
		nodes_in_trees[cluster_origin] = i
		i = i+1

	for infection_source in pi.keys():
		newick_str = trees[nodes_in_trees[infection_source]]
		k_str = '{p0}'
		for pair in pi[infection_source]:
			p1 = pair[1] # should be character
			nodes_in_trees[int(p1)] = nodes_in_trees[infection_source]
			#print(nodes_in_trees)
			p_str = '({p0},' + p1 + ':1):' + str(onset[int(p1)])
			# then replace '_' in newick_str with p_str
			k_str = k_str.format(p0 = p_str)
		k_str = k_str.format(p0 = str(infection_source)+':1')
		newick_str = newick_str.replace(',%s:1' % infection_source, ',%s' % k_str)
		trees[nodes_in_trees[infection_source]] = newick_str

	for i in range(len(trees)):
		#print(trees[i])
		trees[i] += ';'
		trees[i] = trees[i][1:]

	return trees

# merge transmission trees together with longer timeline
def merge_trees(full_trees, ancestral_R0, ancestral_death, ancestral_time, seed):
	random.seed(seed)
	birth_rate = ancestral_R0 * ancestral_death
	ntax = len(full_trees)
	tree = treesim.birth_death_tree(birth_rate=birth_rate, death_rate=ancestral_death,max_time=ancestral_time,num_extant_tips=ntax,rng=random)	

	i = 0
	for leaf in tree.leaf_nodes():
		t = dendropy.Tree.get(data=full_trees[i],schema='newick')
		leaf.add_child(t.seed_node)
		i = i+1
	return(tree)