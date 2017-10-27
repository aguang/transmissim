import os
from numpy import random
from ete3 import Tree
from rpy2.robjects import NA_Integer
import dendropy
from dendropy.simulate import treesim
import subprocess

# currently not working out
# modified from dendropy
# potential bug in dendropy for trees with different rates @ different times
#def generate_discrete_R0_tree(birth_rates, death_rates, time_intervals):
#	assert(all(x >= 0 for x in birth_rates))
#	assert(all(x >= 0 for x in death_rates))
#	assert(all(x >= 0 for x in time_intervals))
#	tree = treesim.birth_death_tree(birth_rates[0],
#		death_rates[0],
#		max_time=time_intervals[0],
#		repeat_until_success=True)
#	print(tree.as_string(schema='newick'))
#	for i in range(1,last_index):
#		tree = treesim.birth_death_tree(birth_rates[i],
#				death_rates[i],
#				max_time=time_intervals[i],
#				tree = tree,
#				repeat_until_success=True)
#		print(tree.as_string(schema='newick'))
#	return tree

# input: sampling time in specific tree, (tree parameters, seed)
# output: individual's viral phylogeny
# simulate viral tree for an individual
#def individual_viral_tree(sampling_time, birth_rates, death_rates, time_intervals, seed):
#	assert len(birth_rates) == len(time_intervals)
#	assert len(birth_rates) == len(death_rates)
#	if len(birth_rates) > 1:
#		last_index = next(x[0] for x in enumerate(time_intervals) if x[1] > sampling_time)
#		time_intervals[last_index] = sampling_time
#		time_intervals[last_index+1:] = [0] * (len(times)-last_index-1)
#	else:
#		time_intervals[0] = sampling_time
#	print("time_interval:", time_intervals)
#	viral_tree = generate_discrete_R0_tree(birth_rates, death_rates, time_intervals)
#	print(viral_tree.as_string(schema="newick"))
#	return(viral_tree)
#	subprocess.run([simphy, "-st f:%s" % sampling_time, "-sb f:%s" % birth_rate, "-sd f:%s" % death_rate, "-cs %s" % seed, "-sp f:10000", "-o %s" % out, "-v 0"], stdout=subprocess.DEVNULL)
#	os.system("%s -st f:%s -sb f:%s -sd f:%s -cs %s -sp f:10000 -o %s -v 0 " % (simphy, sampling_time, birth_rate, death_rate, seed, out))
def individual_viral_tree(sampling_time, birth_rate, death_rate, simphy, seed, out):
#	print("out: %s" % (out))
	os.system("%s -st f:%s -sb f:%s -sd f:%s -cs %s -sp f:10000 -o %s -v 0" % (simphy, sampling_time, birth_rate, death_rate, seed, out))

# input: source time in specific tree, source tree, (seed)
# output: source branch index in specific tree
# pick sequence for transmission
def transmission_source_branch(transmission_time, source_tree, seed):
	height = source_tree.get_distance(source_tree.get_leaves()[0])
#	print("height of tree is %s: " % (height))
	if transmission_time < 0:
		raise ValueError("Transmission time cannot be < 0")
	possibilities = []
	for leaf in source_tree.iter_leaves(is_leaf_fn=lambda n: n.get_distance(source_tree) > transmission_time):
		possibilities.append(leaf)
	random.seed(seed)
#	print("Length of possibilities is: %s" % (len(possibilities)))
	i=random.randint(0, len(possibilities))
	return possibilities[i], len(possibilities)

# input: source tree indices, global transmission times
# return transmission times specific to the source tree
# note: this is for outbreaker, when branch lengths would be in # of mutations as well...
def local_transmission_time(onset, ancestors):
	l = [onset[0]]
	for t,a in zip(onset[1:], ancestors[1:]):
		l.append(t-onset[a])
	return l

# joint together viral tree
def joint_viral_tree(individual_trees, source_branches, ances, local_time):
	# number of individual trees = number of transmission sources
	assert (len(individual_trees) == len(source_branches)), "Number of individual trees should equal number of source branches"

	for i in range(1, len(individual_trees)):
		source = ances[i] # ances contains the indices of the source tree
#		print("source: ", source)
		S = individual_trees[source] 
		# find source node in source tree
		B = source_branches[i]
#		print("source_branches[i]:", B)
		# compute new distances for transmission node
		new_node_dist = local_time[i] - (B.up).get_distance(S)
		branch_dist = B.dist - new_node_dist
		# add transmission node
		T = B.up.add_child(name="T_%s_%s" % (source, i), dist=new_node_dist)
		T.add_features(host=source)
		# detach source node
		B.detach()
		# add source node as a child of transmission node
		T.add_child(B)
		B.dist = branch_dist
		# add infected tree to transmission node
		I = individual_trees[i]
		T.add_child(I)
		I.dist = 0
	vt = individual_trees[0]

	return vt

def make_list_of_individual_viral_trees(sampling_times, birth_rate, death_rate, seed, simphy, out_dir):
	individual_trees = []
	c = 0
	for i in sampling_times:
		if i == 1.0: # deal with -st error
			i = 1.000001
		out = "%s/%s" % (out_dir, c)
		#sampling_time, birth_rate, death_rate, simphy, seed, out
		individual_viral_tree(i, birth_rate, death_rate, simphy, seed, out)
		#viral_tree = individual_viral_tree(i, birth_rates, death_rates, time_intervals, seed)
		#ete_tree = Tree(viral_tree.as_string(schema="newick")[5:])
		t = Tree("%s/1/s_tree.trees" % (out))
		#for leaf in viral_tree.leaf_node_iter():
		#	print(leaf.taxon)
		#	leaf.taxon = "%s-%s" % (c, leaf.taxon)
		#individual_trees.append(viral_tree)
		for leaf in t.iter_leaves():
			leaf.name = "%s-%s" % (c, leaf.name)
		individual_trees.append(t)
		c=c+1
	return(individual_trees)

# run all functions together
def viral(onset, cluster_duration, ancestral_duration, ances, birth_rates, death_rates, seed, simphy, out_dir):
	# todo: make more efficient by collapsing redundant for loops
	for i in range(len(onset)):
		if ances[i] == NA_Integer:
			ances[i] = 0
		else:
			onset[i] = onset[i] + ancestral_duration - cluster_duration
	onset.insert(0,0)
	ances.insert(0,NA_Integer)

	#print(ances)
	#print(onset)
	local_transmission_times = local_transmission_time(onset, ances)
	#print(local_transmission_times)
	# y = x + a - c
	# originally: c - x
	# now: c - y = c - x - a + c
	# => c - y + a - c = a - y = c - x
	sampling_times = list(map(lambda x: ancestral_duration - x, onset))
	assert(local_transmission_times[i+1] <= sampling_times[i] for i in range(len(sampling_times)))
#	print(sampling_times)

	individual_trees = make_list_of_individual_viral_trees(sampling_times, birth_rates, death_rates, seed, simphy, out_dir)
#	print(len(individual_trees))

	source_branches = []
	for i in range(len(local_transmission_times)):
		time = local_transmission_times[i]
		print("-----")
		print("transmission time is: %s " % (time))
		#print("sampling time for tree is: %s " % sampling_times[i])
		print("index for individual tree is: %s " % (ances[i]))
#		print(ances[i])
		tree = individual_trees[ances[i]]
		# todo: convert ete3 to dendropy
#		print(tree.write())
		source_branches.append(transmission_source_branch(time, tree, seed)[0])

	vt = joint_viral_tree(individual_trees, source_branches, ances, local_transmission_times)
	vt.dist = 0
	print("viral tree is done")
	return vt
