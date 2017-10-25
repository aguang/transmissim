import os
from numpy import random
from ete3 import Tree
import dendropy
from dendropy.simulate import treesim
import subprocess

# modified from dendropy
# potential bug in dendropy for trees with different rates @ different times
def generate_discrete_R0_tree(birth_rates, death_rates, time_intervals):
	assert(all(x >= 0 for x in birth_rates))
	assert(all(x >= 0 for x in death_rates))
	assert(all(x >= 0 for x in time_intervals))
	tree = treesim.birth_death_tree(birth_rates[0],
		death_rates[0],
		max_time=time_intervals[0],
		repeat_until_success=True)
#	print(tree.as_string(schema='newick'))
#	for i in range(1,last_index):
#		tree = treesim.birth_death_tree(birth_rates[i],
#				death_rates[i],
#				max_time=time_intervals[i],
#				tree = tree,
#				repeat_until_success=True)
#		print(tree.as_string(schema='newick'))
	return tree

# input: sampling time in specific tree, (tree parameters, seed)
# output: individual's viral phylogeny
# simulate viral tree for an individual
def individual_viral_tree(sampling_time, birth_rates, death_rates, time_intervals, seed):
	assert len(birth_rates) == len(time_intervals)
	assert len(birth_rates) == len(death_rates)
	last_index = next(x[0] for x in enumerate(time_intervals) if x[1] > sampling_time)
	times = time_intervals
	times[last_index] = sampling_time
	times[last_index+1:] = 0
	viral_tree = generate_discrete_R0_tree(birth_rates, death_rates, times)
	return(viral_tree)
#	subprocess.run([simphy, "-st f:%s" % sampling_time, "-sb f:%s" % birth_rate, "-sd f:%s" % death_rate, "-cs %s" % seed, "-sp f:10000", "-o %s" % out, "-v 0"], stdout=subprocess.DEVNULL)
#	os.system("%s -st f:%s -sb f:%s -sd f:%s -cs %s -sp f:10000 -o %s -v 0 " % (simphy, sampling_time, birth_rate, death_rate, seed, out))

# input: source time in specific tree, source tree, (seed)
# output: source branch index in specific tree
# pick sequence for transmission
def transmission_source_branch(transmission_time, source_tree, seed):
	height = source_tree.get_distance(source_tree.get_leaves()[0])
	#print("height of tree is %s: " % (height))
	if transmission_time < 0:
		raise ValueError("Transmission time cannot be < 0")
	possibilities = []
	for leaf in source_tree.iter_leaves(is_leaf_fn=lambda n: n.get_distance(source_tree) > transmission_time):
		possibilities.append(leaf)
	random.seed(seed)
	print("Length of possibilities is: %s" % (len(possibilities)))
	i=random.randint(0, len(possibilities))
	return possibilities[i], len(possibilities)

# input: source tree indices, global transmission times
# return transmission times specific to the source tree
# note: this is for outbreaker, when branch lengths would be in # of mutations as well...
def local_transmission_time(onset, ancestors):
	l = [0]
	for t,a in zip(onset[1:], ancestors[1:]):
		l.append(t-onset[a-1])
	return l

# joint together viral tree
def joint_viral_tree(individual_trees, source_branches, ances, local_time, duration):
	# number of individual trees = number of transmission sources
	assert (len(individual_trees) == len(source_branches)), "Number of individual trees should equal number of source branches"

	for i in range(1, len(individual_trees)):
		source = ances[i]-1 # ances contains the indices of the source tree
		S = individual_trees[source] 
		# find source node in source tree
		B = source_branches[i]
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

def make_list_of_individual_viral_trees(sampling_times, birth_rates, death_rates, time_intervals, seed, out_dir):
	individual_trees = []
	c = 0
	for i in sampling_times:
		viral_tree = individual_viral_tree(i, birth_rates, death_rates, time_intervals, seed)
		for leaf in viral_tree.leaf_node_iter():
			leaf.taxon = "%s-%s" % (c, leaf.taxon)
		individual_trees.append(viral_tree)
		c=c+1
	return(individual_trees)

# run all functions together
def viral(onset, duration, ances, birth_rates, death_rates, seed, out_dir):
	# todo: make more efficient by collapsing redundant for loops
	transmission_times = local_transmission_time(onset, ances)
	sampling_times = list(map(lambda x: duration-x, onset))

	individual_trees = make_list_of_individual_viral_trees(sampling_times, birth_rates, death_rates, simphy, seed, out_dir)

	source_branches = []
	for i in range(len(local_time)):
		time = local_time[i]
		print("-----")
		print("transmission time is: %s " % (time))
		print("sampling time for tree is: %s " % sampling_times[i])
		print("index for individual tree is: %s " % (ances[i]-1))
		tree = individual_trees[ances[i]-1]
		source_branches.append(transmission_source_branch(time, tree, seed)[0])

	vt = joint_viral_tree(individual_trees, source_branches, ances, local_time, duration)
	vt.dist = 0
	print("viral tree is done")
	return vt
