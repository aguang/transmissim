from ete3 import Tree
from rpy2.robjects import NA_Integer
from transmissim import viraltree
import pytest
import os

# dendropy tests
#class TestGenerate:
#	def test_work_if_valid_rates(self):
#		birth_rates = [0.1, 5]
#		death_rates = [0.1, 5]
#		times = [0.1, 5]
#		tree = viraltree.generate_discrete_R0_tree(birth_rates, death_rates, times)

#	def test_break_if_invalid_rates(self):
#		birth_rates = [-2, 5]
#		death_rates = [3, 5]
#		times = [2, 5]
#		with pytest.raises(AssertionError):
#			tree = viraltree.generate_discrete_R0_tree(birth_rates, death_rates, times)

# class TestIndividualViralTree

class TestSource:
	def test_one(self):
		t = Tree("(1:10.00000000,(2:6.68089429,(3:0.48337029,4:0.48337029):6.19752400):3.31910571);")
#		t = viraltree.annotate(t)
#		viraltree.transmission_source(1, t)
		transmission_time = 9.7
		seed = 112233
		source, length = viraltree.transmission_source_branch(transmission_time, t, seed)
		assert length == 4

		transmission_time = 9
		source, length = viraltree.transmission_source_branch(transmission_time, t, seed)
		assert length == 3

		transmission_time = 6.5
		source, length = viraltree.transmission_source_branch(transmission_time, t, seed)
		assert length == 3

		transmission_time = 2
		source, length = viraltree.transmission_source_branch(transmission_time, t, seed)
		assert length == 2

		transmission_time = 0
		source, length = viraltree.transmission_source_branch(transmission_time, t, seed)
		assert length == 2

		transmission_time = -1
		with pytest.raises(ValueError):
			viraltree.transmission_source_branch(transmission_time, t, seed)

	def test_len0(self):
		t = Tree("a;")
		transmission_time=9.7
		seed = 112233
		with pytest.raises(ValueError):
			viraltree.transmission_source_branch(transmission_time, t, seed)

class TestLocalTransmission:
	def test_equals_onset(self):
		onset = [0, 1, 44, 129]
		ances = [NA_Integer, 0, 1, 1]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 1, 43, 128]
	def test_two(self):
		onset = [0, 1, 51, 101]
		ances = [NA_Integer, 0, 1, 2]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 1, 50, 50]
	def test_three(self):
		onset = [0, 2, 12, 22, 32, 42]
		ances = [NA_Integer, 0, 1, 2, 1, 2]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 2, 10, 10, 30, 30]

	def test_multiple(self):
		onset = [0, 3, 8, 15, 40, 42]
		ances = [NA_Integer, 0, 1, 0, 0, 3]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 3, 5, 15, 40, 27]

class TestJoint:
	def test_simple(self):
		t_a = Tree("(1-1:100,(1-2:75,(1-3:25,1-4:25):50):25);")
		t_b = Tree("((2-1:25, 2-2:25):25, 2-3:50);")
		vt = viraltree.joint_viral_tree([t_a, t_b], [NA_Integer, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 0], [0, 50])
		t = Tree("(1-1:100,((((2-1:25, 2-2:25):25, 2-3:50):0.0, 1-2:50):25,(1-3:25, 1-4:25):50):25);")
		assert vt.robinson_foulds(t)[0] == 0 # topological similarity
		# branch length similarity
		assert vt.get_common_ancestor("2-1", "2-2", "2-3").dist == t.get_common_ancestor("2-1", "2-2", "2-3").dist
		assert vt.get_common_ancestor("1-2", "2-1").dist == t.get_common_ancestor("1-2", "2-1").dist

	def test_branch_length_boundary(self):
		t_a = Tree("(1-1:100,(1-2:75,(1-3:25,1-4:25):50):25);")
		s = t_a.get_common_ancestor("1-2", "1-3", "1-4")
		t_b = Tree("((2-1:25, 2-2:25):25, 2-3:50);")
		t_c = Tree("((3-1:15, 3-2:15):10, 3-3:25);")
		vt = viraltree.joint_viral_tree([t_a, t_b, t_c], [NA_Integer, s, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 0, 0], [0, 25, 75])
		t = Tree("(1-1:100,(((2-1:25, 2-2:25):50, 2-3:75):0.0, ((((3-1:15, 3-2:15):10, 3-3:25):0.0, 1-2:25):50, (1-3:25, 1-4:25):50):0.0):25);")
		assert vt.robinson_foulds(t)[0] == 0 # topological similarity
		# branch length similarity
		assert vt.get_common_ancestor("2-1", "2-2", "2-3").dist == t.get_common_ancestor("2-1", "2-2", "2-3").dist
		assert vt.get_common_ancestor("1-2", "2-1").dist == t.get_common_ancestor("1-2", "2-1").dist

#	def test_multiple_clusters(self):
#		t_a = Tree("(1-1:100,(1-2:75,(1-3:25,1-4:25):50):25);")
#		s = t_a.get_common_ancestor("1-2", "1-3", "1-4")
#		t_b = Tree("((2-1:25, 2-2:25):25, 2-3:50);")
#		t_c = Tree("((3-1:15, 3-2:15):10, 3-3:25);")
#		vt = viraltree.joint_viral_tree([t_a, t_b, t_c], [NA_Integer, s, NA_Integer], [NA_Integer, 1, NA_Integer], [0, 25, 75], 100)
#		t = Tree("(1-1:100,(((2-1:25, 2-2:25):50, 2-3:75):0.0, ((((3-1:15, 3-2:15):10, 3-3:25):0.0, 1-2:25):50, (1-3:25, 1-4:25):50):0.0):25);")

	# Todo: write test with a->b->c and internal node infections

# def test_multiple(tmpdir):
# 	cluster_duration = 100
# 	ancestral_duration = 200
# 	onset = [0, 5, 12, 37, 39]
# 	for i in range(len(onset)):
# 		onset[i] = onset[i] + ancestral_duration - cluster_duration
# 	onset.insert(0,0)

# 	ances = [NA_Integer, 1, NA_Integer, NA_Integer, 3]
# 	for i in range(len(ances)):
# 		if ances[i] == NA_Integer:
# 			ances[i] = 0
# 	ances.insert(0,NA_Integer)

# 	seed = 998878
# 	sampling_times = list(map(lambda x: ancestral_duration-x, onset))

# 	assert sampling_times[0] == 200
# 	assert sampling_times[1] == 100
# 	assert sampling_times[2] == 95
# 	assert sampling_times[3] == 88
# 	assert sampling_times[4] == 63
# 	assert sampling_times[5] == 61
# #	birth_rates = [0.5,4,0.5]
# #	death_rates = [0.5,0.5,0.5]
# #	time_intervals = [5,25,float('inf')]
# 	birth_rate = 0.5
# 	death_rate = 0.5
# #	time_intervals = [float('inf')]
# #	print(tmpdir)
# #	viral_trees = viraltree.make_list_of_individual_viral_trees(sampling_times, birth_rates, death_rates, time_intervals, seed, tmpdir)
# #sampling_times, birth_rates, death_rates, seed, simphy, out_dir
# 	simphy = "simphy_mac64"
# 	viral_trees = viraltree.make_list_of_individual_viral_trees(sampling_times,birth_rate,death_rate,seed,simphy,tmpdir)
# 	assert len(viral_trees) == 6

# def test_viral(tmpdir):
# 	cluster_duration = 100
# 	ancestral_duration = 200
# 	onset = [0, 5, 12, 37, 39]
# 	ances = [NA_Integer, 1, NA_Integer, NA_Integer, 3]
# 	seed = 998878
# 	birth_rate = 0.5
# 	death_rate = 0.5
# 	simphy = "simphy_mac64"
# 	#onset, cluster_duration, ancestral_duration, ances, birth_rates, death_rates, time_intervals, seed, out_dir
# 	vt = viraltree.viral(onset, cluster_duration, ancestral_duration, ances, birth_rate, death_rate, seed, simphy, tmpdir)
# 	height = vt.get_distance(vt.get_leaves()[0])
# 	assert(height==ancestral_duration)

#class TestViral:
#	def test_simple(self):
#		onset = [0, 43, 128]
#		ances = [NA_Integer, 1, 1]
#		duration = 350
#		birth_rate = 0.1
#		death_rate = 0.1
#		seed = 112233
#		out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_viral")
#		print(out_dir)
#		vt = viraltree.viral(onset, cluster_duration, ancestral_duration, ances, birth_rate, death_rate, seed, out_dir)
#		print(vt)
#		assert 1