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

class TestTransmissionSourceBranch:
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

	def test_source_branch_is_in_choices(self):
		onset = [0.0, 189.0]
		ances = [NA_Integer, 1]
		seed = 6
		simphy = "simphy_mac64"

		local_transmission_times = viraltree.local_transmission_time(onset, ances)
		assert len(local_transmission_times) == 2
		t1 = Tree('(1:200.00000000,(2:171.66866975,3:171.66866975):28.33133025);')
		t2 = Tree('(1:11.00000000,2:11.00000000);')

		transmission_time1 = local_transmission_times[0]
		assert(transmission_time1 == NA_Integer)

		transmission_time2 = local_transmission_times[1]
		source, len_sources = viraltree.transmission_source_branch(transmission_time2, t1, seed)
		assert(len_sources == 3)
		assert(int(source.name) in list(range(1,len_sources+1)))

class TestLocalTransmissionTime:
	def test_equals_onset(self):
		onset = [0, 43, 128]
		ances = [NA_Integer, 1, 1]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 43, 128]
	def test_two(self):
		onset = [0, 50, 100]
		ances = [NA_Integer, 1, 2]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 50, 50]
	def test_three(self):
		onset = [0, 10, 20, 30, 40]
		ances = [NA_Integer, 1, 2, 1, 2]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 10, 10, 30, 30]

class TestJointViralTree:
	def test_simple(self):
		t_a = Tree("(1-1:100,(1-2:75,(1-3:25,1-4:25):50):25);")
		t_b = Tree("((2-1:25, 2-2:25):25, 2-3:50);")
		vt = viraltree.joint_viral_tree([t_a, t_b], [NA_Integer, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 1], [0, 50])
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
		vt = viraltree.joint_viral_tree([t_a, t_b, t_c], [NA_Integer, s, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 1, 1], [0, 25, 75])
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

class TestMakeListOfIndividualViralTrees:
	@pytest.fixture(autouse=True)
	def setup(self, tmpdir):
		self.tmpdir = tmpdir.strpath

	def test_three(self):
		seed = 998878
		sampling_times = [100, 50, 10]
		birth_rate = 0.5
		death_rate = 0.5
		simphy = "simphy_lnx64"
		out_dir = self.tmpdir
		viral_trees = viraltree.make_list_of_individual_viral_trees(sampling_times, birth_rate, death_rate, seed, simphy, out_dir)
		assert len(viral_trees) == 3

def get_individuals_from_tree(tree):
    individuals = set([])
    for leaf in tree.iter_leaves():
    	individual = leaf.name.split('-')[0]
    	individuals.add(individual)
    return(individuals)

class TestViral:
	@pytest.fixture(autouse=True)
	def setup(self, tmpdir):
		self.tmpdir = tmpdir.strpath

	def test_num_in_joint_equals_num_individuals(self):
		onset = [0.0, 189.0]
		ances = [NA_Integer, 1]
		duration = 200
		birth_rate = 0.001
		death_rate = 0.001
		seed = 6
		simphy = "simphy_lnx64"
		out_dir = self.tmpdir

		vt = viraltree.viral(onset, duration, ances, birth_rate, death_rate, seed, simphy, out_dir)
		individuals = get_individuals_from_tree(vt)
		assert(len(individuals)==2)