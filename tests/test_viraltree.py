from ete3 import Tree
from rpy2.robjects import NA_Integer
from transmissim import viraltree
import pytest
import os

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
		onset = [0, 43, 128]
		ances = [NA_Integer, 1, 1]
		transmission_times = viraltree.local_transmission_time(onset, ances)
		assert transmission_times == [0, 43, 128]
		assert transmission_times == onset
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

class TestJoint:
	def test_simple(self):
		t_a = Tree("(1-1:100,(1-2:75,(1-3:25,1-4:25):50):25);")
		t_b = Tree("((2-1:25, 2-2:25):25, 2-3:50);")
		vt = viraltree.joint_viral_tree([t_a, t_b], [NA_Integer, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 1], [0, 50], 100)
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
		vt = viraltree.joint_viral_tree([t_a, t_b, t_c], [NA_Integer, s, t_a.search_nodes(name="1-2")[0]], [NA_Integer, 1, 1], [0, 25, 75], 100)
		t = Tree("(1-1:100,(((2-1:25, 2-2:25):50, 2-3:75):0.0, ((((3-1:15, 3-2:15):10, 3-3:25):0.0, 1-2:25):50, (1-3:25, 1-4:25):50):0.0):25);")
		assert vt.robinson_foulds(t)[0] == 0 # topological similarity
		# branch length similarity
		assert vt.get_common_ancestor("2-1", "2-2", "2-3").dist == t.get_common_ancestor("2-1", "2-2", "2-3").dist
		assert vt.get_common_ancestor("1-2", "2-1").dist == t.get_common_ancestor("1-2", "2-1").dist

	# Todo: write test with a->b->c and internal node infections

class TestViral:
	def test_simple(self):
		onset = [0, 43, 128]
		ances = [NA_Integer, 1, 1]
		duration = 350
		birth_rate = 0.1
		death_rate = 0.1
		simphy = "simphy_mac64"
		seed = 112233
		out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_viral")
		vt = viraltree.viral(onset, ances, duration, birth_rate, death_rate, simphy, seed, out_dir)
		print(vt)
		assert 1
