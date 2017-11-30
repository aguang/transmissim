from __future__ import print_function
from transmissim import simulate as sim
from rpy2.robjects import NA_Integer
from rpy2.robjects.packages import importr
import os
import pytest
import pyvolve
import readline
from Bio import SeqIO
from ete3 import Tree

def test_transmission_is_reproducible(tmpdir):

	# make sure transmission in simulate is reproducible
    R0 = 2
    w = "blah"        
    n_hosts = 100
    duration = 100
    rate_import_case = 0
    simphy_path = "simphy_lnx64"
    seed = 998877
    birth_rate = 0.1
    death_rate = 0.1
#    sim.transmission(R0, w, n_hosts, duration, rate_import_case, tmpdir, simphy_path, seed, birth_rate, death_rate)

    # should be 4 directories + simulated_tree.tre + simulated_viral.tre
#    assert(len(tmpdir.listdir())) == 6

#    p = tmpdir.join("simulated_tree.tre")
#    sim_tree = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_out/simulated_tree.tre')
#    tt_out = ""
#    with open(sim_tree) as f:
#        tt_out = f.read()
#    assert p.read() == tt_out

#    q = tmpdir.join("simulated_viral.tre")
#    vt_out = ""
#    sim_viral = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_out/simulated_viral.tre')
#    with open(sim_viral) as f:
#        vt_out = f.read()
#    assert q.read() == vt_out

class TestSequence:
    @pytest.fixture(autouse=True)
    def setup(self, tmpdir):
        self.tmpdir = tmpdir.strpath

    # def test_sequence(self):
    #     full_tree_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'simulated_tree_0.tre')
    #     root_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'hiv-db.fasta')
    #     sequence_out = self.tmpdir
    #     with open(full_tree_path, 'r') as f:
    #         full_tree = f.read()
    #     sim.sequence(full_tree, root_file, sequence_out)
    #     assert 1

@pytest.fixture
def yaml_config():
        import yaml
        yaml_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'params.yaml')
        with open(yaml_file, 'r') as f:
            config = yaml.load(f)
        #yield yaml.load(yaml_file)
        return(config)

@pytest.fixture
def yaml_config2():
        import yaml
        yaml_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'params2.yaml')
        with open(yaml_file, 'r') as f:
            config = yaml.load(f)
        #yield yaml.load(yaml_file)
        return(config)

def test_passes(yaml_config):
    #sim.main(yaml_config) # need to figure out how to import simphy and outbreaker
    assert 1

def test_only_transmission_tree(yaml_config):
    yaml_config['phylogeny']['viral'] = 0
    #sim.main(yaml_config2)
    assert 1

def get_individuals_from_fasta(fasta):
    individuals = set([])
    for record in SeqIO.parse(fasta, "fasta"):
        taxon = record.id.split('-')[0]
        individuals.add(taxon)
    return(individuals)

# def test_transmission_tree_sequences_relative_to_viral(yaml_config, yaml_config2, tmpdir):
#     cfg_viral = yaml_config
#     cfg_trans = yaml_config2

#     dir_viral = tmpdir.mkdir("viral")
#     print(dir_viral)
#     cfg_viral['output']['phylogenyout'] = dir_viral.strpath
#     cfg_viral['output']['sequenceout'] = dir_viral.strpath
#     sim.main(cfg_viral)

#     dir_trans = tmpdir.mkdir("trans")
#     print(dir_trans)
#     cfg_trans['output']['phylogenyout'] = dir_trans.strpath
#     cfg_trans['output']['sequenceout'] = dir_trans.strpath
#     cfg_trans['phylogeny']['viral'] = 0
#     sim.main(cfg_trans)

#     # test that dir_viral and dir_trans have the same transmission trees
#     tree_viral = Tree(dir_viral.join("simulated_tree_0.tre").strpath)
#     tree_trans = Tree(dir_trans.join("simulated_tree_0.tre").strpath)
#     assert(tree_viral.robinson_foulds(tree_trans)[0] == 0)

#     # test that dir_viral and dir_trans have the same sequence taxons
#     sequence_viral = get_individuals_from_fasta(dir_viral.join("simulated_alignment.fasta").strpath)
#     sequence_trans = get_individuals_from_fasta(dir_trans.join("simulated_alignment.fasta").strpath)
#     assert(sequence_viral == sequence_trans)

    # test that dir_viral and dir_trans have the same reads
    # ehhhhh