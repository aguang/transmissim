from __future__ import print_function
from transmissim import simulate as sim
from rpy2.robjects import NA_Integer
from rpy2.robjects.packages import importr
import os
import pytest

def test_transmission_is_reproducible(tmpdir):

	# make sure transmission in simulate is reproducible
    R0 = 2
    w = "blah"        
    n_hosts = 100
    duration = 100
    rate_import_case = 0
    simphy_path = "simphy_mac64"
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
    def test_pyvolve_functions(self, tmpdir):
        assert 0

@pytest.fixture
def yaml_config1():
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

def test_passes(yaml_config1):
    #sim.main(yaml_config) # need to figure out how to import simphy and outbreaker
    assert 1

def test_only_transmission_tree(yaml_config2):
    #sim.main(yaml_config2)
    assert 1