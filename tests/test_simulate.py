from __future__ import print_function
from transmissim import simulate as sim
from rpy2.robjects import NA_Integer
import os
import pytest

def test_transmission_is_reproducible(tmpdir):
	# make sure transmission in simulate is reproducible
    R0 = 2
    w = "blah"        
    n_hosts = 100
    duration = 100
    rate_import_case = 0
    simphy_path = "simphy"
    seed = 998877
    birth_rate = 0.1
    death_rate = 0.1
    sim.transmission(R0, w, n_hosts, duration, rate_import_case, tmpdir, simphy_path, seed, birth_rate, death_rate)

    # should be 4 directories + simulated_tree.tre + simulated_viral.tre
    assert(len(tmpdir.listdir())) == 6

    p = tmpdir.join("simulated_tree.tre")
    tt_out = ""
    with open("test_out/simulated_tree.tre") as f:
        tt_out = f.read()
    assert p.read() == tt_out

    q = tmpdir.join("simulated_viral.tre")
    vt_out = ""
    with open("test_out/simulated_viral.tre") as f:
        vt_out = f.read()
    assert q.read() == vt_out