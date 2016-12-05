import transmission as trans
from rpy2.robjects import NA_Integer
from itertools import repeat

class TestGroupAncestor:

    def test_groupancestor(self):
        ances = [NA_Integer, 1, NA_Integer, 1, 1, 5, 5, 1]
        sources = trans.group_ancestors(ances)
        assert sources == {1:[1,3,4,7],5:[5,6]}
        assert NA_Integer not in sources

    def test_groupancestor2(self):
        ances = [NA_Integer, 1, 2, 1, 2, 5, 4, 7, 1, 7, 8]
        sources = trans.group_ancestors(ances)
        assert sources == {1:[1,3,8], 2:[2,4], 5:[5], 4:[6], 7:[7,9], 8:[10]}
        assert NA_Integer not in sources

class TestPairInfected:

    def test_pairinfected(self):
        source = {1:[1,3,4,7],5:[5,6]}
        pi = trans.pair_infected(source)
        assert pi[0] == [(1, '1'), (2, '3'), (3, '4'), ('0', '7')]
        assert pi[4] == [(1, '5'), ('4','6')]

    def test_pairinfected2(self):
        source = {1:[1,3,8], 2:[2,4], 5:[5], 4:[6], 7:[7,9], 8:[10]}
        pi = trans.pair_infected(source)
        assert pi[0] == [(1, '1'), (2, '3'), ('0', '8')]
        assert pi[1] == [(1, '2'), ('1', '4')]
        assert pi[4] == [('4', '5')]
        assert pi[3] == [('3', '6')]
        assert pi[6] == [(1, '7'), ('6', '9')]
        assert pi[7] == [('7', '10')]

class TestFullTree:

    def test_fulltree(self):
        pi = {0: [(1, '1'), (2, '3'), (3, '4'), ('0','7')], 4:[(1, '5'), ('4','6')]}
        nmut = [NA_Integer, 13, NA_Integer, 120, 143, 14, 122, 272]
        t = trans.full_tree(pi, nmut)
        assert t == '((((0:1,7:1):272,((4:1,6:1):122,5:1):14):143,3:1):120,1:1):13'
