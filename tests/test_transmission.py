from transmissim import transmission as trans
from rpy2.robjects import NA_Integer
from itertools import repeat

class TestGroupAncestor:

    def test_groupancestor(self):
        ances = [NA_Integer, 1, NA_Integer, 1, 1, 5, 5, 1]
        sources, NA_list = trans.group_ancestors(ances)
        assert sources == {1:[1,3,4,7],5:[5,6]}
        assert NA_list == {0, 2}

    def test_groupancestor2(self):
        ances = [NA_Integer, 1, 2, 1, 2, 5, 4, 7, 1, 7, 8]
        sources, NA_list = trans.group_ancestors(ances)
        assert sources == {1:[1,3,8], 2:[2,4], 5:[5], 4:[6], 7:[7,9], 8:[10]}
        assert NA_list == {0}

    def test_groupancestor3(self):
        ances = [NA_Integer, 1, 2, 3, 3, 3, 5, 3, 2, 5, 1, 1, 10, 2, 11]
        sources, NA_list = trans.group_ancestors(ances)
        assert sources == {1: [1, 10, 11], 2: [2, 8, 13], 3: [3, 4, 5, 7], 5: [6, 9], 10: [12], 11: [14]}
        assert NA_list == {0}

    def test_multiple_clusters(self): # outbreaker with seed 998878
        ances = [NA_Integer, 1, NA_Integer, NA_Integer, 3, 5, NA_Integer, 6, NA_Integer, NA_Integer, NA_Integer, NA_Integer, NA_Integer, 9, NA_Integer, NA_Integer]
        sources, NA_list = trans.group_ancestors(ances)
        assert sources == {1: [1], 3: [4], 5: [5], 6: [7], 9: [13]}
        assert NA_list == {0,2,3,6,8,9,10,11,12,14,15}

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

    def test_pairinfected3(self):
        source = {1: [1, 10, 11], 2: [2, 8, 13], 3: [3, 4, 5, 7], 5: [6, 9], 10: [12], 11: [14]}
        pi = trans.pair_infected(source)
        assert pi[0] == [(1, '1'), (2, '10'), ('0', '11')]
        assert pi[1] == [(1, '2'), (2, '8'), ('1', '13')]
        assert pi[2] == [(1, '3'), (2, '4'), (3, '5'), ('2', '7')]
        assert 3 not in pi.keys()
        assert pi[4] == [(1, '6'), ('4', '9')]
        assert pi[9] == [('9', '12')]
        assert pi[10] == [('10', '14')]

    def test_multiple_clusters_pi1(self):
        source = {1: [1], 3: [4], 5: [5], 6: [7], 9: [13]}
        pi = trans.pair_infected(source)
        assert pi == {0: [('0', '1')], 8: [('8', '13')], 2: [('2', '4')], 4: [('4', '5')], 5: [('5', '7')]}

class TestFullTree:

    def test_fulltree(self):
        pi = {0: [(1, '1'), (2, '3'), (3, '4'), ('0','7')], 4:[(1, '5'), ('4','6')]}
        nmut = [NA_Integer, 13, NA_Integer, 120, 143, 14, 122, 272]
        NA_set = {0}
        t = trans.full_tree(pi, nmut, NA_set)
        assert t == ['((((0:1,7:1):272,((4:1,6:1):122,5:1):14):143,3:1):120,1:1):13;']

    def test_multiple_ft1(self):
        ances = [NA_Integer, 1, NA_Integer, NA_Integer, 3, 5, NA_Integer, 6, NA_Integer, NA_Integer, NA_Integer, NA_Integer, NA_Integer, 9, NA_Integer, NA_Integer]
        sources, NA_set = trans.group_ancestors(ances)
        pi = trans.pair_infected(sources)
        nmut = [NA_Integer, 7, NA_Integer, NA_Integer, 48, 1, NA_Integer, 11, NA_Integer, NA_Integer, NA_Integer, NA_Integer, NA_Integer, 54, NA_Integer, NA_Integer]
        t = trans.full_tree(pi, nmut, NA_set)
        assert len(t) == 11
        assert(t[0]) == '(0:1,1:1):7;'
        assert(t[1]) == '(2:1,(4:1,(5:1,7:1):11):1):48;'
        assert(t[2]) == '3:1;'
        assert(t[3]) == '6:1;'
        
