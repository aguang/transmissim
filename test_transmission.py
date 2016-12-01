import transmission as trans
from rpy2.robjects import NA_Integer
from itertools import repeat

class TestBinary:
    def test_groupancestor(self):
        ances = [NA_Integer, 1, NA_Integer, 1, 1, 5, 5, 1]
        sources = trans.group_ancestors(ances)
        assert sources == {1:[1,3,4,7],5:[5,6]}
        assert NA_Integer not in sources

    def test_pairinfected(self):
        source = {1:[1,3,4,7],5:[5,6]}
        pi = trans.pair_infected(source)
        assert pi[0] == [(1, '1'), (2, '3'), (3, '4'), ('0', '7')]
        assert pi[4] == [(1, '5'), ('4','6')]

    def test_assignbl(self):
        pi = {0: [(1, '1'), (2, '3'), (3, '4'), ('0','7')], 4:[(1, '5'), ('4','6')]}
        nmut = [NA_Integer, 13, NA_Integer, 120, 143, 14, 122, 272]
        duration = 350
        bl = trans.assign_bl(pi, nmut, duration)
        assert bl[0] == [(337, 13), (230, 107), (207, 23), (78, 129)]
        assert bl[4] == [(193, 14), (85, 108)]

        duration = 500
        bl = trans.assign_bl(pi, nmut, duration)
        assert list(map(sum, bl[0])) == [500, 487, 380, 357]

    def test_fulltree(self):
        bl = {0: [(337, 13), (230, 107), (207, 23), (78, 129)], 4: [(193, 14), (85, 108)]}
        pi = {0: [(1, '1'), (2, '3'), (3, '4'), ('0','7')], 4:[(1, '5'), ('4','6')]}

        t = trans.full_tree(pi, bl)
        assert t == '((((0:78,7:78):129,((4:85,6:85):108,5:193):14):23,3:230):107,1:337):13'
