import transmission as trans
from rpy2.robjects import NA_Integer

class TestBinary:
    def test_groupancestor(self):
        ances = [NA_Integer, 1, NA_Integer, 1, 1, 5, 5, 1]
        sources = trans.group_ancestors(ances)
        assert sources == {1:[1,3,4,7],5:[5,6]}
        assert NA_Integer not in sources

    def test_pairinfected(self):
        source = {1:[1,3,4,7],5:[5,6]}
        pi = trans.pair_infected(source)
        assert pi[0] == [(1, '1'), (2, '3'), (3, '4'), ('7', '0')]
        assert pi[4] == [(1, '5'), ('6','4')]
