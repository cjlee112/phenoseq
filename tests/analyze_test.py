import unittest
import util

try:
    from phenoseq import analyze
except ImportError:
    util.add_root_path()
    from phenoseq import analyze

class Score_Test(unittest.TestCase):
    'basic scoring tests'

    def test_score_pooled(self):
        regionSNPs = {'foo': [0] * 10}
        regionXS = {'foo': 0.2}
        s = analyze.score_pooled(regionSNPs, regionXS)
        assert util.approx_equal(s[0][0], 2.35e-14)

    def test_score_group(self):
        regionSNPs = {'foo': [0] * 10, 'bar': [0,0]}
        regionXS = {'foo': 0.2, 'bar': 1.}
        groups = {'1':('foo', 'bar')}
        s = analyze.score_groups_pooled(regionSNPs, groups, regionXS)
        assert util.approx_equal(s[0][0], 6.17e-09)

    
