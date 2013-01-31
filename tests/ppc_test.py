import unittest
import util

try:
    from phenoseq import analyze
except ImportError:
    util.add_root_path()
    from phenoseq import analyze

vcfFiles = ('data/aligned_s_8_ACAGTG.vcf',
            'data/aligned_s_8_ACTTGA.vcf',
            'data/aligned_s_8_ATCACG.vcf',
            'data/aligned_s_8_CAGATC.vcf',
            'data/aligned_s_8_CGATGT.vcf',
            'data/aligned_s_8_GCCAAT.vcf',
            'data/aligned_s_8_TGACCA.vcf',
            'data/aligned_s_8_TTAGGC.vcf')

class Gene_Test(unittest.TestCase):
    'basic scoring tests'

    def setUp(self, gbfile='data/NC_000913.gbk'):
        annodb, al, genome = analyze.read_genbank_annots(gbfile)
        self.annodb = annodb
        self.al = al
        self.genome = genome
        self.snps = analyze.read_tag_files(vcfFiles)
        self.gsd = analyze.map_snps(self.snps, al, genome)
        self.gsdNS = analyze.filter_nonsyn(self.gsd)

    def test_gbdata(self):
        assert len(self.annodb) == 4244
        assert len(self.genome) == 1
        assert len(self.genome.values()[0]) == 4639675

    def test_snp_reading(self):
        assert len(self.snps) == 2433

    def test_snp_filtering(self):
        assert len(self.gsd) == 1355
        n = sum([len(t) for t in self.gsd.values()])
        assert n == 2157
        assert len(self.gsdNS) == 1016
        n = sum([len(t) for t in self.gsdNS.values()])
        assert n == 1448
        
    def test_gene_scores(self):
        s = analyze.score_genes_pooled(self.gsdNS, genome=self.genome,
                                       annodb=self.annodb)
        assert util.approx_equal(s[0][0], 1.36e-25)
        assert util.approx_equal(s[1][0], 8.30e-14)
        assert util.approx_equal(s[2][0], 4.77e-04)
        assert util.approx_equal(s[3][0], 0.0452)

    def test_simple_group(self):
        groups = {'foo': ['iclR', 'aceK']}
        s = analyze.score_groups_pooled(self.gsdNS, groups,
                                        genome=self.genome,
                                        annodb=self.annodb)
        assert util.approx_equal(s[0][0], 3.58e-42)
