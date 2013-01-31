import unittest
import subprocess
import glob

class Command_Test(unittest.TestCase):
    'trivial command line tests; just checks whether they die!'

    def setUp(self, vcfPattern='data/aligned_*.vcf'):
        self.vcfFiles = glob.glob(vcfPattern)
        assert len(self.vcfFiles) > 0

    def test_analyze(self):
        subprocess.check_output(['phenoseq_analyze', '-g', 'data/NC_000913.gbk']
                                + self.vcfFiles)

    def test_analyze2(self):
        subprocess.check_output(['phenoseq_analyze', '-g', 'data/NC_000913.gbk',
                                 '--fastafile', 'data/NC_000913.fna']
                                + self.vcfFiles)

    def test_cost(self):
        subprocess.check_output(['phenoseq_cost', '20', '30', '1000',
                                 '50', '5'])

    def test_pathways(self):
        subprocess.check_output('phenoseq_pathways -g data/NC_000913.gbk -f data/func-associations.col -t data/genes.col'.split()
                                + self.vcfFiles)

    def test_kaks(self):
        subprocess.check_output('phenoseq_kaks -g data/NC_000913.gbk -f data/func-associations.col -t data/genes.col'.split()
                                + self.vcfFiles)

    def test_hypergeom(self):
        subprocess.check_output('phenoseq_hypergeom -n 30 -g data/NC_000913.gbk -f data/func-associations.col -t data/genes.col'.split()
                                + self.vcfFiles)

        
