from scipy import stats
from math import log
import glob

''' simple usage example:

laptop:doe leec$ python2.6 -i analyze.py [ACGT]*.vcf
>>> len(tagDict)
10
>>> sum([len(list(rs)) for rs in tagDict.values()])/32.
128.09375

'''

def convert_number(s):
    'convert string to number if possible'
    if not s[0].isdigit():
        return s
    elif s.find('.') >= 0 or s.lower().find('e') >= 0:
        return float(s)
    else:
        return int(s)


class SNP(object):
    'represent VCF fields as object attributes, properly unpacked'
    def __init__(self, colnames, fields, add_attrs=None, **kwargs):
        for i,k in enumerate(colnames):
            if i >= len(fields):
                continue
            if k == 'INFO':
                for s in fields[i].split(';'):
                    k,v = s.split('=')
                    l = v.split(',')
                    if len(l) > 1:
                        setattr(self, k.lower(),
                                [convert_number(x) for x in l])
                    else:
                        setattr(self, k.lower(), convert_number(v))
            else:
                setattr(self, k.lower(), convert_number(fields[i]))
        for k,v in kwargs.items():
            setattr(self, k, v)
        if add_attrs is not None:
            add_attrs(self)

def add_snp_attrs(self):
    'add convenience attributes for SNP'
    self.nreads = sum(self.dp4)
    self.nalt = sum(self.dp4[2:])

def read_vcf(path, **kwargs):
    'read VCF file and return list of SNP objects'
    ifile = open(path)
    l = []
    for line in ifile:
        line = line.strip()
        if line[:2] == '##':
            continue
        elif line[0] == '#':
            vcffields = line[1:].split('\t')
        else:
            l.append(SNP(vcffields, line.split('\t'), **kwargs))
    ifile.close()
    return l

def npool_posterior(snps, npools=(3,4)):
    'assess likelihood of data under different pooling factors'
    logP = [0.] * len(npools)
    for snp in snps:
        for i,npool in enumerate(npools):
            pmf = stats.binom(sum(snp.dp4), 1. / npool)
            logP[i] += log(pmf.pmf(sum(snp.dp4[2:])))
    return logP

def basic_snp_filter(snp):
    return snp.qual >= 99 and snp.af1 <= 0.5

def get_binom_cutoff(c, p, s, nmax=1):
    'find minimal cutoff that gives at most nmax false positives per genome'
    pmf = stats.binom(c, p)
    l = 0
    r = c
    while r > l + 1:
        mid = (l + r) / 2
        if pmf.sf(mid - 1) * s > nmax:
            l = mid
        else:
            r = mid
    return r

def binom_filter(snp, p=0.01, s=4e+06):
    'return True if snp passes sequencing error cutoff'
    return snp.nalt >= get_binom_cutoff(snp.nreads, p, s)

class ReplicateSet(object):
    'filters SNP candidates using replicate lanes'
    def __init__(self, mergedFile, replicateFiles, minRep=2, maxF=0.5):
        self.snps = read_vcf(mergedFile, add_attrs=add_snp_attrs)
        self.minRep = minRep
        self.maxF = maxF
        self.replicates = []
        d = {}
        for path in replicateFiles: # read replicate lanes
            snps = read_vcf(path, add_attrs=add_snp_attrs)
            self.replicates.append(snps)
            for snp in snps: # index SNP records from replicate lanes
                d.setdefault((snp.pos, snp.alt), []).append(snp)
        self.repMap = d

    def __getitem__(self, snp):
        'get repMap records matching this snp pos and letter'
        return self.repMap.get((snp.pos, snp.alt), ())

    def __iter__(self):
        'generate all snps matching filter criteria'
        for snp in self.snps:
            if snp.af1 <= self.maxF and len(self[snp]) >= self.minRep:
                yield snp


class TagDict(dict):
    'keys are tags, values are ReplicateSet objects'
    def __init__(self, tagDict):
        dict.__init__(self)
        d = {}
        for tag,args in tagDict.items():
            self[tag] = ReplicateSet(*args)

def read_tag_files(tagFiles, tagfunc=lambda s:s.split('.')[0],
                      replicatefunc=lambda s:glob.glob('*_' + s + '.vcf'),
                      *args):
    'construct TagDict for all these tagFiles'
    d = {}
    for tagFile in tagFiles:
        tag = tagfunc(tagFile)
        repFiles = replicatefunc(tag)
        d[tag] = (tagFile, repFiles) + args
    return TagDict(d)

if __name__ == '__main__':
    import sys
    tagFiles = sys.argv[1:]
    tagDict = read_tag_files(tagFiles)
