from scipy import stats
from math import log

def convert_number(s):
    'convert string to number if possible'
    if not s[0].isdigit():
        return s
    elif s.find('.') >= 0:
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
