from scipy import stats
from math import log
import glob
import re
import sys
try:
    from scipy import stats
except ImportError:
    pass
try:
    from Bio import SeqIO
except ImportError:
    pass
try:
    from pygr import seqdb, annotation, cnestedlist, sequtil
except ImportError:
    pass
    

''' simple usage example:

laptop:doe leec$ python2.6 analyze.py NC_000913.gbk [ACGT]*.vcf
reading gene annotations from NC_000913.gbk
reading tag files: ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']
scoring genes...
top 20 hits: [(6.6292398002665553e-23, 'acrB'), (9.5886861652296713e-09, 'marC'), (1.2420147995288613e-07, 'stfP'), (7.5888980168626557e-07, 'ykgC'), (2.4852332221987147e-06, 'aes'), (1.2076613799972474e-05, 'ampH'), (2.677716300626633e-05, 'paoC'), (2.7423615997657459e-05, 'nfrA'), (3.0707657321643245e-05, 'ydhB'), (8.2316128201056232e-05, 'yaiP'), (0.00011947463446331284, 'acrA'), (0.00017183498597744571, 'xanQ'), (0.00017801741636809529, 'ykgD'), (0.0002470441650715134, 'yegQ'), (0.00024848120894985465, 'yfjJ'), (0.00026020746746671289, 'yagX'), (0.00032216191302056855, 'pstA'), (0.00033538869893896079, 'prpE'), (0.00035039080759966868, 'mltF'), (0.00044366190311229754, 'purE')]

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

    def get_snps(self):
        'get all filtered SNPs sorted in ascending positional order'
        l = []
        for tag, rs in self.items():
            for snp in rs:
                l.append((snp.pos, tag, snp))
        l.sort()
        return l


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


def read_genome_annots(gbfile, iseq=0, featureType='CDS'):
    features = list(SeqIO.parse(gbfile, 'genbank'))[iseq].features
    genome = seqdb.SequenceFileDB(gbfile.split('.')[0] + '.fna')
    seqID = genome.keys()[iseq]
    annodb = annotation.AnnotationDB({}, genome,
                                     sliceAttrDict=dict(id=0, start=1, stop=2,
                                                        orientation=3))
    for f in features:
        if f.type == featureType:
            try:
                name = f.qualifiers['gene'][0]
            except KeyError:
                pass
            else:
                annodb.new_annotation(name,
                    (seqID, f.location.start.position,
                     f.location.end.position, f.strand))
    al = cnestedlist.NLMSA('tmp', 'memory', pairwiseMode=True)
    for a in annodb.itervalues():
        al.addAnnotation(a)
    al.build()
    return annodb, al, genome[seqID]


class GeneSNPDict(dict):
    def __init__(self, tagDict, annodb, al, dna):
        dict.__init__(self)
        self.tagDict = tagDict
        self.annodb = annodb
        self.al = al
        self.dna = dna
        self.nAT = self.nGC = weight = 0
        rc = dict(A='T', C='G', G='C', T='A')
        for pos, tag, snp in tagDict.get_snps():
            try:
                geneSNP = al[dna[pos - 1:pos]].keys()[0]
            except IndexError:
                pass
            else:
                if geneSNP.orientation < 0:
                    geneSNP = -geneSNP
                gene = geneSNP.id
                ipos = geneSNP.start % 3
                codonStart = geneSNP.start - ipos
                codon = str(geneSNP.path.sequence[codonStart:codonStart + 3])
                b = snp.alt # substitution letter
                if geneSNP.sequence.orientation < 0: # must complement
                    b = rc[b.upper()]
                codonAlt = codon[:ipos] + b + codon[ipos + 1:]
                snp.aaRef = sequtil.translate_orf(codon)
                snp.aaAlt = sequtil.translate_orf(codonAlt)
                if snp.aaRef != snp.aaAlt: # only count non-synonymous muts
                    weight = 1
                    self.setdefault(gene, []).append((tag,snp))
                else:
                    weight = 0
            if snp.ref in 'GCgc':
                self.nGC += weight
            else:
                self.nAT += weight

    def get_scores(self, gcTotal=None, atTotal=None, geneGCDict=None):
        'will use pre-computed GC/AT count data if you provide it'
        if gcTotal is None:
            gcTotal, atTotal = calc_gc(str(self.dna))
        results = []
        for k, v in self.items():
            if geneGCDict:
                gcLen, atLen = geneGCDict[k]
            else:
                gcLen, atLen = calc_gene_gc(self.annodb, k)
            pois = stats.poisson(self.nGC * float(gcLen) / gcTotal +
                                 self.nAT * float(atLen) / atTotal)
            results.append((pois.sf(len(v) - 1), k))
        results.sort()
        return results

def calc_gc(s):
    'calc GC and AT counts in string s'
    s = s.upper()
    gcTotal = s.count('G') + s.count('C')
    atTotal = len(s) - gcTotal
    return gcTotal, atTotal

def calc_gene_gc(annodb, geneID):
    'get GC, AT counts for specified gene'
    ann = annodb[geneID]
    return calc_gc(str(ann.sequence))

def generate_subsets(tagFiles, annodb, al, dna):
    'generate results from all possible subsets of tagFiles'
    d = {}
    gcTotal, atTotal = calc_gc(str(dna))
    geneGCDict = {}
    for geneID in annodb:
        geneGCDict[geneID] = calc_gene_gc(annodb, geneID)
    for i in range(1, pow(2, len(tagFiles))):
        print 'subset', i
        l = [tagFile for (j, tagFile) in enumerate(tagFiles)
             if i & pow(2, j)]
        tagDict = read_tag_files(l)
        gsd = GeneSNPDict(tagDict, annodb, al, dna)
        results = gsd.get_scores(gcTotal, atTotal, geneGCDict)
        d[tuple([s.split('.')[0] for s in l])] = results
    return d

if __name__ == '__main__':
    print 'reading gene annotations from', sys.argv[1]
    annodb, al, dna = read_genome_annots(sys.argv[1])
    tagFiles = sys.argv[2:]
    print 'reading tag files:', tagFiles
    tagDict = read_tag_files(tagFiles)
    gsd = GeneSNPDict(tagDict, annodb, al, dna)
    print 'scoring genes...'
    results = gsd.get_scores()
    print 'top 20 hits:', results[:20]
    
