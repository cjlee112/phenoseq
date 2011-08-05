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
    try:
        if s.find('.') >= 0 or s.lower().find('e') >= 0:
            return float(s)
        else:
            return int(s)
    except ValueError:
        return s


class SNP(object):
    'represent VCF fields as object attributes, properly unpacked'
    def __init__(self, colnames, fields, add_attrs=None, **kwargs):
        for i,k in enumerate(colnames):
            if i >= len(fields):
                continue
            if k == 'INFO':
                for s in fields[i].split(';'):
                    try:
                        k,v = s.split('=')
                    except ValueError: # VCF files contain DS with no value?
                        k = s
                        v = None
                        continue
                    l = v.split(',')
                    if len(l) > 1:
                        setattr(self, k.lower(),
                                [convert_number(x) for x in l])
                    else:
                        setattr(self, k.lower(), convert_number(v))
            elif k == 'FORMAT' and i + 2 <= len(fields):
                vals = fields[i + 1].split(':')
                for j,k in enumerate(fields[i].split(':')):
                    setattr(self, k.lower(), convert_number(vals[j]))
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

def read_vcf(path, vcffields=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                              'FILTER', 'INFO', 'FORMAT'), **kwargs):
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

class SNPSet(object):
    'filter SNP candidates from one lane using filter expression'
    def __init__(self, vcfFile,
                 filterExpr='snp.af1 <= 0.5 and getattr(snp, "pv4", (0.,))[0] >= 0.01 and snp.qual > 90'):
        self.snps = read_vcf(vcfFile, add_attrs=add_snp_attrs)
        self.filterExpr = filterExpr

    def __iter__(self):
        for snp in self.snps:
            if eval(self.filterExpr):
                yield snp
        

class ReplicateSet(object):
    'filters SNP candidates using replicate lanes'
    def __init__(self, mergedFile, replicateFiles,
                 filterExpr='snp.af1 <= 0.5 and len(self[snp]) >= 2'):
        self.snps = read_vcf(mergedFile, add_attrs=add_snp_attrs)
        self.filterExpr = filterExpr
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
            if eval(self.filterExpr):
                yield snp


class TagDict(dict):
    'keys are tags, values are ReplicateSet objects'
    def __init__(self, tagDict, containerClass=SNPSet):
        dict.__init__(self)
        d = {}
        for tag,args in tagDict.items():
            self[tag] = containerClass(*args)

    def get_snps(self):
        'get all filtered SNPs sorted in ascending positional order'
        l = []
        for tag, rs in self.items():
            for snp in rs:
                l.append((snp.pos, tag, snp))
        l.sort()
        return l


def get_filestem(s):
    return s.split('.')[0]

def get_replicate_files(s):
    return glob.glob('*_' + s + '.vcf')

def read_tag_files(tagFiles, tagfunc=get_filestem, *args):
    'construct TagDict for multiple tagFiles'
    d = {}
    for tagFile in tagFiles:
        tag = tagfunc(tagFile)
        d[tag] = (tagFile,) + args
    return TagDict(d)

def read_replicate_files(tagFiles, tagfunc=get_filestem,
                         replicatefunc=get_replicate_files,
                         *args):
    'construct TagDict for multiple tags each with replicate files'
    d = {}
    for tagFile in tagFiles:
        tag = tagfunc(tagFile)
        repFiles = replicatefunc(tag)
        d[tag] = (tagFile, repFiles) + args
    return TagDict(d, containerClass=ReplicateSet)


def read_vcf_singleton(vcfFile):
    'load dataset consisting of a single vcf file'
    d = dict(mytag=(vcfFile, (), 0))
    return TagDict(d)


def read_genome_annots(gbfile, fastafile=None, iseq=0, featureType='CDS'):
    'construct annotation DB for gene coding regions in a genome'
    try:
        gbparse = SeqIO.parse(gbfile, 'genbank')
    except TypeError: # SeqIO changed its interface?
        ifile = open(gbfile)
        try:
            gbparse = SeqIO.parse(ifile, 'genbank')
            features = list(gbparse)[iseq].features
        finally:
            ifile.close()
    else:
        features = list(gbparse)[iseq].features
    if fastafile is None:
        fastafile = gbfile.split('.')[0] + '.fna'
    genome = seqdb.SequenceFileDB(fastafile)
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
    def __init__(self, tagDict, annodb, al, dna, count_syn=False):
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
                if count_syn or snp.aaRef != snp.aaAlt: # count this SNP
                    weight = 1
                    self.setdefault(gene, []).append((tag,snp))
                else:
                    weight = 0
            if snp.ref in 'GCgc':
                self.nGC += weight
            else:
                self.nAT += weight

    def get_scores(self, gcTotal=None, atTotal=None, geneGCDict=None,
                   useBonferroni=True):
        'will use pre-computed GC/AT count data if you provide it'
        if useBonferroni:
            correction = len(self)
        else:
            correction = 1.
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
            results.append((correction * pois.sf(len(v) - 1), k))
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


def get_gc_totals(annodb, dna):
    'get GC/AT counts for whole genome and for individual genes in annodb'
    gcTotal, atTotal = calc_gc(str(dna))
    geneGCDict = {}
    for geneID in annodb:
        geneGCDict[geneID] = calc_gene_gc(annodb, geneID)
    return gcTotal, atTotal, geneGCDict


def gc_gene_sizes(annodb, dna, gcBias=36.):
    'return GC-bias weighted sizes for whole genome and each gene'
    gcTotal, atTotal, geneGCDict = get_gc_totals(annodb, dna)
    total = gcTotal * gcBias + atTotal 
    l = []
    for gc, at in geneGCDict.values():
        l.append(gc * gcBias + at)
    return total, l

def generate_subsets(tagFiles, annodb, al, dna, nbest=None, *args):
    'generate results from all possible subsets of tagFiles'
    d = {}
    gcTotal, atTotal, geneGCDict = get_gc_totals(annodb, dna)
    for i in range(1, pow(2, len(tagFiles))):
        print 'subset', i
        l = [tagFile for (j, tagFile) in enumerate(tagFiles)
             if i & pow(2, j)]
        tagDict = read_tag_files(l, *args)
        gsd = GeneSNPDict(tagDict, annodb, al, dna)
        results = gsd.get_scores(gcTotal, atTotal, geneGCDict)
        d[tuple([s.split('.')[0] for s in l])] = results[:nbest]
    return d


def get_pool_tags(poolFile):
    'extract tags from filenames of the form tag_*_tag_*...'
    l = poolFile.split('_')
    return tuple([l[i] for i in range(0, len(l), 2)])


def process_pools(poolFiles, annodb, al, dna, tagFunc=get_pool_tags,
                  nbest=None, *args):
    'generate results for a list of pool files'
    d = {}
    gcTotal, atTotal, geneGCDict = get_gc_totals(annodb, dna)
    for poolFile in poolFiles:
        tagDict = read_vcf_singleton(poolFile)
        gsd = GeneSNPDict(tagDict, annodb, al, dna)
        results = gsd.get_scores(gcTotal, atTotal, geneGCDict)
        d[tagFunc(poolFile)] = results[:nbest]
    return d


def enumerate_pool_subsets(poolFiles, n, tagFunc=get_pool_tags, d=None):
    'recursively generate all disjoint pool combinations containing n pools'
    if d is None:
        d = set()
    for i,poolFile in enumerate(poolFiles):
        tags = get_pool_tags(poolFile)
        if [1 for t in tags if t in d]: # check for overlap
            continue # can't use this pool, due to overlap
        if n == 1: # terminate the recursion
            yield (poolFile,)
        else: # recurse
            d.update(tags)
            for result in enumerate_pool_subsets(poolFiles[i + 1:], n - 1,
                                                 tagFunc, d):
                yield (poolFile,) + result
            for tag in tags:
                d.remove(tag)


def generate_pool_subsets(poolFiles, npools, annodb, al, dna,
                          tagfunc=lambda s:s,
                          replicatefunc=lambda s:(), minRep=0,
                          nbest=None, *args):
    'run phenoseq analysis on all disjoint pool combinations'
    d = {}
    gcTotal, atTotal, geneGCDict = get_gc_totals(annodb, dna)
    for i in npools:
        for pfiles in enumerate_pool_subsets(poolFiles, i):
            print 'analyzing', pfiles
            tagDict = read_tag_files(pfiles, tagfunc, replicatefunc,
                                     minRep, *args)
            gsd = GeneSNPDict(tagDict, annodb, al, dna)
            results = gsd.get_scores(gcTotal, atTotal, geneGCDict)
            k = reduce(lambda x,y:x + y, [get_pool_tags(t) for t in pfiles])
            d[k] = results[:nbest]
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
    
