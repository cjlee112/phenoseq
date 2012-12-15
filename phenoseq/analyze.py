from math import log
import glob
import re
import sys
import warnings
try:
    from scipy import stats
except ImportError:
    warnings.warn('scipy.stats not found.  Some phenoseq functions will not work.')
try:
    from Bio import SeqIO
except ImportError:
    warnings.warn('biopython not found.  Some phenoseq functions will not work.')
try:
    from pygr import seqdb, annotation, cnestedlist, sequtil
except ImportError:
    warnings.warn('pygr not found.  Some phenoseq functions will not work.')
    

''' simple usage example:

laptop:doe leec$ python2.6 analyze.py NC_000913.gbk [ACGT]*.vcf
reading gene annotations from NC_000913.gbk
reading tag files: ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']
scoring genes...
top 20 hits: [(6.6292398002665553e-23, 'acrB'), (9.5886861652296713e-09, 'marC'), (1.2420147995288613e-07, 'stfP'), (7.5888980168626557e-07, 'ykgC'), (2.4852332221987147e-06, 'aes'), (1.2076613799972474e-05, 'ampH'), (2.677716300626633e-05, 'paoC'), (2.7423615997657459e-05, 'nfrA'), (3.0707657321643245e-05, 'ydhB'), (8.2316128201056232e-05, 'yaiP'), (0.00011947463446331284, 'acrA'), (0.00017183498597744571, 'xanQ'), (0.00017801741636809529, 'ykgD'), (0.0002470441650715134, 'yegQ'), (0.00024848120894985465, 'yfjJ'), (0.00026020746746671289, 'yagX'), (0.00032216191302056855, 'pstA'), (0.00033538869893896079, 'prpE'), (0.00035039080759966868, 'mltF'), (0.00044366190311229754, 'purE')]

'''


#######################################################################
# VCF reading functions and SNP object class

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


def get_attrname(k):
    'cast all-uppercase attr names to lowercase'
    if k == k.upper():
        return k.lower()
    else:
        return k

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
                    k = get_attrname(k)
                    l = v.split(',')
                    if len(l) > 1:
                        setattr(self, k, [convert_number(x) for x in l])
                    else:
                        setattr(self, k, convert_number(v))
            elif k == 'FORMAT' and i + 2 <= len(fields):
                vals = fields[i + 1].split(':')
                for j,k in enumerate(fields[i].split(':')):
                    setattr(self, get_attrname(k), convert_number(vals[j]))
            else:
                setattr(self, get_attrname(k), convert_number(fields[i]))
        for k,v in kwargs.items():
            setattr(self, k, v)
        if add_attrs is not None:
            add_attrs(self)

    def __repr__(self):
        return '<SNP ' + str(self.chrom) + ':' + str(self.pos) + ':' + self.ref \
            + ':' + self.alt + '>'

    def __hash__(self):
        return hash((self.chrom, self.pos, self.alt))

    def __cmp__(self, other):
        return cmp((self.chrom, self.pos, self.alt),
                   (other.chrom, other.pos, other.alt))


def add_snp_attrs(self):
    'add convenience attributes for SNP'
    self.nreads = sum(self.dp4)
    self.nalt = sum(self.dp4[2:])

def read_vcf(path, vcffields=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                              'FILTER', 'INFO', 'FORMAT'), **kwargs):
    'read VCF file and return list of SNP objects'
    ifile = open(path)
    l = []
    try:
        for line in ifile:
            line = line.strip()
            if line[:2] == '##':
                continue
            elif line[0] == '#':
                vcffields = line[1:].split('\t')
            else:
                l.append(SNP(vcffields, line.split('\t'), **kwargs))
    finally:
        ifile.close()
    return l




#######################################################################
# snp false positive filtering functions

def filter_snps(snps,
                filterExpr='snp.af1 <= 0.5 and getattr(snp, "pv4", (0.,))[0] >= 0.01 and snp.qual > 90'):
    'filter SNPs using attribute expressions'
    for snp in snps:
        if eval(filterExpr):
            yield snp

def filter_snps_repfiles(snps, replicateFiles,
                         filterExpr='snp.af1 <= 0.5 and len(get_replicates(snp)) >= 2'):
    'filters SNP candidates using replicate lanes'
    replicates = []
    repMap = {}
    for path in replicateFiles: # read replicate lanes
        repSnps = read_vcf(path, add_attrs=add_snp_attrs)
        replicates.append(repSnps)
        for snp in repSnps: # index SNP records from replicate lanes
            repMap.setdefault(snp, []).append(snp)
    def get_replicates(s):
        return repMap.get(s, ())
        
    for snp in snps:
        if eval(self.filterExpr):
            yield snp




#######################################################################
# utilities for processing multiple tagged libraries

def get_filestem(s):
    return s.split('.')[0]

def get_replicate_files(s):
    return glob.glob('*_' + s + '.vcf')

def read_tag_files(tagFiles, tagfunc=get_filestem, filterFunc=filter_snps,
                   replicatefunc=None, add_attrs=add_snp_attrs, *args, **kwargs):
    'merge SNPs from multiple tagFiles, with optional filtering'
    l = []
    for tagFile in tagFiles:
        tag = tagfunc(tagFile)
        snps = read_vcf(tagFile, tag=tag, add_attrs=add_attrs, **kwargs)
        if replicatefunc:
            repFiles = replicatefunc(tag)
            snps = list(filterFunc(snps, repFiles, *args))
        elif filterFunc:
            snps = list(filterFunc(snps, *args))
        l += snps
    return l






#######################################################################
# gene annotation reading

def read_genbank_annots(gbfile, fastafile=None, iseq=0, featureType='CDS'):
    '''construct annotation DB for gene CDS intervals.
    NB: this assumes each gene consists of ONE interval.
    This cannot be used for multi-exon genes!'''
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


def read_known_genes(path, sep='\t'):
    '''read UCSC knownGene file to obtain exon dict whose values are gene ID
    (first transcript ID in that gene)'''
    d = {}
    genes = {}
    trLen = {}
    ifile = open(path)
    try:
        for line in ifile:
            line = line.strip()
            t = line.split(sep)
            trID, seqID = t[:2]
            if t[2] == '-':
                ori = -1
            else:
                ori = 1
            geneID = trID
            stops = t[9].split(',')
            exons = []
            length = 0
            for i,start in enumerate(t[8].split(',')[:-1]):
                start = int(start)
                stop = int(stops[i])
                length += stop - start
                exons.append((start, stop))
                try: # see if exon already stored
                    geneID = d[(seqID, ori, start, stop)]
                except KeyError:
                    pass
            for exon in exons: # save exons as belonging to this gene
                d[(seqID, ori, exon[0], exon[1])] = geneID
            genes.setdefault(geneID, []).append(trID)
            trLen[trID] = length
    finally:
        ifile.close()
    return d, genes, trLen

def get_gene_maxlengths(genes, trLen):
    'get dict of gene lengths based on longest transcript'
    d = {}
    for gene, transcripts in genes.items():
        d[gene] = max([trLen[txID] for txID in transcripts])
    return d

def read_exon_annots(genome, genesFile='knownGene.txt'):
    '''read multi-exon transcript set and build exon annotation db
    and exon-to-gene mapping'''
    exonDict, genes, trLen = read_known_genes(genesFile)
    geneLengths = get_gene_maxlengths(genes, trLen)
    totalSize = sum(geneLengths.values())
    annodb = annotation.AnnotationDB({}, genome,
                                     sliceAttrDict=dict(id=0, orientation=1,
                                                        start=2, stop=3))
    al = cnestedlist.NLMSA('tmp', 'memory', pairwiseMode=True,
                           maxlen=1000000000)
    i = 0
    exonGene = {}
    for t,geneID in exonDict.iteritems():
        a = annodb.new_annotation(i, t)
        exonGene[i] = geneID
        i += 1
        al.addAnnotation(a)
    al.build()
    return annodb, al, exonGene, totalSize, geneLengths



#######################################################################
# basic composition analysis functions

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




#######################################################################
# snp mapping and impact filtering functions

def map_snps(snps, al, genome, exonGene):
    '''map snps to genes, returning dict of {gene:[snp1, snp2, ...]}.
    NB: this assumes al contains mapping to exon annotations, whose
    IDs can be mapped to gene IDs by exonGene.'''
    geneSNP = {}
    for snp in snps:
        ival = genome[snp.chrom][snp.pos - 1:snp.pos]
        try:
            a = al[ival].keys()[0]
        except IndexError:
            pass
        else:
            geneID = exonGene[a.id]
            geneSNP.setdefault(geneID, []).append(snp)
    return geneSNP


def map_snps_chrom1(snps, al, dna):
    '''map snps to genes, returning dict of {gene:[snp1, snp2, ...]}.
    NB: this assumes all snps map on a single chromosome!'''
    d = {}
    for snp in snps:
        try:
            geneSNP = al[dna[snp.pos - 1:snp.pos]].keys()[0]
        except IndexError:
            pass
        else:
            if geneSNP.orientation < 0:
                geneSNP = -geneSNP
            gene = geneSNP.id
            snp.geneSNP = geneSNP # record its position in the gene
            d.setdefault(gene, []).append(snp)
    return d

rcDict = dict(A='T', C='G', G='C', T='A')

def is_nonsyn(snp):
    'test whether snp is non-synonymous, assuming frame 0 CDS annotation'
    geneSNP = snp.geneSNP
    ipos = geneSNP.start % 3
    codonStart = geneSNP.start - ipos
    codon = str(geneSNP.path.sequence[codonStart:codonStart + 3])
    b = snp.alt # substitution letter
    if geneSNP.sequence.orientation < 0: # must complement
        b = rcDict[b.upper()]
    codonAlt = codon[:ipos] + b + codon[ipos + 1:]
    snp.aaRef = sequtil.translate_orf(codon)
    snp.aaAlt = sequtil.translate_orf(codonAlt)
    return snp.aaRef != snp.aaAlt
        
def filter_nonsyn(geneSNPdict):
    'filter gene snps to just non-synonymous subset'
    d = {}
    for gene,snps in geneSNPdict.items():
        snps = filter(is_nonsyn, snps)
        if snps: # only save non-empty list
            d[gene] = snps
    return d

def get_gene_na_ns(geneSNPdict):
    'get dict with counts (Na,Ns) for each gene'
    d = {}
    for gene,snps in geneSNPdict.items():
        Na = len(filter(is_nonsyn, snps))
        d[gene] = (Na, len(snps) - Na)
    return d


#######################################################################
# scoring functions

def score_pooled(regionSNPs, regionXS, useBonferroni=True):
    '''Uses poisson scoring for pooled library data (i.e. where
    multiple samples are pooled in each library).'''
    if useBonferroni:
        correction = len(regionSNPs)
    else:
        correction = 1.
    results = []
    for k,snps in regionSNPs.items():
        pois = stats.poisson(regionXS[k])
        results.append((correction * pois.sf(len(snps) - 1), k))
    results.sort()
    return results


class GeneCrossSection(object):
    'object acts as dict for computing gene cross-section'
    def __init__(self, geneSNPdict, gcTotal=None, atTotal=None,
                 geneGCDict=None, dnaseq=None, annodb=None):
        self.geneSNPdict = geneSNPdict
        self.geneGCDict = geneGCDict
        self.annodb = annodb
        self.dnaseq = dnaseq
        if gcTotal is None:
            gcTotal, atTotal = calc_gc(str(dnaseq))
        nGC = nAT = 0
        for snps in geneSNPdict.values(): # GC vs AT SNP totals
            for snp in snps:
                if snp.ref in 'GCgc':
                    nGC += 1
                else:
                    nAT += 1
        self.gcTotal = gcTotal
        self.atTotal = atTotal
        self.nGC = nGC
        self.nAT = nAT

    def __getitem__(self, k):
        'Compute effective mutation cross-section for a gene'
        try:
            gcLen, atLen = self.geneGCDict[k]
        except (TypeError, KeyError):
            gcLen, atLen = calc_gene_gc(self.annodb, k)
        return float(self.nGC * gcLen) / self.gcTotal \
                + float(self.nAT * atLen) / self.atTotal


def score_genes_pooled(geneSNPdict, useBonferroni=True,
                       **kwargs):
    '''Uses poisson scoring for pooled library data (i.e. where
    multiple samples are pooled in each library).
    Will use pre-computed GC/AT count data if you provide it'''
    geneXS = GeneCrossSection(geneSNPdict, **kwargs)
    return score_pooled(geneSNPdict, geneXS,
                        useBonferroni=useBonferroni)

def get_group_xs_snps(geneSNPdict, geneXS, groupDict):
    'combine snps and cross-sections for genes in each group'
    groupXS = {}
    groupSNPs = {}
    for k,genes in groupDict.items():
        snps = []
        xs = 0.
        for gene in genes:
            try:
                xs += geneXS[gene]
            except KeyError:
                warnings.warn('ignoring missing gene annotation: '
                              + gene)
            else:
                snps += geneSNPdict.get(gene, [])
        groupXS[k] = xs
        groupSNPs[k] = snps
    return groupSNPs, groupXS


def score_groups_pooled(geneSNPdict, groupDict, geneXS=None,
                        useBonferroni=True, **kwargs):
    'poisson p-value for gene groups'
    if not geneXS:
        geneXS = GeneCrossSection(geneSNPdict, **kwargs)
    groupSNPs, groupXS = get_group_xs_snps(geneSNPdict, geneXS,
                                           groupDict)
    return score_pooled(groupSNPs, groupXS,
                        useBonferroni=useBonferroni)

def score_genes(geneSNPdict, nsample, totalSize, geneLengths,
                useBonferroni=True):
    '''Uses binomial scoring for unpooled library data (i.e. each library
    tag is a single sample.'''
    if useBonferroni:
        correction = len(geneSNPdict)
    else:
        correction = 1.
    nSNP = sum([len(l) for l in geneSNPdict.values()]) # total in all samples
    mu = float(nSNP) / (nsample * totalSize) # SNP density per exome nucleotide
    geneScores = []
    for geneID, snps in geneSNPdict.items():
        d = {}
        for snp in snps: # only count one snp per sample
            d[snp.tag] = 0
        pois = stats.poisson(mu * geneLengths[geneID]) # pmf for one sample
        pmf = stats.binom(nsample, pois.sf(0))
        geneScores.append((correction * pmf.sf(len(d) - 1), geneID))
    geneScores.sort()
    return geneScores



#######################################################################
# sub-experiment analysis utilities
# for generating subsets of the experimental data and analyzing
# results that would be obtained from those subsets.

def generate_subsets(tagFiles, annodb, al, dna, nbest=None, *args):
    'generate results from all possible subsets of tagFiles'
    d = {}
    gcTotal, atTotal, geneGCDict = get_gc_totals(annodb, dna)
    for i in range(1, pow(2, len(tagFiles))):
        print 'subset', i
        l = [tagFile for (j, tagFile) in enumerate(tagFiles)
             if i & pow(2, j)]
        snps = read_tag_files(l, *args)
        results = analyze_nonsyn(snps, annodb, al, dna,
                                 gcTotal, atTotal, geneGCDict)
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
        snps = read_vcf_singleton(poolFile)
        results = analyze_nonsyn(snps, annodb, al, dna,
                                 gcTotal, atTotal, geneGCDict)
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
            snps = read_tag_files(pfiles, tagfunc, replicatefunc,
                                     minRep, *args)
            results = analyze_nonsyn(snps, annodb, al, dna,
                                     gcTotal, atTotal, geneGCDict)
            k = reduce(lambda x,y:x + y, [get_pool_tags(t) for t in pfiles])
            d[k] = results[:nbest]
    return d



#######################################################################
# simple examples of how to run a complete analysis

def analyze_nonsyn(snps, annodb, al, dna, **kwargs):
    'scores genes based on non-synonymous SNPs only'
    gsd = map_snps_chrom1(snps, al, dna)
    gsd = filter_nonsyn(gsd)
    return score_genes_pooled(gsd, dnaseq=dna, annodb=annodb,
                              **kwargs)

def analyze_nonsyn_groups(groups, snps, annodb, al, dna,
                          **kwargs):
    'scores groups based on non-synonymous SNPs only'
    gsd = map_snps_chrom1(snps, al, dna)
    gsd = filter_nonsyn(gsd)
    return score_groups_pooled(gsd, groups, dnaseq=dna,
                               annodb=annodb, **kwargs)


def analyze_monodom(genome, vcfFiles=glob.glob('*.vcf'),
                    genesFile='knownGene.txt'):
    snps = read_tag_files(vcfFiles, filterFunc=None, add_attrs=None)
    annodb, al, exonGene, totalSize, geneLengths = \
            read_exon_annots(genome, genesFile)
    gsd = map_snps(snps, al, genome, exonGene)
    return score_genes(gsd, len(vcfFiles), totalSize, geneLengths)

def main():
    print 'reading gene annotations from', sys.argv[1]
    annodb, al, dna = read_genbank_annots(sys.argv[1])
    tagFiles = sys.argv[2:]
    print 'reading tag files:', tagFiles
    snps = read_tag_files(tagFiles)
    print 'scoring genes...'
    results = analyze_nonsyn(snps, annodb, al, dna)
    print 'top 20 hits:', results[:20]

if __name__ == '__main__':
    main()
    
