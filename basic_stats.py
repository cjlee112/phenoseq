
import random
import numpy
from scipy import stats
from math import exp, log

def sample_genes(ngene=4200, nmut=100):
    """Generate counts of genes with 1, 2 or more mutations,
    treating each gene as equally likely."""
    d = {}
    for i in range(nmut):
        k = random.randint(0, ngene)
        d[k] = d.get(k, 0) + 1
    c = {}
    for i in d.values():
        c[i] = c.get(i, 0) + 1
    return c

def waiting_times(n, mu):
    'generates infinite sequence of exponential waiting times'
    pdf = stats.expon(0, mu)
    i = 0
    data = pdf.rvs(n)
    while True:
        try:
            yield data[i]
        except IndexError:
            i = 0
            data = pdf.rvs(n)
            yield data[i]
        i += 1

def sample_maxhit(n, nstrain, ntarget, ngene, lam):
    '''Compute fraction of best hits that are real targets in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.
    n: number of replicates
    nstrain: number of strains to use for each replicate
    ntarget: number of target genes that cause the phenotype
    ngene: number of total genes in the genome
    lam: the Poisson lambda parameter representing mean number of
    mutations per gene.  E.g. for 50 mutations in 4000 genes
    lam = 50./4000'''
    mu = 1. / lam # mean waiting time parameter
    f = 1. - exp(-ntarget / mu)
    r = numpy.random.random_sample(n * nstrain) # random x from 0 - 1
    start = (-mu) * numpy.log(1. - f * r) # location of 1st target mutation
    waitTime = waiting_times(2. * n * nstrain * ntarget * lam, 1. / lam)
    pois = stats.poisson(nstrain * lam)
    istart = 0
    ngood = nok = 0
    for i in xrange(n):
        realhits = {}
        for j in xrange(nstrain):
            x = start[istart]
            istart += 1
            while x < ntarget: # record target mutations in this strain
                k = int(x)
                realhits[k] = realhits.get(k, 0) + 1
                x += waitTime.next()
        nontargets = pois.rvs(ngene - ntarget) # #hits per gene
        maxreal = max(realhits.values())
        maxbad = nontargets.max()
        if maxreal > maxbad: # real target hit more than any non-target
            ngood += 1
        elif maxreal == maxbad: # real target tied with best non-target
            nok += 1
    return float(ngood) / n, float(ngood + nok) / n

