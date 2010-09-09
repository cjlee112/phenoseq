
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

def expon_screen(nstrain, ntarget, start, waitTime, pFail):
    '''screen based on exponential waiting times; good for low densities,
    very slow for high densities'''
    realhits = {}
    for j in xrange(nstrain):
        x = start.next()
        while x < ntarget: # record target mutations in this strain
            if random.random() >= pFail: # mutation detected
                k = int(x)
                realhits[k] = realhits.get(k, 0) + 1
            x += waitTime.next()
    if len(realhits) > 0:
        return numpy.array(realhits.values())
    else: # no hits at all
        return numpy.array((0,))


def poisson_screen(nstrain, ntarget, pois):
    'screen based on poisson; good for high densities'
    realhits = numpy.core.zeros(ntarget)
    for j in xrange(nstrain):
        while True:
            hits = pois.rvs(ntarget)
            if hits.max() > 0:
                break
        realhits += hits
    return realhits

def fast_screen(ntarget, pois):
    'assumes every strain will pass screen, so model all strains as a group'
    return pois.rvs(ntarget)

def sample_maxhit_gen(n, nstrain, ntarget, ngene, lambdaTrue,
                      lambdaError=0., pFail=0., nwait=1000, pPoisson=0.1,
                      pFast=0.1):
    '''generate lists of hit counts in real targets vs. nontargets in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.
    n: number of replicates
    nstrain: number of strains to use for each replicate
    ntarget: number of target genes that cause the phenotype
    ngene: number of total genes in the genome
    lambdaTrue: the Poisson lambda parameter representing mean number of
    true mutations per gene.  E.g. for 50 mutations in 4000 genes
    lambdaTrue = 50./4000
    lambdaError: Poisson parameter representing mean number of sequencing
    errors (false mutation calls) per gene.
    pFail: probability that true mutation detection will fail'''
    lam = lambdaTrue * (1. - pFail) + lambdaError
    pois = stats.poisson(nstrain * lam)
    p = stats.poisson(ntarget * lam).pmf(0)
    usePoisson = useFast = False
    if p * nstrain < pFast: # use fast screen approximation
        useFast = True
    elif p < pPoisson: # use Poisson screen
        usePoisson = True
        poisTarget = stats.poisson(lam)
    else: # use waiting times
        lam = lambdaTrue + lambdaError
        mu = 1. / lam # mean waiting time parameter
        f = 1. - exp(-ntarget / mu)
        r = numpy.random.random_sample(n * nstrain) # random x from 0 - 1
        start = (-mu) * numpy.log(1. - f * r) # location of 1st target mutation
        startIter = iter(start)
        waitTime = waiting_times(nwait, mu)
    for i in xrange(n): # analyze replicates
        if useFast:
            realhits = fast_screen(nstrain, pois)
        elif usePoisson:
            realhits = poisson_screen(nstrain, ntarget, poisTarget)
        else:
            realhits = expon_screen(nstrain, ntarget, startIter, waitTime,
                                    pFail)
        nontargets = pois.rvs(ngene - ntarget) # #hits per gene
        yield realhits, nontargets


def sample_maxhit(n, nstrain, ntarget, ngene, lambdaTrue,
                  lambdaError=0., pFail=0., nwait=1000, pPoisson=0.1,
                  pFast=0.1):
    '''Compute fraction of best hits that are real targets in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.'''
    ngood = nok = 0
    for realhits, nontargets in sample_maxhit_gen(n, nstrain, ntarget,
                                                  ngene, lambdaTrue,
                                                  lambdaError, pFail,
                                                  nwait, pPoisson, pFast):
        maxreal = realhits.max()
        maxbad = nontargets.max()
        if maxreal > maxbad: # real target hit more than any non-target
            ngood += 1
        elif maxreal == maxbad: # real target tied with best non-target
            nok += 1
    return float(ngood) / n, float(ngood + nok) / n


def sample_maxhit_rank(n, nstrain, ntarget, ngene, lambdaTrue,
                       lambdaError=0., pFail=0., nwait=1000, pPoisson=0.1,
                       pFast=0.1, rankCut=0.67):
    '''Compute distribution of real targets detected above rankCut rank, in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.'''
    d = {}
    p = 1. / n
    for realhits, nontargets in sample_maxhit_gen(n, nstrain, ntarget,
                                                  ngene, lambdaTrue,
                                                  lambdaError, pFail,
                                                  nwait, pPoisson, pFast):
        realhits.sort() # sort top hits to end of array
        nontargets.sort()
        ireal = -1
        ibad = lastReal = 0
        realend = -len(realhits)
        end = -(ngene - ntarget)
        while ireal >= realend:
            nhits = realhits[ireal]
            while ireal - 1 >= realend and \
                  realhits[ireal - 1] == nhits: # count targets w/ nhits
                ireal -= 1
            while ibad - 1 >= end and nontargets[ibad - 1] >= nhits:
                ibad -= 1 # count nontargets with >= nhits
            if float(ibad) / (ibad + ireal) <= rankCut: # still acceptable
                lastReal = -ireal # save as latest count of real hits
            else: # over threshold, so exit
                break
            ireal -= 1 
        d[lastReal] = d.get(lastReal, 0.) + p # update probability dist'n
    return d

def expectation(d):
    e = 0.
    for x, p in d.items():
        e += x * p
    return e

class CutoffSuccess(object):
    def __init__(self, ngene, c=75, npool=4, epsilon=0.01,
                 cutoffs=range(1,26), ntotal=4e+06):
        self.pmf = stats.binom(c, (1. - epsilon) / npool)
        self.pmf2 = stats.binom(c, epsilon)
        self.pFail = self.pmf.cdf(numpy.array(cutoffs) - 1) # cdf() offset by one
        self. lambdaError = (float(ntotal) / ngene) \
                            * self.pmf2.sf(numpy.array(cutoffs) - 1)
        self.ngene = ngene

    def calc_success(self, n, nstrain, ntarget, ntrue=50):
        '''compute probability that top hit(s) is a real target, as a
        function of cutoff threshold'''
        l = []
        for i in range(len(self.pFail)):
            l.append(sample_maxhit(n, nstrain, ntarget, self.ngene,
                                   float(ntrue) / self.ngene,
                                   self.lambdaError[i], self.pFail[i]))
        return l

    def calc_ranks(self, n, nstrain, ntarget, ntrue=50, rankCut=0.67):
        '''compute expected number of real hit(s) based on rankCut, as a
        function of cutoff threshold'''
        l = []
        for i in range(len(self.pFail)):
            d = sample_maxhit_rank(n, nstrain, ntarget, self.ngene,
                                   float(ntrue) / self.ngene,
                                   self.lambdaError[i], self.pFail[i],
                                   rankCut=rankCut)
            l.append(expectation(d))
        return l

def calc_pooling_success(n, nstrain, ntarget, ngene, c=75,
                         pools=range(1, 6) + range(7, 20, 2),
                         epsilon=0.01, cutoffs=range(5,20),
                         ntrue=50, ntotal=4e+06):
    results = []
    for npool in pools:
        print 'npool=', npool
        c = CutoffSuccess(ngene, c, npool, epsilon, cutoffs, ntotal)
        l = c.calc_success(n, nstrain, ntarget, ntrue)
        results.append(max([t[0] for t in l]))
    return results
