
import random
import numpy
from scipy import stats
from math import exp, log

class HitCountList(object):
    'efficient container for sparse lists, i.e. where many entries are 0'
    def __init__(self, counts):
        self.max = counts.max()
        n = 0
        total = len(counts) # total number of genes
        d = {}
        for c in counts:
            if c not in d:
                i = (counts == c).sum()
                d[c] = i # save this hit count
                n += i
                if n >= total: # counted all genes
                    break
        self.d = d

    def __len__(self):
        return self.max + 1

    def __getitem__(self, k):
        try:
            return self.d[k]
        except KeyError:
            return 0

class PoissonRange(object):
    'represents confidence interval of poisson dist, above residual'
    def __init__(self, lam, residual=0.001, nmin=0):
        self.poisson = pois = stats.poisson(lam)
        nmax = nmin
        while pois.sf(nmax) >= residual: # find max plausible hits
            if pois.cdf(nmax) < residual:
                nmin = nmax
            nmax += 1
        self.nmin = nmin # sample values will fall within this range
        self.nmax = nmax
        self.pmf = pois.pmf(numpy.arange(nmin, nmax))

    def renormalize(self):
        'renormalize self.pmf to sum to 100%'
        self.pmf /= self.pmf.sum()

    def sample_totals(self, nsample, n):
        'get n draws of total counts in sample of size nsample'
        counts = numpy.random.multinomial(nsample, self.pmf, size=n)
        totals = counts * numpy.arange(self.nmin, self.nmax)
        return totals.sum(axis=1)


def sample_maxhit_gen(n, nstrain, ntarget, ngene, prDecoy, lambdaTrue,
                      lambdaError=0., pFail=0., residual=0.001,
                      targetProbs=None):
    '''generate lists of hit counts in real targets vs. nontargets in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.
    decoyCounts vectors are offset by poissonRange.nmin'''
    prTarget = PoissonRange(ntarget * lambdaTrue, residual / nstrain, 1)
    prTarget.renormalize() # condition on at least one hit in target region
    targetN = prTarget.sample_totals(nstrain, n)
    if targetProbs: # user-supplied target size vector
        targetProbs = [(p * (1. - pFail)) for p in targetProbs] + [pFail]
    else: # uniform target sizes
        targetProbs = ((1. - pFail) / ntarget,) * ntarget + (pFail,)
    targetProbs = numpy.array(targetProbs)
    if lambdaError:
        poisErr = stats.poisson(lambdaError)
        targetNoise = poisErr.rvs((n, ntarget))
    decoyCounts = numpy.random.multinomial(ngene - ntarget,
                                           prDecoy.pmf, size=n)
    for i in xrange(n): # analyze replicates
        targetHits = numpy.random.multinomial(targetN[i], targetProbs)
        if lambdaError:
            hits = targetHits[:ntarget] + targetNoise[i]
        else:
            hits = targetHits[:ntarget]
        nmax = hits.max()
        if nmax < ntarget: # list all counts from 0 to nmax
            targetCounts = [(hits == j).sum() for j in xrange(nmax + 1)]
        else: # use efficient container for sparse list
            targetCounts = HitCountList(hits)
        yield targetCounts, decoyCounts[i]


def sample_maxhit(n, nstrain, ntarget, ngene, lambdaTrue, lambdaError=0.,
                  pFail=0., residual=0.001, targetProbs=None):
    '''Compute fraction of best hits that are real targets in
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
    ngood = nok = 0
    lam = lambdaTrue * (1. - pFail) + lambdaError
    prDecoy = PoissonRange(nstrain * lam, residual / ngene)
    for realhits, nontargets in sample_maxhit_gen(n, nstrain, ntarget, ngene,
                        prDecoy, lambdaTrue, lambdaError, pFail, residual,
                        targetProbs):
        maxreal = len(realhits) - 1 - prDecoy.nmin
        if nontargets[max(maxreal + 1, 0):].sum() == 0: # no better nontarget
            if maxreal < len(nontargets) and nontargets[maxreal]:
                nok += 1 # real target tied with best non-target
            else: # real target hit more than any non-target
                ngood += 1
    return float(ngood) / n, float(ngood + nok) / n


def sample_maxhit_rank(n, nstrain, ntarget, ngene, lambdaTrue,
                       lambdaError=0., pFail=0., fdr=0.67, residual=0.001,
                       targetProbs=None):
    '''Compute distribution of real targets detected above fdr, in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.'''
    d = {}
    p = 1. / n
    lam = lambdaTrue * (1. - pFail) + lambdaError
    prDecoy = PoissonRange(nstrain * lam, residual / ngene)
    for realhits, nontargets in sample_maxhit_gen(n, nstrain, ntarget, ngene,
                        prDecoy, lambdaTrue, lambdaError, pFail, residual,
                        targetProbs):
        nhit = lastReal = 0
        for i in range(len(realhits) - 1, 0, -1):
            nhit += realhits[i]
            nbad = nontargets[max(i - prDecoy.nmin, 0):].sum()
            if float(nbad) / (nbad + nhit) > fdr:
                break
            lastReal = nhit
        d[lastReal] = d.get(lastReal, 0.) + p # update probability dist'n
    return d

def expectation(d):
    e = 0.
    for x, p in d.items():
        e += x * p
    return e

class CutoffSuccess(object):
    def __init__(self, ngene, c=75, npool=4, epsilon=0.01,
                 cutoffs=range(1,26), ntotal=4e+06, mutFac=1.):
        self.pmf = stats.binom(c, (1. - mutFac * epsilon) / npool)
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

    def calc_ranks(self, n, nstrain, ntarget, ntrue=50, fdr=0.67):
        '''compute expected number of real hit(s) based on rankCut, as a
        function of cutoff threshold'''
        l = []
        for i in range(len(self.pFail)):
            d = sample_maxhit_rank(n, nstrain, ntarget, self.ngene,
                                   float(ntrue) / self.ngene,
                                   self.lambdaError[i], self.pFail[i],
                                   fdr=fdr)
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

def optimal_nstrains(ntarget, ngene, nmut=50, n=1000, goal=0.8, r=None):
    'find number of strains that achieves desired yield'
    lam = float(nmut) / ngene
    l = ntarget
    while r is None or r - l > 1:
        if r is None: # expand upper bound
            m = 2 * l
        else: # midpoint of interval
            m = (l + r) / 2
        d = sample_maxhit_rank(n, m, ntarget, ngene, lam)
        v = expectation(d) # compute average yield
        if v > goal * ntarget: # exceeded goal
            r = m
        else: # below goal
            l = m
    return r

def optimal_yield(nstrain, npool=4, c=75, epsilon=0.01, n=1000, nmut=50,
                  ntarget=20, ngene=4200, ntotal=4.64e+06, threshold=0.98,
                  nwait=4, mutFac=1.):
    'find optimal yield by scanning mutation call cutoff values'
    lam = float(nmut) / ngene
    genesize = ntotal / ngene
    pmfMut = stats.binom(c, (1. - mutFac * epsilon) / npool)
    pmfErr = stats.binom(c, epsilon)
    for i in range(c):
        if pmfErr.sf(i) * ntotal < 2 * nmut: # noise starting to fall below signal
            break
    yieldMax = cutoff = None
    while yieldMax is None or yieldLast >= yieldMax * threshold:
        d = sample_maxhit_rank(n, nstrain, ntarget, ngene, lam,
                               pmfErr.sf(i) * genesize, pmfMut.cdf(i))
        yieldLast = expectation(d)
        if yieldMax is None or yieldLast > yieldMax:
            yieldMax = yieldLast
            cutoff = i + 1
        elif cutoff is not None and i - cutoff > nwait:
            break
        i += 1
    return yieldMax, cutoff


def optimal_npool(nstrain, minFrac=0.98, l=1, r=None, *args, **kwargs):
    'find maximum npool whose yield >= minFrac * yieldMax'
    yieldMax, i = optimal_yield(nstrain, l, *args, **kwargs) # baseline
    while r is None or r - l > 1:
        if r is None: # expand upper bound
            m = 2 * l
        else: # midpoint of interval
            m = (l + r) / 2
        y, i = optimal_yield(nstrain, npool=m, *args, **kwargs)
        if y >= yieldMax * minFrac: # still acceptable
            l = m
            yieldLast = y
        else: # beyond acceptable range
            r = m
    return yieldLast, l
