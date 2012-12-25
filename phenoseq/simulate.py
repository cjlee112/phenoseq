from optparse import OptionParser
import random
import numpy
from scipy import stats
from math import exp, log, ceil

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
    def __init__(self, lam, residual=1e-06, nmin=0, ngene=1, calcSF=False):
        self.ngene = ngene
        self.poisson = pois = stats.poisson(lam)
        nmax = nmin
        if nmin > 0:
            pTotal = pois.sf(nmin - 1) # sum prob for [nmin,infty)
        else:
            pTotal = 1.
        while pois.sf(nmax) / pTotal >= residual: # find max plausible hits
            if (pois.cdf(nmax) - 1. + pTotal) / pTotal < residual:
                nmin = nmax
            nmax += 1
        self.nmin = nmin # sample values will fall within this range
        self.nmax = nmax + 1
        self.pmf = pois.pmf(numpy.arange(nmin, nmax + 1))
        if calcSF: # save precomputed survival probabilities
            self.sf = pois.sf(numpy.arange(nmin - 1, nmax))

    def renormalize(self):
        'renormalize self.pmf to sum to 100%'
        self.pmf /= self.pmf.sum()

    def sample_totals(self, nsample, n):
        'get n draws of total counts in sample of size nsample'
        counts = numpy.random.multinomial(nsample, self.pmf, size=n)
        totals = counts * numpy.arange(self.nmin, self.nmax)
        return totals.sum(axis=1)


def gen_hits_base(n, nstrain, ntarget, ngene, prDecoy, lambdaTrue,
                  lambdaError=0., pFail=0., residual=0.001,
                  targetProbs=None):
    'slow test version of generate_hits()'
    assert prDecoy.nmin == 0
    i = n
    targetPois = stats.poisson(lambdaTrue)
    decoyPois = stats.poisson(lambdaTrue * nstrain)
    for j in xrange(n):
        targetCounts = numpy.zeros(ntarget)
        for k in xrange(nstrain):
            while True:
                i += 1
                if i >= n:
                    targetHits = targetPois.rvs((n, ntarget))
                    i = 0
                if targetHits[i].max() > 0: # select strain that hits target
                    break
            targetCounts += targetHits[i]
        decoyCounts = decoyPois.rvs(ngene - ntarget)
        kmax = max(targetCounts.max(), decoyCounts.max())
        targetRes = [((targetCounts == k).sum()) for k in range(kmax + 1)]
        decoyRes = [((decoyCounts == k).sum()) for k in range(kmax + 1)]
        yield numpy.array(targetRes), numpy.array(decoyRes)


def generate_hits(n, nstrain, ntarget, ngene, prDecoy, lambdaTrue,
                  lambdaError=0., pFail=0., residual=1e-06,
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
        poisErr = stats.poisson(nstrain * lambdaError)
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


def generate_hits_varsize(n, nstrain, ntarget, targetSizes, decoyBins,
                          lambdaTrue, lambdaError=0., pFail=0.,
                          residual=1e-06):
    '''generate lists of hit counts in real targets vs. nontargets in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.
    decoyCounts vectors are offset by poissonRange.nmin'''
    lam = lambdaTrue * (1. - pFail) + lambdaError
    total = sum(targetSizes)
    targetScores = [stats.poisson(s * nstrain * lam).
                    sf(numpy.arange(-1, int(nstrain * max(1, 2 * s * lambdaTrue))))
                    for s in targetSizes]
    prTarget = PoissonRange(total * lambdaTrue, residual / nstrain, 1)
    prTarget.renormalize() # condition on at least one hit in target region
    targetN = prTarget.sample_totals(nstrain, n)
    targetProbs = [(s / total * (1. - pFail)) for s in targetSizes] \
                  + [pFail]
    targetProbs = numpy.array(targetProbs)
    if lambdaError:
        poisErrs = [stats.poisson(s * nstrain * lambdaError)
                    for s in targetSizes]
        poisCounts = [pois.rvs((n, 1)) for pois in poisErrs]
        targetNoise = numpy.concatenate(poisCounts, axis=1)
    decoyCounts = [numpy.random.multinomial(pr.ngene, pr.pmf, size=n)
                   for pr in decoyBins]
    for i in xrange(n): # analyze replicates
        targetHits = numpy.random.multinomial(targetN[i], targetProbs)
        if lambdaError:
            hits = targetHits[:ntarget] + targetNoise[i]
        else:
            hits = targetHits[:ntarget]
        scores = [(targetScores[j][k], 1, 0)
                  for j,k in enumerate(hits) if k > 0]
        for j, pr in enumerate(decoyBins):
            scores += [(pr.sf[k], 0, g)
                       for k,g in enumerate(decoyCounts[j][i])
                       if g > 0 and k + pr.nmin > 0]
        scores.sort()
        yield scores


def generate_hits_varsize2(n, nstrain, ntarget, geneSizes, decoyBins,
                           lambdaTrue, lambdaError=0., pFail=0.,
                           residual=1e-06):
    '''this variant draws a set of different target genes for each of the
    n samples.'''
    lam = lambdaTrue * (1. - pFail) + lambdaError
    decoyCounts = [numpy.random.multinomial(pr.ngene, pr.pmf, size=n)
                   for pr in decoyBins]
    poisMax = stats.poisson(sum(geneSizes[-ntarget:]) * lambdaTrue)
    kMax = poisMax.isf(residual / nstrain) # biggest plausible k for target
    for i in xrange(n): # analyze replicates
        targetSizes = random.sample(geneSizes, ntarget) # draw target genes
        total = sum(targetSizes)
        poisTarget = stats.poisson(total * lambdaTrue)
        cvec = numpy.arange(kMax)
        pvec = poisTarget.pmf(cvec)
        counts = numpy.random.multinomial(nstrain, pvec[1:] / (1. - pvec[0]))
        targetN = (counts * cvec[1:]).sum() # total hits in target region
        targetProbs = [(s / total * (1. - pFail)) for s in targetSizes] \
                      + [pFail]
        targetProbs = numpy.array(targetProbs)
        targetHits = numpy.random.multinomial(targetN, targetProbs)
        if lambdaError:
            poisErrs = [stats.poisson(s * nstrain * lambdaError).rvs(1)[0]
                        for s in targetSizes]
            hits = targetHits[:ntarget] + numpy.array(poisErrs)
        else:
            hits = targetHits[:ntarget]
        scores = [(stats.poisson(targetSizes[j] * nstrain * lam).sf(k - 1),
                   1, 0) for j,k in enumerate(hits) if k > 0]
        for j, pr in enumerate(decoyBins):
            scores += [(pr.sf[k], 0, g)
                       for k,g in enumerate(decoyCounts[j][i])
                       if g > 0 and k + pr.nmin > 0]
        scores.sort()
        yield scores


def sample_maxhit(n, nstrain, ntarget, ngene, lambdaTrue, lambdaError=0.,
                  pFail=0., residual=1e-06, targetProbs=None,
                  genFunc=generate_hits):
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
    for realhits, nontargets in genFunc(n, nstrain, ntarget, ngene,
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
                       lambdaError=0., pFail=0., fdr=0.67, residual=1e-06,
                       targetProbs=None, genFunc=generate_hits):
    '''Compute distribution of real targets detected above fdr, in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.'''
    pYield = numpy.zeros(ntarget + 1)
    p = 1. / n
    lam = lambdaTrue * (1. - pFail) + lambdaError
    prDecoy = PoissonRange(nstrain * lam, residual / ngene)
    for realhits, nontargets in genFunc(n, nstrain, ntarget, ngene,
                        prDecoy, lambdaTrue, lambdaError, pFail, residual,
                        targetProbs):
        nhit = lastReal = 0
        for i in range(len(realhits) - 1, 0, -1):
            nhit += realhits[i]
            nbad = nontargets[max(i - prDecoy.nmin, 0):].sum()
            if float(nbad) / (nbad + nhit) > fdr:
                break
            lastReal = nhit
        pYield[lastReal] += p # update probability dist'n
    return pYield

def calc_yield_varsize(n, nstrain, ntarget, geneSizes, lambdaTrue, decoyBins,
                       lambdaError=0., pFail=0., fdr=0.67, residual=1e-06,
                       targetProbs=None, genFunc=generate_hits_varsize,
                       **kwargs):
    '''Compute distribution of real targets detected above fdr, in
    multihit screen, using probability of hitting target
    gene(s) conditioned on requiring at least one such hit.'''
    pYield = numpy.zeros(ntarget + 1)
    p = 1. / n
    for scores in genFunc(n, nstrain, ntarget, geneSizes, decoyBins,
                          lambdaTrue, lambdaError, pFail, residual,
                          **kwargs):
        nhit = nbad = 0
        for pval, targetGenes, decoyGenes in scores:
            if decoyGenes > 0:
                nbad += decoyGenes
                if float(nbad) / (nbad + nhit) > fdr:
                    break
            else:
                nhit += targetGenes
        pYield[nhit] += p # update probability dist'n
    return pYield


def build_decoy_bins(geneSizes, nstrain, lam, nbin=10, **kwargs):
    'construct nbin PoissonRange objects for variable gene sizes'
    geneSizes.sort()
    n = float(len(geneSizes)) / nbin # number of genes per bin, float
    l = []
    start = 0
    for i in range(1, nbin + 1):
        end = int(i * n)
        avgSize = sum(geneSizes[start:end]) / float(end - start)
        l.append(PoissonRange(avgSize * nstrain * lam, calcSF=True,
                              ngene=end - start, **kwargs))
        start = end # start of the next interval
    return l


def expectation(a):
    'compute expectation value from an array of probabilities'
    return (a * numpy.arange(len(a))).sum()

class CutoffSuccess(object):
    def __init__(self, ngene, c=75, npool=4, epsilon=0.01,
                 cutoffs=range(1,26), ntotal=4e+06, mutFac=3.):
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

def find_cutoff(l, r, func):
    'find v where func(v) switches from False to True'
    while r is None or r - l > 1:
        if r is None:
            m = l * 2
        else:
            m = (l + r) / 2
        if func(m):
            r = m
        else:
            l = m
    return l


def get_yield_params(npool=4, c=75, epsilon=0.01, nmut=50, ngene=4200,
                     ntotal=4.64e+06, mutFac=3, **kwargs):
    pmfMut = stats.binom(c, (1. - mutFac * epsilon) / npool)
    pmfErr = stats.binom(c, epsilon)
    def no_false_positives(cut):
        return pmfErr.sf(cut) * ntotal < 0.05
    i = find_cutoff(0, c, no_false_positives) + 1
    def too_many_false_negatives(cut):
        return pmfMut.cdf(cut) * nmut > 0.05
    j = find_cutoff(0, c, too_many_false_negatives)
    return pmfMut, pmfErr, i, j


def optimal_yield(nstrain, npool=4, c=75, epsilon=0.01, n=1000, nmut=50,
                  ntarget=20, ngene=4200, ntotal=4.64e+06, **kwargs):
    'find optimal yield by scanning mutation call cutoff values'
    pmfMut, pmfErr, i, j = \
         get_yield_params(npool, c, epsilon, nmut, ngene, ntotal)
    lam = float(nmut) / ngene
    genesize = ntotal / ngene
    if i < j: # well-separated cutoffs, easy to discriminate
        i = (i + j) / 2 # just point between the two cutoffs
        d = sample_maxhit_rank(n, nstrain, ntarget, ngene, lam,
                               pmfErr.sf(i) * genesize, pmfMut.cdf(i),
                               **kwargs)
        return expectation(d), i + 1
    else:
        l = []
        for i in range(j, i + 1): # search usable range
            pFail = pmfMut.cdf(i)
            d = sample_maxhit_rank(n, nstrain, ntarget, ngene, lam,
                                   pmfErr.sf(i) * genesize, pFail,
                                   **kwargs)
            l.append((expectation(d), i + 1))
            if pFail > 0.95: # hopeless past this point
                break
        l.sort()
        return l[-1]


def optimal_npool(nstrain, minFrac=0.98, l=1, r=None, totalCov=323,
                  *args, **kwargs):
    'find maximum npool whose yield >= minFrac * yieldMax'
    yieldMax, i = optimal_yield(nstrain, l, *args, **kwargs) # baseline
    while r is None or r - l > 1:
        if r is None: # expand upper bound
            m = 2 * l
        else: # midpoint of interval
            m = (l + r) / 2
        nlib = int(ceil(float(nstrain) / m))
        c = int(totalCov / nlib)
        y, i = optimal_yield(nstrain, m, c, *args, **kwargs)
        if y >= yieldMax * minFrac: # still acceptable
            l = m
            yieldLast = y
        else: # beyond acceptable range
            r = m
    return yieldLast, l


def get_linear_devs(data):
    x0, y0 = data[0]
    w = data[-1][0] - x0
    h = data[-1][1] - y0
    a = h / w
    model = [(y0 + a * (t[0] - x0)) for t in data[1:-1]]
    return [(t[1] - model[i]) / abs(model[i])
            for i,t in enumerate(data[1:-1])]
    

def min_coverage(nstrain, npool, minFrac=0.98, l=10,
                 minSep=5, *args, **kwargs):
    'find minimum coverage whose yield >= minFrac * yieldMax'
    def well_separated(cov):
        pmfMut, pmfErr, i, j = \
         get_yield_params(npool, cov, **kwargs)
        return i + minSep < j
    r = find_cutoff(l, None, well_separated) + 1 # cov w/ good sig-noise sep

    yieldMax, i = optimal_yield(nstrain, npool, r, *args, **kwargs) # baseline
    holder = [None]
    def acceptable_yield(c):
        y, i = optimal_yield(nstrain, npool, c, *args, **kwargs)
        if y >= yieldMax * minFrac:
            holder[0] = y # record the yield
            return True
    i = find_cutoff(l, r, acceptable_yield) + 1
    return holder[0], i
    

def min_cost(nstrain=32, laneCov=4300, laneCost=800., libCost=50.,
             ntarget=5, **kwargs):
    'find minimum cost protocol for analyzing nstrain mutant strains'
    l = []
    for nlib in range(1, nstrain / 2 + 1):
        if nstrain % nlib == 0:
            npool = nstrain / nlib
            hits, cov = min_coverage(nstrain, npool, ntarget=ntarget,
                                     **kwargs)
            cost = nlib * libCost + float(laneCost) * nlib * cov / laneCov
            print nlib, cov, cost
            l.append((cost, nlib, cov, hits))
    l.sort()
    return l[0]

def main():
    usage = "usage: %prog nstrain laneCov laneCost libCost ntarget"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    print '''Strains: %s
Coverage per lane: %s
Cost per lane: $%s
Cost per library: $%s
Target genes: %s
''' % tuple(args)
    t = min_cost(*[int(s) for s in args)
    print '''Minimum experiment cost: $%0.2f
Total libraries: %d
Coverage per library: %1.0f
Expected hits: %1.2f
''' % t

if __name__ == '__main__':
    main()
