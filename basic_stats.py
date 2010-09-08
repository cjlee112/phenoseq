
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
    try:
        return max(realhits.values())
    except ValueError: # no hits at all
        return 0


def poisson_screen(nstrain, ntarget, pois):
    'screen based on poisson; good for high densities'
    realhits = numpy.core.zeros(ntarget)
    for j in xrange(nstrain):
        while True:
            hits = pois.rvs(ntarget)
            if hits.max() > 0:
                break
        realhits += hits
    return realhits.max()


def sample_maxhit(n, nstrain, ntarget, ngene, lambdaTrue,
                  lambdaError=0., pFail=0., nwait=1000, pPoisson=0.1):
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
    lam = lambdaTrue * (1. - pFail) + lambdaError
    pois = stats.poisson(nstrain * lam)
    if stats.poisson(ntarget * lam).pmf(0) > pPoisson: # use waiting times
        usePoisson = False
        lam = lambdaTrue + lambdaError
        mu = 1. / lam # mean waiting time parameter
        f = 1. - exp(-ntarget / mu)
        r = numpy.random.random_sample(n * nstrain) # random x from 0 - 1
        start = (-mu) * numpy.log(1. - f * r) # location of 1st target mutation
        startIter = iter(start)
        waitTime = waiting_times(nwait, mu)
    else: # use Poisson screen
        usePoisson = True
        poisTarget = stats.poisson(lam)
    ngood = nok = 0
    for i in xrange(n):
        if usePoisson:
            maxreal = poisson_screen(nstrain, ntarget, poisTarget)
        else:
            maxreal = expon_screen(nstrain, ntarget, startIter, waitTime,
                                   pFail)
        nontargets = pois.rvs(ngene - ntarget) # #hits per gene
        maxbad = nontargets.max()
        if maxreal > maxbad: # real target hit more than any non-target
            ngood += 1
        elif maxreal == maxbad: # real target tied with best non-target
            nok += 1
    return float(ngood) / n, float(ngood + nok) / n


def calc_cutoff_success(n, nstrain, ntarget, ngene, c=75, npool=4,
                        epsilon=0.01, cutoffs=range(1,26),
                        ntrue=50, ntotal=4e+06):
    '''compute probability that top hit(s) is a real target, as a
    function of cutoff threshold'''
    pmf = stats.binom(c, (1. - epsilon) / npool)
    pmf2 = stats.binom(c, epsilon)
    pFail = pmf.cdf(numpy.array(cutoffs) - 1) # cdf() offset by one
    lambdaError = (float(ntotal) / ngene) * pmf2.sf(numpy.array(cutoffs) - 1)
    return [sample_maxhit(n, nstrain, ntarget, ngene, float(ntrue) / ngene,
                          lambdaError[i], pFail[i])
            for i in range(len(cutoffs))]

