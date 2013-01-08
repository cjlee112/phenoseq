
===========================
Phenotype Sequencing Theory
===========================

Mutation Count Distributions in Target and Non-Target Genes
-----------------------------------------------------------

To model independent phenotype selection events, we need probability
distributions for the number of mutations in *target genes* (genes
where mutations can cause the desired phenotype) and for non-target genes.
First we consider a simple model in which genes are assumed to
have uniform size, and then extend it to variable gene sizes.
Note that we treat "size" as a general parameter combining the many
factors that affect the probability of observing a mutation in 
a given region, including not only its length in the genomic sequence,
but all other factors such as its base composition, mutational biases,
and selection biases.

Independent mutations occurring over a genome are commonly modeled
using the Poisson distribution.  Specifically, if the expectation value
for the number of mutations expected in a region is :math:`\lambda`, 
the probability of observing exactly :math:`k` mutations in that region is
given by

.. math:: p(k|\lambda) = \frac{e^{-\lambda}\lambda^k}{k!}

Now consider the following simple model of a phenotype selection screen.
Assume that mutations at a subset of sites in the genome can cause
the desired phenotype; call this the "target region" and
designate its total size as :math:`\tau`.  Defining the density of
mutations resulting from mutagenesis as :math:`\mu`, the expected
number of mutations in the target region is :math:`\lambda=\mu\tau`.
For convenience we express :math:`\mu` in terms of the number of 
mutations per gene, and :math:`\tau` as simply the number of target genes.
For clarity we will designate the number of hits in the whole target region
as :math:`h` (to distinguish it from the number of hits in a 
single gene :math:`k`).
To model the effect of the phenotype selection screen, we require
that at least one mutation be present in the target region, 
which alters the conditional mutation probability:

.. math:: p(h|h \ge 1, \mu,\tau) = \frac{p(h|\mu,\tau)}{1-p(h=0|\mu,\tau)}
          = \frac{e^{-\mu\tau}(\mu\tau)^h}{(1-e^{-\mu\tau})h!}
          = \frac{(\mu\tau)^h}{(e^{\mu\tau}-1)h!}

Now consider the distribution of counts in a single target gene.
We simply calculate the likelihood of individual counts from the whole
target region to be assigned to a single gene, treating this as
a binomial event with probability :math:`\theta=1/\tau`.

.. math:: p(k|\mu,\tau) = \sum_{h \ge k}{p(k|h)p(h|h \ge 1, \mu, \tau)}
          = \sum_{h=k}^{\infty}{{h \choose k}\tau^{-k}(1-\tau^{-1})^{h-k}
          \frac{(\mu\tau)^h}{(e^{\mu\tau}-1)h!}}

.. math:: =\sum_{h=k}^{\infty}{\frac{h!}{(h-k)!k!}\frac{\mu^k}{e^{\mu\tau}-1}
          \frac{\mu^{h-k}\tau^{h-k}(1-\tau^{-1})^{h-k}}{h!}}

.. math:: = \frac{\mu^k}{(e^{\mu\tau}-1)k!}
          \sum_{i=0}^{\infty}{\frac{(\mu\tau(1-\tau^{-1}))^i}{i!}}
          = \frac{e^{\mu(\tau-1)}\mu^k}{(e^{\mu\tau}-1)k!}

for :math:`k \ge 1`.
We have to treat :math:`k=0` as a special case, since :math:`p(h=0)=0`:

.. math:: p(k=0|\mu,\tau) = \frac{e^{\mu(\tau-1)} -1}{e^{\mu\tau}-1}

Note that while this is the correct marginal distribution for a single
target gene, it should *not* be used to simulate real phenotype
sequencing experiments.  
Simply applying this marginal distribution independently
to each target gene is not strictly valid, since the :math:`k` values for
target genes in a phenotype sequencing experiment
are *not* independent (due to the :math:`h \ge 1`
condition for the entire target region)!

By contrast, for a non-target gene, the count distribution is just:

.. math:: p(k|\mu) = \frac{e^{-\mu}\mu^k}{k!}

Simulating a Phenotype Sequencing Experiment
--------------------------------------------

Target Genes
............

For a set of :math:`s` independent mutant strains that pass
the phenotype screen, the distribution of the 
total number of mutations in the target region simply follows
the sum of :math:`s` independent draws from the conditional
distribution :math:`p(h|h \ge 1, \mu,\tau)`.  
We model this as follows: we extract the vector
of values :math:`p_h = p(h|h \ge 1, \mu,\tau)`
for a confidence interval :math:`h_{min},h_{max}` such that
:math:`\sum_{h=h_{min}}^{h_{max}}{p_h} \ge 1-\delta`
for a stringent confidence threshold :math:`\delta`, 
construct a multinomial
distribution from this probability vector, and draw samples
of :math:`s` counts each from this multinomial.  Specifically,
each draw is a vector of :math:`\{n_h\}` observation counts for
each possible outcome :math:`h`, such that :math:`\sum_h{n_h}=s`.
This yields a sample distribution for the total number of mutations 
:math:`m=\sum_h{hn_h}` observed in the target region.

Given :math:`m` mutations in the target region, we model the distribution
of mutation counts in individual target genes as follows.  Assuming
that there are :math:`\tau` total genes in the target region,
we construct a multinomial based on a probability vector of
uniform gene probabilities :math:`p_i = 1/\tau`, and draw a sample
of :math:`m` counts, i.e. a vector :math:`\{n_i\}` such that
:math:`\sum_i{n_i}=m`.  The :math:`n_i` represent the individual
mutation counts in each target gene.  
We then count the number of target genes :math:`g_k`
with a specified number of mutations :math:`k`.
We sample their distribution
by generating :math:`n=1000` replicates of the above process,
for any specific set of input parameters (:math:`\mu, \tau, s`, etc.).

Non-target Genes
................

We modeled the distribution of mutation counts in non-target genes by
a similar methodology.  If :math:`\mu` is the expected number of
mutations per non-target gene in a single mutant strain, the
distribution of total mutations per gene in :math:`s` independent
strains is itself just a Poisson with mean :math:`s\mu`:

.. math:: p(k|s, \mu) = \frac{e^{-s\mu}(s\mu)^k}{k!}

Again, we extract a probability vector of values
:math:`p_k = p(k|s,\mu)`
for a confidence interval :math:`k_{min},k_{max}` such that
:math:`\sum_{k=k_{min}}^{k_{max}}{p_k} \ge 1-\delta`,
and construct a multinomial distribution from this probability
vector.  The distribution of the number of genes :math:`g_k'` that 
contain exactly :math:`k` mutations is given by drawing :math:`g-\tau`
counts from this multinomial, where :math:`g` is the total number
of genes in the genome, and :math:`\tau` is the number of target
genes in the genome.


