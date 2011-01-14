
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
To model the effect of the phenotype selection screen, we require
that at least one mutation be present in the target region, 
which alters the conditional mutation probability:

.. math:: p(k|k \ge 1, \mu,\tau) = \frac{p(k|\mu,\tau)}{1-p(k=0|\mu,\tau)}
          = \frac{e^{-\mu\tau}(\mu\tau)^k}{(1-e^{-\mu\tau})k!}
          = \frac{(\mu\tau)^k}{(e^{\mu\tau}-1)k!}

