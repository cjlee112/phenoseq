
######################################
Phenotype Sequencing Simulation Module
######################################

.. module:: simulate
   :synopsis: Phenotype sequencing simulator
.. moduleauthor:: Christopher Lee <leec@chem.ucla.edu>
.. sectionauthor:: Christopher Lee <leec@chem.ucla.edu>

Background
----------

Phenotype sequencing seeks to identify some set of *target genes*
(i.e. the set of genes where mutations can cause a given phenotype).
A mathematical description of how we model phenotype sequencing
is presented in :doc:`/theory`.  Here we just outline a few key 
aspects of how to use the simulation module:

* It uses the Poisson distribution (i.e. uniform hit probability)
  to generate random samples of mutation counts in each gene in
  the genome.  Each mutant strain is required to have at least
  one mutation in the target region; i.e. at least one target
  gene must be mutated, to cause the phenotype.

* Two distinct *gene size models* are available: *uniform* gene size
  (all genes in the genome have the exact same size); *variable*
  gene size (you specify a fixed number of size classes that
  represent the gene size distribution).  Note that our simulations
  with the *E. coli* gene size distribution showed that using
  10 gene size classes gave nearly identical results to using
  the exact *E. coli* gene size distribution (i.e. one size class
  for each of the 4221 genes), and that using a uniform gene
  size model gave very similar, but slightly more conservative
  results.  See Fig. 6 in
  `the paper <http://www.plosone.org/article/info:doi/10.1371/journal.pone.0016517>`_.

  Of course, the uniform size model is faster than the variable
  size model, and a smaller number of size classes is faster than
  a large number of size classes.

* Two basic measures of success are available: the *probability*
  that the top scoring hit is a true target; and the *number* of true
  targets discovered on average at a given false discovery rate (FDR).

* The basic parameters that you must provide for a simulation experiment
  are:

  * **the mutagenesis density**, typically expressed in average number
    of hits per gene.

  * **the number of true target genes**

  * **the number of total genes in the genome**

  * **the number of mutant strains** analyzed in your experiment

  * **the number of simulated experiments** you want to run.  Phenoseq
    can run about 1000 simulations (of 80 *E. coli* strains each) in about
    0.1 second on my MacbookPro.

  * **the mutation false positive density**: if non-zero, this allows
    you to simulate the effect of false positive errors (i.e. if your
    experimental data included a significant number of reported 
    mutations that are actually just sequencing errors).

  * **the mutation false negative rate**: if non-zero, this allows
    you to simulate the effect of failing to detect some fraction
    of the real mutations present in each strain.


Examples
--------

Modeling "Ideal" Phenotype Sequencing Yields
............................................

The discovery yield of phenotype sequencing depends crucially on
the number of mutant strains sequenced.  To analyze this, it is
convenient to begin in an idealized setting where the effects 
of sequencing error and coverage are ignored, i.e. assuming
that every mutation is detected, and there are no false positives.

Example: compute the probability of finding the target gene as
the top scoring hit, out of 25,000 mouse exome genes with a mutation
density of 0.01 mutations / gene, using one to eight mutant genomes::

   pYield = [simulate.sample_maxhit(1000, i, 1, 25000, .01)[0] for i in range(1, 9)]


Example: compute the probability of finding the target gene among
the top five hits, out of 25,000 mouse exome genes with a mutation
density of 0.1 mutations / gene, using one to eight mutant genomes::

   pYield = [simulate.sample_maxhit_rank(1000, i, 1, 25000, .1, fdr=.8)[1] for i in range(1, 9)]


Modeling the effects of sequencing error and coverage
.....................................................

Example: compute the yield for a single target gene in a mouse exome
of 25,000 genes at mutation density 0.01 mutation / gene, with a
total exome size of 50 Mb, for 8 mutant genomes with a pooling factor
of 4, 100x coverage, and default (1%) sequencing error rate::

   yieldMax, cutoff = optimal_yield(8, 4, 100, nmut=250, ntarget=1, ngene=25000, ntotal=50e+06)


Convenience Functions
---------------------

.. function:: sample_maxhit(n, nstrain, ntarget, ngene, lambdaTrue, lambdaError=0., pFail=0., residual=1e-06, targetProbs=None, genFunc=generate_hits)

   Runs phenotype sequencing simulation to compute probability that the
   highest scoring gene will be a true target (i.e. a gene that causes
   the phenotype), using a *uniform* gene size model for non-target
   genes (You *can* assign non-uniform gene sizes to the target genes).

   *n*: the number of simulated experiments to run.

   *nstrain*: the number of mutant strains to generate in each experiment.

   *ntarget*: the number of target genes that cause the phenotype.

   *ngene*: number of total genes in the genome

   *lambdaTrue*: the Poisson lambda parameter representing mean number of
   true mutations per gene.  E.g. for 50 mutations in 4000 genes
   lambdaTrue = 50./4000

   *lambdaError*: Poisson parameter representing mean number of sequencing
   errors (false positive mutation calls) per gene.

   *pFail*: probability that a given true mutation will not be detected
   (false negative).

   *residual*: cutoff threshold for truncating the Poisson distribution.
   That is, it enumerates all values of the hit-count *k* whose p-value
   is above the specified *residual*.

   *targetProbs*, if not None, must be a list of probabilities
   that a given target region mutation will be in each of the target
   genes.  (This must sum to 1.0).  If None, the target probabilities
   are uniform (equal for each target gene).

   *genFunc*: the function to use for generating a sample of simulations.

   Returns tuple *(pBest, pTie)*, where *pBest* is the probability
   that a target gene scores better than *any* non-target gene,
   and *pTie* is the probability
   that a target gene scores at least as high as *any* non-target gene.

.. function:: sample_maxhit_rank(n, nstrain, ntarget, ngene, lambdaTrue, lambdaError=0., pFail=0., fdr=0.67, residual=1e-06, targetProbs=None, genFunc=generate_hits)

   Identical to :func:`sample_maxhit`, except that it returns a
   probability distribution *pYield* representing the probability
   that at least *i* true targets will be detected at the specified
   false discovery rate.  In other words, ``pYield[0]`` is the probability
   that no true targets are discovered in one experiment,
   and ``pYield[i]`` is the probability
   that ``i`` true targets are discovered in one experiment.

   *fdr* is the false discovery rate, i.e. the maximum fraction of
   non-target genes allowed before truncating the list of top-scoring
   genes.

.. function:: calc_yield_varsize(n, nstrain, ntarget, geneSizes, lambdaTrue, decoyBins, lambdaError=0., pFail=0., fdr=0.67, residual=1e-06, targetProbs=None, genFunc=generate_hits_varsize, **kwargs)

   Similar to :func:`sample_maxhit_rank`, except that it uses a
   variable gene size model for non-target genes.  You can choose
   between two possible *genFunc* functions for generating
   a sample of simulations:
   :func:`generate_hits_varsize` uses a fixed set of target sizes
   for all the simulated experiments.  In other words, *geneSizes* must
   be a list whose length matches *ntarget*.  Alternatively, 
   :func:`generate_hits_varsize2` draws a different random sample of
   target sizes for each of the simulated experiments.  In this 
   case, *geneSizes* should just be a list of the sizes of all possibly relevant
   genes (typically, all genes in the genome).

   *kwargs*, if any, are passed as keyword arguments to *genFunc*.

.. function:: optimal_yield(nstrain, npool=4, c=75, epsilon=0.01, n=1000, nmut=50, ntarget=20, ngene=4200, ntotal=4.64e+06, threshold=0.98, nwait=4, mutFac=3., **kwargs)
   
   Model the achievable yield for a given phenotype sequencing experiment
   design, consisting of the following factors:

   * *nstrain*: the number of mutant genomes
   * *npool*: the pooling factor, i.e. number of strains per library
   * *c*: the total coverage used for all strains, i.e. average number
     of reads covering any given nucleotide position.  Thus, the average
     number of reads per strain at any given position is ``c/npool``.
   * *epsilon*: sequencing error rate per base per read.
   * *nmut*: total number of mutations per genome
   * *ntarget*: number of target genes
   * *ngene*: total number of genes
   * *ntotal*: total effective size of the genome for your experiment.
     If you are sequencing the whole genome it's simply the total nucleotides
     in your genome.  If you're sequencing an exome, it's the size of
     the exome you are sequencing.
   * *kwargs*, if any, are passed to :func:`simulate.sample_maxhit_rank()`.

   Returns a tuple ``yieldMax, cutoff`` where ``yieldMax`` is the best
   achievable yield from :func:`simulate.sample_maxhit_rank()`, and
   ``cutoff`` is the number of reads required to report a mutation
   at a given nucleotide.


