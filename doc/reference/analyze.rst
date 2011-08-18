
============================================
Phenoseq Data Analysis Classes and Functions
============================================

.. module:: analyze
   :synopsis: Phenotype sequencing experimental data analysis
.. moduleauthor:: Christopher Lee <leec@chem.ucla.edu>
.. sectionauthor:: Christopher Lee <leec@chem.ucla.edu>


Overview
--------

The **analyze** module provides functions for

* reading SNP datasets and gene annotation data;
* filtering SNPs to remove false positives or based on impact criteria
  (e.g. non-synonymous mutations only).
* mapping SNPs to genes in a variety of ways.
* scoring genes to identify the likely causes of the phenotype.

These functions use just a few standard "interfaces" to pass data:

* SNP datasets are simply passed as lists of :class:`SNP` objects, each with
  attributes describing its location, etc.

* gene:SNP mappings are passed as a dictionary whose keys are gene IDs,
  and whose associated values are the lists of SNPs in each gene.


A simple example of using this module to score a set of genes
as potential causes of a phenotype::

    from phenoseq import analyze
    print 'reading gene annotations from', sys.argv[1]
    annodb, al, dna = analyze.read_genbank_annots(sys.argv[1])
    tagFiles = sys.argv[2:]
    print 'reading tag files:', tagFiles
    snps = analyze.read_tag_files(tagFiles)
    print 'scoring genes...'
    results = analyze.analyze_nonsyn(snps, annodb, al, dna)
    print 'top 20 hits:', results[:20]



Convenience Functions
---------------------

Data Reading Functions
......................

.. function:: read_vcf(path, vcffields=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'), **kwargs)

   Reads a VCF file from *path*
   and returns a list of :class:`SNP` objects, with
   attributes provided by the VCF fields.  The 
   `list of VCF data fields <http://samtools.sourceforge.net/mpileup.shtml>`_
   is available here.

   Optional *kwargs* are passed directly to the :class:`SNP` constructor. 

.. function:: read_tag_files(tagFiles, tagfunc=get_filestem, filterFunc=filter_snps, replicatefunc=None, add_attrs=add_snp_attrs, *args, **kwargs)

   Merges SNPs from multiple library files, returning a list of 
   :class:`SNP` objects each with a ``tag`` attribute representing the
   specific library it was observed in.  Each tag could correspond to 
   a single sample, or to a pool of multiple samples depending on your
   experimental design.

   *tagFiles*: a list of VCF files each representing a SNP dataset for
   one tag (library).

   *tagfunc*: a function that takes a tag file path and returns a tag ID.

   *filterFunc*: a generator function that takes a set of SNPs and
   generates the subset that pass some desired criteria.  Can be 
   set to ``None``.

   *replicatefunc*: a function that takes a tag ID and returns a list
   of VCF file paths representing the replicate datasets for that
   tag.  If you specify this, you should also set
   ``filterFunc=analyze.filter_snps_repfiles`` to make it filter
   the SNPs using the replicate datasets.

   *add_attrs*: a function that adds attributes to each
   :class:`SNP` object.  Can be set to ``None``.

   *args*: optional additional arguments to be passed to the
   ``filterFunc`` filter function for filtering SNPs.

   *kwargs*: optional arguments to be passed to the :func:`read_vcf`
   function for reading each library file.

.. function:: read_genbank_annots(gbfile, fastafile=None, iseq=0, featureType='CDS')

   Constructs a (gene) annotation database from a Genbank genome file.
   NB: this assumes each gene consists of **one** interval.
   This cannot be used for multi-exon genes!
   It takes the following arguments:

   *gbfile* must be a path to a Genbank format file containing
   gene feature records.

   *fastafile*, if not None, must be a path to a FASTA format file
   containing the genome sequence.  If None, a path will be constructed
   automatically by replacing the *gbfile* suffix with ``.fna``.

   *iseq* must be the index of the genome sequence within both *gbfile*
   and *fastafile*, i.e. ``iseq=0`` means the first sequence in the file.

   *featureType* specifies the Genbank feature type to extract for
   constructing annotations.  By default, it extracts the coding sequence
   regions.

   Returns three values: 

   * an annotation database

   * an alignment container storing the alignment of the annotations to
     the genome sequence.

   * the genome sequence object.

   **note**: this function requires both the BioPython ``SeqIO`` module,
   and the Pygr ``seqdb`` and ``annotation`` modules.

.. function:: read_exon_annots(genome, genesFile='knownGene.txt')

   Constructs a (gene) annotation database from a UCSC knownGene
   transcript set.  It builds an exon annotation database from this
   transcript set, and merges transcripts that share a common exon
   into a single gene.  This mapping is suitable for use with 
   multi-exon genes.

   *genome* must be dictionary-like object that maps sequence IDs to
   sequence objects.

   *genesFile*: the path to the UCSC knownGene text file.  The coordinates
   in this file must be valid for use on the specified *genome*.

   Returns:

   * an annotation database

   * an alignment container storing the alignment of the annotations to
     the genome sequence.

   * a dictionary mapping exon annotation ID values to gene ID values.

   * the total size of the annotated exome.

   * a dictionary mapping gene ID values to their associated maximum
     transcript length.  This gives the effective target size of each
     gene in the exome.

   **note**: this function requires the Pygr ``seqdb`` and ``annotation`` modules.

SNP false-positive filtering functions
......................................

These take a list of SNPs as input, and generate a subset of the SNPs
that pass some specified criteria.

.. function:: filter_snps(snps, filterExpr='snp.af1 <= 0.5 and getattr(snp, "pv4", (0.,))[0] >= 0.01 and snp.qual > 90')

   *snps*: a list of :class:`SNP` objects.

   *filterExpr*: a valid Python expression to be evaluated for each SNP.
   The SNP passes the filter if the expression is True.

.. function:: filter_snps_repfiles(snps, replicateFiles, filterExpr='snp.af1 <= 0.5 and len(get_replicates(snp)) >= 2')

   *snps*: a list of :class:`SNP` objects.

   *repFiles*: a list of VCF files representing replicate lanes to be
   used for testing whether the SNP replicated adequately.

   *filterExpr*: a valid Python expression to be evaluated for each SNP.
   The SNP passes the filter if the expression is True.


SNP - Gene mapping functions
............................

.. function:: map_snps(snps, al, genome, exonGene)

   Maps SNPs to genes using the specified alignment to exon annotations,
   and the exon to gene mapping given by *exonGene*.  This function is
   suitable for multi-exon gene data read by :func:`read_exon_annots`.
   Returns a gene:snps dictionary.

   *snps*: a list of :class:`SNP` objects.

   *al*: an alignment of genome sequence intervals to exon annotations.

   *genome*: the genome on which the exon annotations are aligned.

   *exonGene*: a dictionary whose keys are exon annotation IDs, 
   and whose associated values are gene IDs.

.. function:: map_snps_chrom1(snps, al, dna)

   Maps SNPs to genes on a single chromosome.  This function is suitable
   for single-interval gene data read by :func:`read_genbank_annots`.
   Returns a gene:snps dictionary.

   *snps*: a list of :class:`SNP` objects.

   *al*: an alignment of genome sequence intervals to exon annotations.

   *dna*: the chromosome sequence object on which the exon annotations are aligned.


SNP Impact filtering functions
..............................

.. function:: filter_nonsyn(geneSNPdict)

   Filters gene SNPs to just the subset that are non-synonymous.
   *geneSNPdict* should be a gene:snp dict from :func:`map_snps_chrom1`.
   Returns a filtered gene:snp dictionary.


Gene Scoring functions
......................

.. function:: score_genes_pooled(geneSNPdict, gcTotal=None, atTotal=None, geneGCDict=None, useBonferroni=True, dnaseq=None, annodb=None)

   Scores genes based on Poisson p-value for the total number of hits
   in each gene.  This is suitable for data where multiple samples were
   pooled in each library (making it impossible to tell which SNPs occurred
   in exactly which sample).

   *geneSNPdict* should be a gene:snp dict.

   *gcTotal,atTotal,geneGCDict* should be GC/AT base counts generated by
   :func:`get_gc_totals`.  These are optional: the function will compute
   them itself if they are not provided.

   *gcTotal*: if not None, the count of G and C bases in the genome
   sequence.  If None, it and *atTotal* are automatically calculated from the
   *dna* sequence passed to the constructor.

   *atTotal*: the count of A and T bases in the genome sequence.

   *geneGCDict*: if not None, must be a dictionary whose keys are
   geneIDs, and whose associated values are tuples giving
   ``(gcTotal, atTotal)`` specifically for each gene.  If None,
   the values are calculated automatically from the *annodb*
   passed to the constructor.

   *useBonferroni=True* forces the function to multiply the p-values
   by the total number of genes being tested.

   *dnaseq* is used to compute GC/AT base count totals if not provided 
   (see above).

   *annodb* is used to compute GC/AT base counts for each gene if not
   provided (see above).

   Returns a sorted list of ``(p_value, geneID)`` giving the 
   phenotype sequencing scores for all genes in which SNPs were reported.
   Since *smallest* p-values are the most significant, the top hits
   are at the *beginning* of this list.  The score calculations
   are described in :doc:`/theory`.

.. function:: score_genes(geneSNPdict, nsample, totalSize, geneLengths, useBonferroni=True)

   Scores genes based on the number of samples in which each gene was
   hit.  This is suitable for unpooled data, where each library contains
   only one sample.

   *geneSNPdict* should be a gene:snp dict.

   *nsample* should be the total number of samples.

   *totalSize* should represent the total size (in nucleotides) of the
   exome being analyzed.

   *geneLengths* should be a dictionary whose keys are gene IDs, and
   whose values are the size (in nucleotides) of that gene's exonic region
   (its effective target size for this analysis).

   *useBonferroni=True* forces the function to multiply the p-values
   by the total number of genes being tested.

   Returns a sorted list of ``(p_value, geneID)`` giving the 
   phenotype sequencing scores for all genes in which SNPs were reported.
   Since *smallest* p-values are the most significant, the top hits
   are at the *beginning* of this list.  The score calculations
   are described in :doc:`/theory`.


SNP Dataset Classes
-------------------

The SNP object
..............

.. class:: SNP(colnames, fields, add_attrs=None, **kwargs)

   A generic object for representing a SNP, simply as a Python
   object with attributes.  It provides no methods.  It is typically
   initialized from VCF data and carries all the VCF information
   fields as its attributes; attribute names are lower-case.  E.g.
   to access the AF1 field, use the SNP object's *af1* attribute.

   * *colnames*: a list of attribute names to use for the data in *fields*

   * *fields*: a list of data to bind as attributes to this SNP object.

   * *add_attrs*: an optional function that will add extra attributes
     to this object, derived from the attributes provided by *fields*

   * *kwargs*: optional key=value pairs to bind as additional attributes
     for this SNP object.

