
============================================
Phenoseq Data Analysis Classes and Functions
============================================

.. module:: analyze
   :synopsis: Phenotype sequencing experimental data analysis
.. moduleauthor:: Christopher Lee <leec@chem.ucla.edu>
.. sectionauthor:: Christopher Lee <leec@chem.ucla.edu>

Convenience Classes
-------------------

* :class:`TagDict` represents a total SNP dataset for an entire set
  of mutant strains, typically constructed from datasets for multiple
  tagged libraries.  Its :meth:`TagDict.get_snps` method provides
  a convenient way to iterate over all the :class:`SNP` objects in ascending
  positional order.

* :class:`GeneSNPDict` represents the mapping of a :class:`TagDict` 
  variant dataset onto a set of genes.  Its :meth:`GeneSNPDict.get_scores`
  method provides a convenient way to get phenotype sequencing scores
  for each gene.

A simple example of using this module to score a set of genes
as potential causes of a phenotype::

    from phenoseq import analyze
    print 'reading gene annotations from', sys.argv[1]
    annodb, al, dna = analyze.read_genome_annots(sys.argv[1])
    tagFiles = sys.argv[2:]
    print 'reading tag files:', tagFiles
    tagDict = analyze.read_tag_files(tagFiles)
    gsd = analyze.GeneSNPDict(tagDict, annodb, al, dna)
    print 'scoring genes...'
    results = gsd.get_scores()
    print 'top 20 hits:', results[:20]



Convenience Functions
---------------------

.. function:: read_vcf(path, **kwargs)

   Reads a VCF file from *path*
   and returns a list of :class:`SNP` objects, with
   attributes provided by the VCF fields.  The 
   `list of VCF data fields <http://samtools.sourceforge.net/mpileup.shtml>`_
   is available here.

   Optional *kwargs* are passed directly to the :class:`SNP` constructor. 

.. function:: read_tag_files(tagFiles, tagfunc=get_filestem, replicatefunc=get_replicate_files, *args)

   Constructs and returns a :class:`TagDict` for a set of tags,
   each with replicate datasets.

   *tagFiles*: a list of VCF files each representing a SNP dataset for
   one tag.

   *tagfunc*: a function that takes a tag file path and returns a tag ID.

   *replicatefunc*: a function that takes a tag ID and returns a list
   of VCF file paths representing the replicate datasets for that
   tag.

   *args*: optional additional arguments to be passed to the
   :class:`ReplicateSet` constructor for each tagID.

.. function:: read_vcf_singleton(vcfFile)

   Constructs and returns a :class:`TagDict` for a single VCF file.

.. function:: read_genome_annots(gbfile, fastafile=None, iseq=0, featureType='CDS')

   Constructs a (gene) annotation database from a Genbank genome file.
   Returns three values: 

   * an annotation database

   * an alignment container storing the alignment of the annotations to
     the genome sequence.

   * the genome sequence object.

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

   **note**: this function requires both the BioPython ``SeqIO`` module,
   and the Pygr ``seqdb`` and ``annotation`` modules.

SNP Dataset Classes
-------------------

The SNP object
..............

.. class:: SNP(colnames, fields, add_attrs=None, **kwargs)

   A generic object for representing a SNP, simply as a Python
   object with attributes.  It provides no methods.

   * *colnames*: a list of attribute names to use for the data in *fields*

   * *fields*: a list of data to bind as attributes to this SNP object.

   * *add_attrs*: an optional function that will add extra attributes
     to this object, derived from the attributes provided by *fields*

   * *kwargs*: optional key=value pairs to bind as additional attributes
     for this SNP object.

Comparing data from replicate runs
..................................

.. class:: ReplicateSet(mergedFile, replicateFiles, minRep=2, maxF=0.5)

   A convenience class for reading data for two or more replicate
   SNP datasets.

   *mergedFile*: a VCF file path containing merged SNP data from all the 
   replicate datasets.

   *replicateFiles*: a list of file paths representing each of the 
   replicate SNP datasets, each in VCF format.

   *minRep*: the minimum number of times a given :class:`SNP` must
   replicate (i.e. the number of replicate datasets that report it)
   in order for it to be reported by :class:`ReplicateSet`.

   *maxF*: the maximum reported frequency that a given :class:`SNP`
   can have in the merged dataset, in order for it to be reported by
   :class:`ReplicateSet`.

.. method:: ReplicateSet.__iter__()

   Iterate over all :class:`SNP` objects that pass the above criteria.

.. method:: ReplicateSet.__getitem__(snp)

   The *snp*: key must be a :class:`SNP` object returned by the iterator.
   Its associated value is a list of the corresponding :class:`SNP` objects
   found in the replicate datasets (or an empty tuple if none).

Union of SNP data from multiple tagged libraries
................................................

.. class:: TagDict(tagDict)

   *tagDict* must be a dictionary whose keys are tag IDs (typically
   the six nucleotide tag sequence itself), and whose associated
   values are lists of arguments for constructing a :class:`ReplicateSet`
   for that tag.

   The resulting :class:`TagDict` object is itself a dictionary, whose
   keys are tag IDs, and whose associated values are :class:`ReplicateSet`
   objects represent the SNP data for that specific tag.

.. method:: TagDict.get_snps()

   Returns a list of all filtered SNPs (from all tags), sorted in ascending
   positional order.  Specifically, it is a list of tuples, each consisting
   of three values: ``(snp_pos, tag, snp)``, where ``snp_pos`` is the 
   nucleotide position of the SNP, ``tag`` is the tagID of the dataset
   it was observed in, and ``snp`` is the :class:`SNP` object providing
   all the SNP's attributes.

   Note that if the same variant is observed in different tag datasets,
   it will be reported once for each of those tag datasets.

Gene to SNP mapping and scoring
...............................

.. class:: GeneSNPDict(tagDict, annodb, al, dna, count_syn=False)

   Stores a mapping of SNPs to genes, filters them for non-synonymous
   substitutions, and returns gene scores based on the density
   of SNP hits within each gene.  It acts as a dictionary
   whose keys are gene IDs, and whose associated values are
   lists of SNPs found in that gene.  Each reported SNP in a gene
   is represented by a tuple of the form ``(tagID, snp)``, where
   the ``tagID`` indicates the library in which the SNP was reported,
   and ``snp`` is a :class:`SNP` object describing the SNP.

   It takes the following arguments:

   *tagDict*: a :class:`TagDict` representing the set of all variants
   found in all mutant strains.

   *annodb*: an annotation database representing the set of genes
   in the genome, specifically their coding regions.  *annodb*
   must be a dictionary whose keys are gene IDs, and whose associated
   value is a sequence interval representing the coding region.
   Typically this is obtained from :func:`read_genome_annots`.

   *al*: an alignment of the genes to the genome.
   Typically this is obtained from :func:`read_genome_annots`.

   *dna*: the genome sequence on which the annotations (and alignment)
   are mapped.
   Typically this is obtained from :func:`read_genome_annots`.

   *count_syn*: if True, includes synonymous mutations in the 
   analysis.  By default they are excluded from the analysis.

.. method:: GeneSNPDict.get_scores(gcTotal=None, atTotal=None, geneGCDict=None)

   Returns a sorted list of ``(p_value, geneID)`` giving the 
   phenotype sequencing scores for all genes in which SNPs were reported.
   Since *smallest* p-values are the most significant, the top hits
   are at the *beginning* of this list.  The score calculations
   are described in :doc:`/theory`.

   *gcTotal*: if not None, the count of G and C bases in the genome
   sequence.  If None, it and *atTotal* are automatically calculated from the
   *dna* sequence passed to the constructor.

   *atTotal*: the count of A and T bases in the genome sequence.

   *geneGCDict*: if not None, must be a dictionary whose keys are
   geneIDs, and whose associated values are tuples giving
   ``(gcTotal, atTotal)`` specifically for each gene.  If None,
   the values are calculated automatically from the *annodb*
   passed to the constructor.

