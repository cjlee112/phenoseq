
===========================
Basic Phenoseq User Guide
===========================

Obtaining Phenoseq
------------------

* You can either download the latest release of Phenoseq from its
  GitHub `downloads <https://github.com/cjlee112/phenoseq/tags>`_ page;

* or you can get the very latest code commits using Git, directly
  from the Phenoseq's `GitHub <https://github.com/cjlee112/phenoseq>`_ page.
  (This requires using `Git <http://git-scm.com>`_).

Installation
------------

Phenoseq is a pure Python package; you need 
`Python <http://python.org>`_ to run it.  You install Phenoseq
by running its ``setup.py`` install script::

  sudo python setup.py install

(leave out ``sudo`` if you are installing to a location you have
write privileges for).

Notes:

* it installs both the Phenoseq module and several
  convenient command-line scripts for running standard
  analyses (for details, see the examples below).
* if you have Python
  `setuptools <http://pypi.python.org/pypi/setuptools>`_,
  and `pip <http://pypi.python.org/pypi/pip>`_
  installed, the Phenoseq installer should be able to 
  automatically install its dependencies for you (see below
  for details).  Note that this can take several minutes!

Dependencies
............

(Automatically installed, if possible.  If you need to install
them manually, use your system's software installer tool or
`pip <http://pypi.python.org/pypi/pip>`_.)

* `numpy <http://numpy.scipy.org/>`_ and
  `scipy <http://www.scipy.org/>`_
  are used both for scoring experimental
  results, and for simulation.
  
Other dependencies are optional, depending on what phenoseq
functions you wish to use:

* Biopython's `SeqIO <http://www.biopython.org/wiki/SeqIO>`_
  module is used for reading
  Genbank CDS annotations for a genome, in
  :func:`analyze.read_genbank_annots()`.
  
* `Pygr <https://code.google.com/p/pygr/>`_
  is used to represent gene annotation intervals
  on a genome, for scoring experimental results.


Analyzing a Phenotype Sequencing Experiment
-------------------------------------------

The following command line scripts are installed automatically
when you install ``phenoseq`` (i.e. via ``python setup.py install``).

Variant call data (in the form of VCF files) 
can be processed by phenoseq in various ways. For single gene results, 
run the ``phenoseq_analyze`` command::

    phenoseq_analyze -g NC_000913.gbk [ACGT]*.vcf

It takes the following command options:

* ``-g GBFILEPATH``: required genbank file containing CDS annotations.  Also
  it uses this filename to look for a FASTA file containing the
  reference genome (if a FASTA file path is not explicitly provided
  by the ``--fastafile`` option).  Specifically, it will replace the
  suffix of the **GBFILEPATH** with the suffix ``.fna``.
* ``--fastafile FASTAFILEPATH``: (optional) path to FASTA file containing
  the reference genome sequence.
* ``--gene-qualifier GENEFIELD``: (optional) field name to use for
  extracting gene IDs from Genbank annotations.  Its default value is "gene";
  however, some Genbank files lack "gene" labels.  Other qualifiers that
  may be present if "gene" is missing are: "locus_tag"; or "protein_id".

Note that ``phenoseq_analyze`` extracts gene annotations (specifically,
CDS annotations) from the supplied Genbank file, and reports their
phenoseq scores using their associated gene IDs.  If a CDS annotation
lacks a gene ID it will simply be reported as "unlabeled_CDS_#"
(where # is an integer); in that case you can use the 
``--gene-qualifier`` option (see above) to tell it how to look 
for gene IDs.

For pathway phenoseq analysis, additional 
`EcoCyc <http://ecocyc.org>`_ files are required. Specifically, the file 
"func-associations.col" is necessary to supply the functionally associated groups, and "genes.col" 
are required to map the gene names to those used by the genbank annotation file. Then phenoseq can 
be run on functionally associated groups by::

    phenoseq_pathways -g NC_000913.gbk -f func-associations.col -t genes.col [ACGT]*.vcf

Similarly, Ka/Ks ratios and p-values can be computed as follows::

    phenoseq_kaks -g NC_000913.gbk -f func-associations.col -t genes.col [ACGT]*.vcff

and hypergeometric p-values with an additional parameter for the number of top genes to 
test for enrichment in::

    phenoseq_hypergeom -n 30 -g NC_000913.gbk -f func-associations.col -t genes.col [ACGT]*.vcf


Simulations
------------

For experimental design purposes, simulations are available to determine various expected yields and costs.
Invoke the simulation by::

    phenoseq_cost num_strains mean_coverage_per_lane cost_per_lane cost_per_library num_targets

For example,::

    user@home$ phenoseq_cost 20 30 1000 50 5
    Strains: 20
    Coverage per lane: 30
    Cost per lane: $1000
    Cost per library: $50
    Target genes: 5

    Minimum experiment cost: $7166.67
    Total libraries: 10
    Coverage per library: 20
    Expected hits: 3.71


Processing Raw Data
-------------------

NOTE: the following scripts are provided only as examples,
since most groups have their own preferred pipelines for processing
nextgen sequencing data and generating variant calls.
The scripts shown in these examples are located in the 
``phenoseq/examples`` directory, and were designed to work
with Illumina tagged single-end read data.

To process raw sequencer data, we must complete the following steps

* Demixing
* Aligning
* Analyzing

Demixing
--------

The first task is to demultiplex the raw pooled sequence fragments. 
Let us assume the data filenames have the Illumina Native format ``s_L_{1,2}_n_qseq.txt``. 
The number L is the lane, and the number n is the tile number. 
For each tile, there are 2 files. The fragments are in the file with a 
1 after the second underscore and the tags associated to those fragments 
are in the file with the 2 after the underscore. The tags used can be 
found in the file "SampleSheet.csv" and should be given directly in the command line since the file may noto be available. Because the tags are designed with 
error correction in mind, if a reported tag differs from one of the tags 
in the sample sheet file by just one base, we assume that it was tagged 
with that tag.

The script demultiplex.py can process the raw files and separate the reads 
by tag::

    python demultiplex.py data_directory AACTCG ATGTGC GTCATT ...

The output is a collection of files of the form ``s_(lane)_(tag)_qseq.txt`` .

Aligning
--------

The resulting demixed reads can be aligned with novoalign. First we must build an indexed version of a refence genome. If the reference genome file is named ``NC_000913.fna``, the command to build the index is::

    ./novocraft/novoindex ecoli.nix NC_000913.fna

which creates the index ``ecoli.nix``. 

Alignment then proceeds with the command::

    ./novoalign/novocraft/novoalign -d ecoli.nix -f s_1_ACTTGA_qseq.txt

which aligns the reads in the file ``s_1_ACTTGA_qseq.txt`` against the index. Since there will be several such files, we can automate this process with a simple script align_batch.py::

    python align_batch.py <novoalign_executable> <index_filename>

e.g.

    python align_batch.py novoalign/novocraft/novoalign ssuis.nix

This script executes commands such as

    ./novoalign/novocraft/novoalign -d ssuis.nix -o SAM -f s_1_CGATGT_qseq.txt > aligned_s_1_CGATGT.sam

The ``-o SAM`` option outputs the data in SAMTOOLS format and the aligned reads are in aligned_s_1_CGATGT.sam .

See 

http://www.novocraft.com/wiki/tiki-index.php?page=Getting+Started&structure=Novocraft+Technologies&page_ref_id=70

for more information on using novoalign.




Custom Analysis in the Python Interpreter
-----------------------------------------

If the default usage is not sufficent, basic access to the processed data is easy.

Loading Data
............

Initially, the data must be loaded from the processed files. First, the annotated reference genome is needed to determine if mutations are synonymous and in coding regions::

	>>> from phenoseq.analyze import *
	>>> annotated_genome_filename = "NC_000913.gbk"
	>>> annodb, al, genome = read_genbank_annots(annotated_genome_filename)

This might take a couple of minutes on modest hardware. Next, read in the data from the VCF files.  In python, use::

	>>> import glob
	>>> tag_files = glob.glob('*.vcf')
        >>> tag_files
        ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']                                                
	>>> snps = read_tag_files(tag_files)


Analyzing Data
..............

The result from any SNP reading function such as :func:`analyze.read_vcf`
or :func:`analyze.read_tag_files` is a list of :class:`analyze.SNP` objects.
We can inspect the first few::

        >>> snps[:5]
        [<SNP chr1:7682394:G:C>, <SNP chr1:23847535:A:G>, ...]

The next step in analysis is to map the SNPs to genes, using the alignment
object obtained above, which maps sequence intervals to gene CDS intervals.  
Here's a simple example that assumes all the SNPs map on one DNA sequence 
(e.g. a microbial genome)::

        >>> gsd = map_snps(snps, al, genome)

The result is a gene:snp dictionary, whose keys are gene IDs,
and whose values are lists of SNPs found in that gene::

        >>> gsd
        {'fugA':[<SNP chr1:343652:T:C>], ...}

We can filter these results to just nonsynonymous SNPs::

        >>> gsd = filter_nonsyn(gsd)



Scoring Mutations
.................

Finally, we score the genes for significant p-values::

        >>> scores = score_genes_pooled(gsd, genome=genome, annodb=annodb)
	>>> for hit in scores:
	...     print hit
	... 
	(6.7585463507686869e-23, 'acrB')
	(9.6429750487530072e-09, 'marC')
	(1.2481477231487551e-07, 'stfP')
	(7.6301063544178727e-07, 'ykgC')
	(2.4971914594342781e-06, 'aes')
	(1.2133651191132762e-05, 'ampH')
	(2.6930241003283795e-05, 'paoC')
	(2.7593050733850882e-05, 'nfrA')
	(3.0833069533854329e-05, 'ydhB')
	(8.2645380374133238e-05, 'yaiP')
	(0.00011995056941060593, 'acrA')
	(0.00017251088960147507, 'xanQ')
	(0.0001786206550615194, 'ykgD')
	(0.0002480120870963014, 'yegQ')
	(0.00024916389158152248, 'yfjJ')
	(0.00026148314689727225, 'yagX')
	(0.00032324465826595041, 'pstA')
	(0.0003368649972321227, 'prpE')
	(0.00035174665129372739, 'mltF')
	(0.00044489155029703195, 'purE')

