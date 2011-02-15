
===========================
Basic Phenoseq User Guide
===========================

------------
Dependencies
------------

Phenoseq depends on pygr (https://code.google.com/p/pygr/) and BioPython (http://www.biopython.org/wiki/SeqIO) for processing and parsing, and SciPy (http://www.scipy.org/) for statistics.

-------------------
Processing Raw Data
-------------------

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
found in the file "SampleSheet.csv". Because the tags are designed with 
error correction in mind, if a reported tag differs from one of the tags 
in the sample sheet file by just one base, we assume that it was tagged 
with that tag.

The script demultiplex.py can process the raw files and separate the reads 
by tag::

    python demultiplex.py SampleSheet.csv data_directory

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

Analyzing
---------

SNP analysis is handled by mkpileup.py and analyze.py. The first script produces vcf files from the .sam files from the alignment::

    python mkpileup.py NC_000913.gbk [ACGT]*.vcf. 

These files can then be processed by analyze.py::

    python analyze.py NC_000913.gbk [ACGT]*.vcf


---------------------------
Basic Classes and Functions
---------------------------

If the default usage is not sufficent, basic access to the processed data is easy.

Loading Data
------------

Initially, the data must be loaded from the processed files. First, the annotated reference genome is needed to determine if mutations are synonymous and in coding regions::

	>>> from analyze import *
	>>> annotated_genome_filename = "NC_000913.gbk"
	>>> annodb, al, dna = read_genome_annots(annotated_genome_filename)

This might take a couple of minutes on modest hardware. Next, read in the data from the *.vcf files. In Bash, the argument ``[ACGT]*.vcf`` (from above) expands to something like::

	user@home:~/phenoseq$ echo [ACGT]*.vcf
	ACAGTG.vcf ACTTGA.vcf ATCACG.vcf CAGATC.vcf CGATGT.vcf CTTGTA.vcf GATCAG.vcf GCCAAT.vcf TGACCA.vcf TTAGGC.vcf 

depending on your tags. In python, use::

	>>> tag_files = ['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']                                                
	>>> tag_dict = read_tag_files(tag_files)

or obtain the list from ``sys.argv`` by passing in the parameters shown above for analyze.py::

	user@home:~/phenoseq$ python2.5 - NC_000913.gbk [ACGT]*.vcf
	Python 2.5.5 (r255:77872, Sep 14 2010, 16:22:46) 
	[GCC 4.4.5] on linux2
	Type "help", "copyright", "credits" or "license" for more information.
	>>> import sys
	>>> sys.argv
	['-', 'NC_000913.gbk', 'ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']
	>>> annotated_genome_filename = sys.argv[1]
	>>> tag_files = sys.argv[2:]
	>>> annotated_genome_filename
	'NC_000913.gbk'
	>>> tag_files
	['ACAGTG.vcf', 'ACTTGA.vcf', 'ATCACG.vcf', 'CAGATC.vcf', 'CGATGT.vcf', 'CTTGTA.vcf', 'GATCAG.vcf', 'GCCAAT.vcf', 'TGACCA.vcf', 'TTAGGC.vcf']


Analyzing Data
--------------

Use the parsed data to initialize a ``GeneSNPDict``. 

	>>> gsd = GeneSNPDict(tag_dict, annodb, al, dna)

GeneSNPDict has an important boolean parameter ``count_syn`` (default is False) that tells the class to discard synonymous mutations.

An instance of GeneSNPDict has as keys coding regions in which mutations were detected, and might have a fairly large number of keys::

	>>> gsd.keys()
	['yfjK', 'speF', 'speG', 'speD', 'ybcL', 'ybcM', 'ybcO', 'ybcJ', 'ybcK', 'ybcF', 'nrdR', 'ybcY', 'pgpA', 'ypeA'...
	>>> len(gsd.keys())
	1428

The values of GeneSNPDict are lists of tuples consisting of a tag and a ``SNP`` object::

	>>> key = gsd.keys()[0]
	>>> gsd[key]
	[('GCCAAT', <analyze.SNP object at 0xb2aec0c>), ('TGACCA', <analyze.SNP object at 0xb936f6c>)]

indicating the library that the SNP was found in (via the tag). The associated SNP object provides information about the mutation. For instance, the following attributes indicate the position, quality, original base, new base, and the original amino acid and new amino acid of the codon::

	>>> snp = gsd[key][0][1]
	>>> print snp.pos, snp.qual, snp.ref, snp.alt, snp.aaRef, snp.aaAlt 
	2760136 99 G A P L

---------------------
Scoring Mutations
---------------------

GeneSNPDict has a method ``get_scores`` for scoring mutations. Calling this function returns a list of tuples of p-values and genes in ascending order::

	>>> scores = gsd.get_scores()[:20]
	>>> for score in scores:
	...     print score
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

